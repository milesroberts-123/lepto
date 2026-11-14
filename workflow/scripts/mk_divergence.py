#!/usr/bin/env python3
"""Count fixed synonymous/nonsynonymous differences between 1:1 orthologs.

Walks each codon alignment from mk_align_pairs.py codon by codon. Codons
with a gap or stop in either sequence are skipped; codons whose ingroup
positions include a site that is polymorphic in any pool are excluded
(shared polymorphism is not divergence, following the DnaSP/MKado
convention). Remaining differing codons are scored with
classify_codon_change (minimum-replacement paths averaged), giving
fractional Dn/Ds for multi-change codons.

Output TSV: orthogroup gene dn ds n_diff_codons n_excluded_poly
            n_excluded_gap n_stop
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mk_common import (classify_codon_change, gene_of, read_fasta,  # noqa: E402
                       split_alignment, translate)


def load_map(path):
    positions = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            tx, cdspos, chrom, pos = f[0], int(f[2]), f[3], int(f[4])
            positions.setdefault(tx, []).append((cdspos, chrom, pos))
    for tx in positions:
        positions[tx].sort()
    return positions


def count_gene(ing, out, tx, map_positions, poly):
    """Return (dn, ds, n_diff, n_poly, n_gap, n_stop) for one alignment."""
    pos_list = map_positions.get(tx)
    if pos_list is None:
        return None
    dn = ds = 0.0
    n_diff = n_poly = n_gap = n_stop = 0
    in_idx = 0  # number of ingroup bases consumed so far (1-based cdspos)
    length = len(ing)
    for start in range(0, length - 2, 3):
        a_codon = ing[start:start + 3]
        b_codon = out[start:start + 3]
        if "-" in a_codon:
            in_idx += sum(1 for ch in a_codon if ch != "-")
            n_gap += 1
            continue
        if "-" in b_codon or "N" in a_codon or "N" in b_codon:
            in_idx += 3
            n_gap += 1
            continue
        if translate(a_codon) == "*" or translate(b_codon) == "*" \
                or translate(a_codon) == "X" or translate(b_codon) == "X":
            in_idx += 3
            n_stop += 1
            continue
        if a_codon != b_codon:
            diff_positions = [k for k in range(3) if a_codon[k] != b_codon[k]]
            genomic = []
            for k in diff_positions:
                cdspos = in_idx + k + 1
                if cdspos > len(pos_list):
                    genomic = None
                    break
                genomic.append(pos_list[cdspos - 1])
            if genomic is None:
                n_gap += 1
            elif any((g[1], g[2]) in poly for g in genomic):
                n_poly += 1
            else:
                syn, nonsyn = classify_codon_change(a_codon, b_codon)
                ds += syn
                dn += nonsyn
                n_diff += 1
        in_idx += 3
    return dn, ds, n_diff, n_poly, n_gap, n_stop


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignments", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--polymorphic-sites", required=True)
    ap.add_argument("--ingroup", required=True)
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--failures", required=True)
    args = ap.parse_args()

    map_positions = load_map(args.map)
    poly = set()
    with open(args.polymorphic_sites) as fh:
        fh.readline()
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                poly.add((f[0], int(f[1])))

    rows = []
    failures = []
    for name in sorted(os.listdir(args.alignments)):
        if not name.endswith(".fna"):
            continue
        og = name[:-4]
        seqs = read_fasta(os.path.join(args.alignments, name))
        split = split_alignment(seqs, args.ingroup, args.outgroup)
        if split is None:
            failures.append((og, "alignment_names"))
            continue
        ing, out, tx = split
        counted = count_gene(ing, out, tx, map_positions, poly)
        if counted is None:
            failures.append((og, "no_cds_map"))
            continue
        dn, ds, n_diff, n_poly, n_gap, n_stop = counted
        rows.append((og, gene_of(tx), dn, ds, n_diff, n_poly, n_gap, n_stop))

    with open(args.out, "w") as out:
        out.write("orthogroup\tgene\tdn\tds\tn_diff_codons\t"
                  "n_excluded_poly\tn_excluded_gap\tn_stop\n")
        for og, gene, dn, ds, n_diff, n_poly, n_gap, n_stop in rows:
            out.write(f"{og}\t{gene}\t{dn:.4f}\t{ds:.4f}\t{n_diff}\t"
                      f"{n_poly}\t{n_gap}\t{n_stop}\n")

    with open(args.failures, "w") as out:
        out.write("orthogroup\treason\n")
        for og, reason in failures:
            out.write(f"{og}\t{reason}\n")

    print(f"{len(rows)} genes with divergence counts, "
          f"{len(failures)} skipped")


if __name__ == "__main__":
    main()
