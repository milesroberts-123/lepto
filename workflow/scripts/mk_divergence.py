#!/usr/bin/env python3
"""Count fixed synonymous/nonsynonymous differences between 1:1 orthologs.

Walks each codon alignment from mk_align_pairs.py codon by codon. Codons
with a gap or stop in either sequence are skipped; codons whose ingroup
positions include a site that is polymorphic in any pool are excluded
(shared polymorphism is not divergence, following the DnaSP/MKado
convention). Remaining differing codons are scored with
classify_codon_change (minimum-replacement paths averaged), giving
fractional Dn/Ds for multi-change codons.

Memory design: the CDS map is streamed once, transcript block by
transcript block (mk_extract_cds.py writes it transcript-sorted), and
each block is a few thousand rows; only the block whose transcript has
an alignment is kept. The polymorphic-site set is held as per-scaffold
sorted `array('q')` vectors with bisect lookups rather than a set of
(chrom, pos) tuples. Peak RSS is therefore bounded by the polymorphism
table alone (~8 B/site), independent of the genome-wide CDS map size —
the first version loaded every CDS position into a dict of tuples and
was OOM-killed on the real annotation (412 M positions).

Output TSV: orthogroup gene dn ds n_diff_codons n_excluded_poly
            n_excluded_gap n_stop
"""

import argparse
import bisect
import os
import sys
from array import array

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mk_common import (classify_codon_change, gene_of, read_fasta,  # noqa: E402
                       split_alignment, translate)


def load_polymorphic(path):
    """polymorphic_sites.tsv -> {chrom: sorted array('q') of positions}."""
    acc = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.split()
            if len(f) < 2:
                continue
            acc.setdefault(f[0], array("q")).append(int(f[1]))
    for chrom, arr in acc.items():
        acc[chrom] = array("q", sorted(set(arr)))
    return acc


def is_poly(poly, chrom, pos):
    arr = poly.get(chrom)
    if not arr:
        return False
    i = bisect.bisect_left(arr, pos)
    return i < len(arr) and arr[i] == pos


def iter_map_blocks(path):
    """Yield (transcript, [(cdspos, chrom, pos), ...]) in file order."""
    with open(path) as fh:
        fh.readline()
        cur = None
        block = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            tx = f[0]
            if tx != cur:
                if cur is not None:
                    yield cur, block
                cur = tx
                block = []
            block.append((int(f[2]), f[3], int(f[4])))
        if cur is not None:
            yield cur, block


def alignment_headers(path, ingroup):
    """First record's name if this file looks like an ingroup alignment."""
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                if not name.startswith(ingroup + "|"):
                    return None
                return name.split("|", 1)[1]
    return None


def count_gene(ing, out, positions, poly):
    """Return (dn, ds, n_diff, n_poly, n_gap, n_stop) for one alignment."""
    dn = ds = 0.0
    n_diff = n_poly = n_gap = n_stop = 0
    in_idx = 0  # ingroup bases consumed so far (1-based cdspos)
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
        if translate(a_codon) in ("*", "X") or translate(b_codon) in ("*", "X"):
            in_idx += 3
            n_stop += 1
            continue
        if a_codon != b_codon:
            diff_positions = [k for k in range(3) if a_codon[k] != b_codon[k]]
            genomic = []
            for k in diff_positions:
                hit = positions.get(in_idx + k + 1)
                if hit is None:
                    genomic = None
                    break
                genomic.append(hit)
            if genomic is None:
                n_gap += 1
            elif any(is_poly(poly, chrom, pos) for chrom, pos in genomic):
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

    poly = load_polymorphic(args.polymorphic_sites)
    print(f"polymorphic sites: "
          f"{sum(len(a) for a in poly.values())} in {len(poly)} scaffolds")

    # pre-pass: transcript -> alignment file (only cheap header reads)
    wanted = {}
    failures = []
    for name in sorted(os.listdir(args.alignments)):
        if not name.endswith(".fna"):
            continue
        og = name[:-4]
        path = os.path.join(args.alignments, name)
        tx = alignment_headers(path, args.ingroup)
        if tx is None:
            failures.append((og, "alignment_names"))
            continue
        wanted[tx] = (og, path)
    print(f"{len(wanted)} alignment files, {len(failures)} unusable")

    rows = []
    seen = set()
    for tx, block in iter_map_blocks(args.map):
        if tx not in wanted:
            continue
        seen.add(tx)
        og, path = wanted[tx]
        positions = {cdspos: (chrom, pos) for cdspos, chrom, pos in block}
        seqs = read_fasta(path)
        split = split_alignment(seqs, args.ingroup, args.outgroup)
        if split is None:
            failures.append((og, "alignment_names"))
            continue
        ing, out, _ = split
        dn, ds, n_diff, n_poly, n_gap, n_stop = count_gene(
            ing, out, positions, poly)
        rows.append((og, gene_of(tx), dn, ds, n_diff, n_poly, n_gap, n_stop))

    for tx, (og, path) in wanted.items():
        if tx not in seen:
            failures.append((og, "no_cds_map"))

    rows.sort()
    failures.sort()
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
