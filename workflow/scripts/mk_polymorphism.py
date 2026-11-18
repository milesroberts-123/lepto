#!/usr/bin/env python3
"""Count Pn/Ps per gene and per pool from pool allele frequencies.

Classification is by translation, not by degeneracy class: for each
polymorphic CDS site the reference codon is read from the CDS fasta,
the alternative codon is that codon with the observed ALT base at the
site's codon position, and the site counts as synonymous (Ps) when both
codons encode the same amino acid, nonsynonymous (Pn) otherwise
(stop-gain counts as Pn). 2-fold sites therefore count too, classified
by the allele actually observed, and no degeneracy annotation is needed.

Inputs:
  --map   CDS map from mk_extract_cds.py (transcript gene cdspos chrom pos)
  --cds   CDS fasta from mk_extract_cds.py (transcript -> CDS)
  --freq  ID=path of a per-pool allele-frequency table (repeatable)

The frequency table may be the raw grenedalf CSV (header + 8 columns:
CHROM,POS,REF,ALT,...,ALT_CNT,DEPTH,FREQ -- the ldx.smk/freqk.smk
positional contract) or the 7-column extract written by mk.smk's
mk_frequency_filter (chrom,pos,ref,alt,alt_count,depth,freq).

A site is polymorphic in a pool when depth >= --min-depth, ALT count >=
--min-alt-count and 0 < ALT frequency < 1. The pooled table counts a
site once if it is polymorphic in any pool. **Multi-allelic sites** —
the same position carrying different ALT alleles in different pools —
are dropped from the Pn/Ps counts (the synonymous/nonsynonymous call
would depend on which allele you pick) but are kept in
polymorphic_sites.tsv so the divergence step still excludes that codon.

Memory design: only polymorphic sites are indexed. The union of
per-pool polymorphic positions is held per scaffold as numpy int64
position vectors (8 B/site) plus a 1 B allele code; the map is streamed
once in transcript blocks to classify those sites; a second pass over
the pool tables accumulates counts. Memory scales with the number of
polymorphic sites, not with the genome-wide CDS site count (85 M at
this annotation).

Outputs (in --outdir):
  <ID>_polymorphism.tsv   gene pn ps
  pooled_polymorphism.tsv gene pn ps
  polymorphic_sites.tsv   chrom pos  (union of pools, any CDS site)
"""

import argparse
import os
import sys
from array import array

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mk_common import gene_of, read_fasta, translate  # noqa: E402

ALT_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
CODE_ALT = {v: k for k, v in ALT_CODE.items()}
UNCLASSIFIED = 255
INVALID = 254


def parse_record(fields):
    """(chrom, pos, ref, alt, alt_count, depth, freq) or None."""
    if not fields or fields[0] == "CHROM":
        return None
    try:
        if len(fields) == 7:      # mk_frequency_filter extract
            return (fields[0], int(fields[1]), fields[2].upper(),
                    fields[3].upper(), int(fields[4]), int(fields[5]),
                    float(fields[6]))
        if len(fields) >= 8:      # raw grenedalf frequency CSV
            return (fields[0], int(fields[1]), fields[2].upper(),
                    fields[3].upper(), int(fields[5]), int(fields[6]),
                    float(fields[7]))
    except ValueError:
        return None
    return None


def iter_records(path, min_depth, min_alt):
    with open(path) as fh:
        for line in fh:
            rec = parse_record(line.rstrip("\n").split(","))
            if rec is None:
                continue
            chrom, pos, ref, alt, ac, depth, freq = rec
            if depth < min_depth or ac < min_alt:
                continue
            if not (0.0 < freq < 1.0):
                continue
            yield chrom, pos, ref, alt


def iter_map_blocks(path):
    """Yield (transcript, [(chrom, cdspos, pos, gene), ...]) in file order."""
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
            block.append((f[3], int(f[2]), int(f[4]), f[1]))
        if cur is not None:
            yield cur, block


def classify(codon, alt, codon_pos):
    """0 = synonymous, 1 = nonsynonymous, INVALID if unclassifiable.

    The REF base is the codon base itself at this CDS position (the
    frequency table's REF column is cross-checked separately), so the
    mutant codon is unambiguous.
    """
    if len(codon) != 3 or any(b not in ALT_CODE for b in codon):
        return INVALID
    if alt not in ALT_CODE:
        return INVALID
    mutant = codon[:codon_pos] + alt + codon[codon_pos + 1:]
    ref_aa = translate(codon)
    alt_aa = translate(mutant)
    if ref_aa == "X" or alt_aa == "X":
        return INVALID
    return 0 if ref_aa == alt_aa else 1


class SiteIndex:
    """Per-scaffold union of polymorphic positions plus classification."""

    def __init__(self):
        self._acc = {}            # chrom -> [array('q'), bytearray()]

    def add(self, chrom, pos, alt):
        entry = self._acc.get(chrom)
        if entry is None:
            entry = self._acc[chrom] = [array("q"), bytearray()]
        entry[0].append(pos)
        entry[1].append(ALT_CODE.get(alt, 255))

    def finalize(self):
        self.pos = {}
        self.alt = {}
        self.gene = {}
        self.cls = {}
        self.n_multiallelic = 0
        for chrom, (pos_arr, alt_arr) in self._acc.items():
            pos = np.frombuffer(pos_arr, dtype=np.int64).copy()
            alt = np.frombuffer(bytes(alt_arr), dtype=np.uint8).copy()
            if len(pos) > 1:
                order = np.lexsort((alt, pos))   # by position, then allele
                pos = pos[order]
                alt = alt[order]
                starts = np.flatnonzero(
                    np.concatenate(([True], pos[1:] != pos[:-1])))
                ends = np.concatenate((starts[1:], [len(pos)]))
                first_alt = alt[starts]
                last_alt = alt[ends - 1]
                multi = first_alt != last_alt
                pos = pos[starts]
                alt = first_alt
                self.n_multiallelic += int(multi.sum())
            else:
                multi = np.zeros(len(pos), dtype=bool)
            self.pos[chrom] = pos
            self.alt[chrom] = alt
            self.gene[chrom] = np.full(len(pos), -1, dtype=np.int32)
            cls = np.full(len(pos), UNCLASSIFIED, dtype=np.uint8)
            cls[multi] = INVALID
            self.cls[chrom] = cls
        self._acc = {}

    def lookup(self, chrom, pos):
        arr = self.pos.get(chrom)
        if arr is None or len(arr) == 0:
            return None, -1
        i = int(np.searchsorted(arr, pos))
        if i >= len(arr) or arr[i] != pos:
            return None, -1
        return arr, i


def write_counts(path, counts, genes):
    with open(path, "w") as out:
        out.write("gene\tpn\tps\n")
        for gid in sorted(counts, key=genes.__getitem__):
            pn, ps = counts[gid]
            out.write(f"{genes[gid]}\t{pn}\t{ps}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--freq", action="append", required=True,
                    help="ID=path of a per-pool allele-frequency table")
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--min-alt-count", type=int, default=3)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    pools = []
    for spec in args.freq:
        pool, sep, path = spec.partition("=")
        if not sep or not pool or not path:
            sys.exit(f"bad --freq value: {spec!r} (expected ID=path)")
        pools.append((pool, path))

    cds = read_fasta(args.cds)
    print(f"{len(cds)} CDS sequences loaded")

    # pass A: union of polymorphic sites across pools
    index = SiteIndex()
    for pool, path in pools:
        n = 0
        for chrom, pos, ref, alt in iter_records(
                path, args.min_depth, args.min_alt_count):
            index.add(chrom, pos, alt)
            n += 1
        print(f"pool {pool}: {n} polymorphic records")
    index.finalize()
    n_union = sum(len(v) for v in index.pos.values())
    print(f"{n_union} polymorphic sites in the union "
          f"({index.n_multiallelic} multi-allelic, dropped from Pn/Ps)")

    # pass B: classify union sites from the streamed CDS map
    genes = []
    gene_ids = {}
    n_matched = 0
    for tx, block in iter_map_blocks(args.map):
        seq = cds.get(tx)
        if seq is None:
            continue
        gid = gene_ids.get(tx)
        if gid is None:
            gid = len(genes)
            gene_ids[tx] = gid
            genes.append(gene_of(tx))
        for chrom, cdspos, pos, _gene in block:
            arr, i = index.lookup(chrom, pos)
            if arr is None or index.cls[chrom][i] != UNCLASSIFIED:
                continue
            start = ((cdspos - 1) // 3) * 3
            codon = seq[start:start + 3].upper()
            alt = CODE_ALT.get(int(index.alt[chrom][i]), "")
            index.cls[chrom][i] = classify(codon, alt, (cdspos - 1) % 3)
            index.gene[chrom][i] = gid
            n_matched += 1
    print(f"{n_matched} sites matched to a CDS")

    n_classified = 0
    for chrom, cls in index.cls.items():
        cls[cls == UNCLASSIFIED] = INVALID
        n_classified += int(((cls == 0) | (cls == 1)).sum())
    print(f"{n_classified} of {n_union} sites classifiable")

    os.makedirs(args.outdir, exist_ok=True)

    # pass C: per-pool counts (multi-allelic sites are already INVALID)
    for pool, path in pools:
        counts = {}
        n = 0
        for chrom, pos, ref, alt in iter_records(
                path, args.min_depth, args.min_alt_count):
            arr, i = index.lookup(chrom, pos)
            if arr is None:
                continue
            cls = int(index.cls[chrom][i])
            if cls > 1:
                continue
            gid = int(index.gene[chrom][i])
            if gid < 0:
                continue
            entry = counts.get(gid)
            if entry is None:
                entry = counts[gid] = [0, 0]
            # counts are [Pn, Ps]; cls 0 = synonymous -> Ps, 1 = nonsyn -> Pn
            entry[0 if cls == 1 else 1] += 1
            n += 1
        write_counts(os.path.join(args.outdir, f"{pool}_polymorphism.tsv"),
                     counts, genes)
        print(f"pool {pool}: {n} sites counted in {len(counts)} genes")

    # pooled table: every union site counts once
    pooled = {}
    for chrom, cls in index.cls.items():
        gids = index.gene[chrom]
        for code in (0, 1):
            sel = np.flatnonzero(cls == code)
            if len(sel) == 0:
                continue
            sel = sel[gids[sel] >= 0]
            uniq, cnt = np.unique(gids[sel], return_counts=True)
            slot = 0 if code == 1 else 1     # 1 = nonsyn -> Pn
            for gid, c in zip(uniq.tolist(), cnt.tolist()):
                entry = pooled.get(gid)
                if entry is None:
                    entry = pooled[gid] = [0, 0]
                entry[slot] += c
    write_counts(os.path.join(args.outdir, "pooled_polymorphism.tsv"),
                 pooled, genes)

    with open(os.path.join(args.outdir, "polymorphic_sites.tsv"), "w") as out:
        out.write("chrom\tpos\n")
        for chrom in sorted(index.pos):
            arr = index.pos[chrom]
            if len(arr) == 0:
                continue
            for start in range(0, len(arr), 200000):
                chunk = arr[start:start + 200000].tolist()
                out.write("".join(f"{chrom}\t{p}\n" for p in chunk))
    print(f"{n_union} polymorphic CDS sites (union of pools)")


if __name__ == "__main__":
    main()
