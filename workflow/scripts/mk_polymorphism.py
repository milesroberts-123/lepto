#!/usr/bin/env python3
"""Count Pn/Ps per gene and per pool from grenedalf allele frequencies.

Inputs:
  --map    hap1 CDS map from mk_extract_cds.py (transcript gene cdspos chrom pos)
  --zero   degenotate 0-fold BED (nonsynonymous sites, prefixed coords)
  --four   degenotate 4-fold BED (synonymous sites, prefixed coords)
  --freq   ID=path of a grenedalf frequency CSV (repeatable)

grenedalf frequency CSV contract (per-sample run with
--write-sample-counts/--write-sample-read-depth/--write-sample-alt-freq,
see vg_diversity.smk grenedalf_frequency): col1 CHROM, col2 POS (1-based),
col6 ALT count, col7 read depth, col8 ALT frequency. This is the same
positional contract used by ldx.smk and freqk.smk.

A site is polymorphic in a pool when depth >= --min-depth,
ALT count >= --min-alt-count and 0 < ALT frequency < 1. Sites are
classified by degenotate degeneracy (0-fold -> Pn, 4-fold -> Ps; 2-fold
and non-CDS sites are dropped). The pooled table counts a site once if
it is polymorphic in any pool.

Outputs (in --outdir):
  <ID>_polymorphism.tsv   gene pn ps
  pooled_polymorphism.tsv gene pn ps
  polymorphic_sites.tsv   chrom pos  (union, for divergence exclusion)
"""

import argparse
import os
import sys


def read_bed_sites(path):
    sites = set()
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                sites.add(f"{f[0]}:{int(f[1]) + 1}")  # BED 0-based -> 1-based pos
    return sites


def build_site_index(map_path, zero, four):
    # key "chrom:pos" (1-based) -> (gene, 'n'|'s'), 0/4-fold sites only
    index = {}
    cds_sites = set()      # every CDS site, any degeneracy (exclusion set)
    with open(map_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            gene, chrom, pos = f[1], f[3], f[4]
            key = f"{chrom}:{pos}"
            cds_sites.add(key)
            if key in index:
                continue
            if key in zero:
                index[key] = (gene, "n")
            elif key in four:
                index[key] = (gene, "s")
    return index, cds_sites


def count_pool(path, index, cds_sites, min_depth, min_alt):
    counts = {}          # gene -> [pn, ps] (0/4-fold sites only)
    poly_keys = set()    # any CDS site polymorphic in this pool
    n_poly = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split(",")
            if len(f) < 8:
                continue
            chrom = f[0]
            if chrom == "CHROM":
                continue
            try:
                pos = int(f[1])
                alt = int(f[5])
                depth = int(f[6])
                freq = float(f[7])
            except ValueError:
                continue
            if depth < min_depth or alt < min_alt:
                continue
            if not (0.0 < freq < 1.0):
                continue
            key = f"{chrom}:{pos}"
            if key not in cds_sites:
                continue
            poly_keys.add(key)
            n_poly += 1
            hit = index.get(key)
            if hit is None:
                continue
            gene, fold = hit
            entry = counts.setdefault(gene, [0, 0])
            entry[0 if fold == "n" else 1] += 1
    return counts, poly_keys, n_poly


def write_counts(path, counts):
    with open(path, "w") as out:
        out.write("gene\tpn\tps\n")
        for gene in sorted(counts):
            pn, ps = counts[gene]
            out.write(f"{gene}\t{pn}\t{ps}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map", required=True)
    ap.add_argument("--zero", required=True)
    ap.add_argument("--four", required=True)
    ap.add_argument("--freq", action="append", required=True,
                    help="ID=path of a grenedalf frequency CSV")
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--min-alt-count", type=int, default=1)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    zero = read_bed_sites(args.zero)
    four = read_bed_sites(args.four)
    index, cds_sites = build_site_index(args.map, zero, four)
    if not index:
        sys.exit("no CDS sites overlap the degeneracy BEDs")
    print(f"{len(index)} 0/4-fold CDS sites indexed "
          f"({len(cds_sites)} CDS sites total)")

    os.makedirs(args.outdir, exist_ok=True)
    pooled_keys = set()

    for spec in args.freq:
        pool, sep, path = spec.partition("=")
        if not sep or not pool or not path:
            sys.exit(f"bad --freq value: {spec!r} (expected ID=path)")
        counts, keys, n_poly = count_pool(
            path, index, cds_sites, args.min_depth, args.min_alt_count)
        write_counts(os.path.join(args.outdir, f"{pool}_polymorphism.tsv"), counts)
        pooled_keys |= keys
        print(f"pool {pool}: {n_poly} polymorphic CDS sites in "
              f"{len(counts)} genes")

    # pooled table: every site counts once, even if polymorphic in several
    # pools; only 0/4-fold sites enter Pn/Ps
    pooled_counts = {}
    for key in pooled_keys:
        hit = index.get(key)
        if hit is None:
            continue
        gene, fold = hit
        entry = pooled_counts.setdefault(gene, [0, 0])
        entry[0 if fold == "n" else 1] += 1
    write_counts(os.path.join(args.outdir, "pooled_polymorphism.tsv"),
                 pooled_counts)

    with open(os.path.join(args.outdir, "polymorphic_sites.tsv"), "w") as out:
        out.write("chrom\tpos\n")
        for key in sorted(pooled_keys, key=lambda k: (k.split(":")[0],
                                                      int(k.split(":")[1]))):
            chrom, pos = key.split(":")
            out.write(f"{chrom}\t{pos}\n")
    print(f"{len(pooled_keys)} polymorphic CDS sites (union of pools)")


if __name__ == "__main__":
    main()
