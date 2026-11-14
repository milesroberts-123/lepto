#!/usr/bin/env python3
"""Count Pn/Ps per gene and per pool from grenedalf allele frequencies.

Inputs:
  --map    hap1 CDS map from mk_extract_cds.py (transcript gene cdspos chrom pos)
  --zero   degenotate 0-fold BED (nonsynonymous sites, prefixed coords)
  --four   degenotate 4-fold BED (synonymous sites, prefixed coords)
  --freq   ID=path of a per-pool allele-frequency table (repeatable)

The frequency table may be either the raw grenedalf CSV (header + 8
columns: CHROM, POS, ..., ALT count col6, depth col7, ALT frequency
col8 -- the ldx.smk/freqk.smk positional contract) or the compact
5-column polymorphic-only extract written by mk.smk's
mk_frequency_filter (chrom, pos, alt, depth, freq). Both are accepted;
the site thresholds are applied here either way, so the two paths give
identical counts.

A site is polymorphic in a pool when depth >= --min-depth,
ALT count >= --min-alt-count and 0 < ALT frequency < 1. Sites are
classified by degenotate degeneracy (0-fold -> Pn, 4-fold -> Ps; 2-fold
and non-CDS sites are dropped from the counts but any polymorphic CDS
site joins the exclusion set). The pooled table counts a site once if
it is polymorphic in any pool.

Memory design: the site index is held as per-scaffold sorted
`array('q')` position vectors with small parallel class/gene arrays and
bisect lookups, not as sets of "chrom:pos" strings. Peak RSS is roughly
20 bytes per indexed CDS site (~1-2 GB at genome scale) instead of the
several GB per string-keyed container the first version used, so the
step runs inside a couple of cores' worth of derived memory.

Outputs (in --outdir):
  <ID>_polymorphism.tsv   gene pn ps
  pooled_polymorphism.tsv gene pn ps
  polymorphic_sites.tsv   chrom pos  (union of all CDS sites, any class)
"""

import argparse
import bisect
import os
import sys
from array import array

POS_BITS = 40
POS_MASK = (1 << POS_BITS) - 1


def read_bed_arrays(path):
    """BED -> {chrom: sorted array('q') of 1-based positions}."""
    acc = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 3:
                continue
            acc.setdefault(f[0], array("q")).append(int(f[1]) + 1)
    for chrom, arr in acc.items():
        acc[chrom] = array("q", sorted(set(arr)))
    return acc


def contains(arr, pos):
    if not arr:
        return False
    i = bisect.bisect_left(arr, pos)
    return i < len(arr) and arr[i] == pos


def build_index(map_path, zero, four):
    """Split the CDS map into 0/4-fold and any-degeneracy per-scaffold arrays."""
    folded_pos, folded_cls, folded_gene = {}, {}, {}
    cds_pos = {}
    genes = []            # gene id -> name
    gene_ids = {}
    n_rows = 0
    with open(map_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos = f[3], int(f[4])
            cds_pos.setdefault(chrom, array("q")).append(pos)
            n_rows += 1
            if contains(zero.get(chrom), pos):
                cls = 0
            elif contains(four.get(chrom), pos):
                cls = 1
            else:
                continue
            gene = f[1]
            gid = gene_ids.get(gene)
            if gid is None:
                gid = len(genes)
                gene_ids[gene] = gid
                genes.append(gene)
            folded_pos.setdefault(chrom, array("q")).append(pos)
            folded_cls.setdefault(chrom, array("B")).append(cls)
            folded_gene.setdefault(chrom, array("i")).append(gid)

    for chrom, arr in cds_pos.items():
        cds_pos[chrom] = array("q", sorted(arr))
    # transcript order is not coordinate order: sort each scaffold's
    # parallel folded arrays by position (per-scaffold transient only)
    for chrom, arr in folded_pos.items():
        order = sorted(range(len(arr)), key=arr.__getitem__)
        folded_pos[chrom] = array("q", (arr[i] for i in order))
        folded_cls[chrom] = array("B", (folded_cls[chrom][i] for i in order))
        folded_gene[chrom] = array("i", (folded_gene[chrom][i] for i in order))
    return folded_pos, folded_cls, folded_gene, cds_pos, genes, n_rows


def parse_freq_line(fields):
    """Return (chrom, pos, alt, depth, freq) or None for header/bad rows."""
    if not fields or fields[0] == "CHROM":
        return None
    try:
        if len(fields) == 5:      # mk_frequency_filter extract
            return (fields[0], int(fields[1]), int(fields[2]),
                    int(fields[3]), float(fields[4]))
        if len(fields) >= 8:      # raw grenedalf frequency CSV
            return (fields[0], int(fields[1]), int(fields[5]),
                    int(fields[6]), float(fields[7]))
    except ValueError:
        return None
    return None


def count_pool(path, index, min_depth, min_alt, chrom_ids):
    folded_pos, folded_cls, folded_gene, cds_pos, _, _ = index
    counts = {}           # gene id -> [pn, ps]
    poly_keys = set()     # encoded (chrom id << 40 | pos), any CDS class
    n_poly = 0
    with open(path) as fh:
        for line in fh:
            rec = parse_freq_line(line.rstrip("\n").split(","))
            if rec is None:
                continue
            chrom, pos, alt, depth, freq = rec
            if depth < min_depth or alt < min_alt:
                continue
            if not (0.0 < freq < 1.0):
                continue
            if not contains(cds_pos.get(chrom), pos):
                continue
            cid = chrom_ids.get(chrom)
            if cid is None:
                cid = chrom_ids[chrom] = len(chrom_ids)
            key = (cid << POS_BITS) | pos
            if key in poly_keys:
                continue
            poly_keys.add(key)
            n_poly += 1
            arr = folded_pos.get(chrom)
            if not arr:
                continue
            i = bisect.bisect_left(arr, pos)
            if i >= len(arr) or arr[i] != pos:
                continue
            gid = folded_gene[chrom][i]
            entry = counts.get(gid)
            if entry is None:
                entry = counts[gid] = [0, 0]
            entry[0 if folded_cls[chrom][i] == 0 else 1] += 1
    return counts, poly_keys, n_poly


def write_counts(path, counts, genes):
    with open(path, "w") as out:
        out.write("gene\tpn\tps\n")
        for gid in sorted(counts, key=genes.__getitem__):
            pn, ps = counts[gid]
            out.write(f"{genes[gid]}\t{pn}\t{ps}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map", required=True)
    ap.add_argument("--zero", required=True)
    ap.add_argument("--four", required=True)
    ap.add_argument("--freq", action="append", required=True,
                    help="ID=path of a per-pool allele-frequency table")
    ap.add_argument("--min-depth", type=int, default=10)
    ap.add_argument("--min-alt-count", type=int, default=1)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    pools = []
    for spec in args.freq:
        pool, sep, path = spec.partition("=")
        if not sep or not pool or not path:
            sys.exit(f"bad --freq value: {spec!r} (expected ID=path)")
        pools.append((pool, path))

    zero = read_bed_arrays(args.zero)
    four = read_bed_arrays(args.four)
    index = build_index(args.map, zero, four)
    del zero, four
    folded_pos, folded_cls, folded_gene, cds_pos, genes, n_rows = index
    n_folded = sum(len(a) for a in folded_pos.values())
    n_cds = sum(len(a) for a in cds_pos.values())
    if not n_folded:
        sys.exit("no CDS sites overlap the degeneracy BEDs")
    print(f"{n_folded} 0/4-fold CDS sites in {len(genes)} genes "
          f"({n_cds} CDS sites from {n_rows} map rows)")

    os.makedirs(args.outdir, exist_ok=True)
    pooled_keys = set()
    chrom_ids = {}        # assigned lazily, only scaffolds with polymorphic sites

    for pool, path in pools:
        counts, keys, n_poly = count_pool(
            path, index, args.min_depth, args.min_alt_count, chrom_ids)
        write_counts(os.path.join(args.outdir, f"{pool}_polymorphism.tsv"),
                     counts, genes)
        pooled_keys |= keys
        print(f"pool {pool}: {n_poly} polymorphic CDS sites in "
              f"{len(counts)} genes")

    # pooled table: every site counts once, even if polymorphic in several
    # pools; only 0/4-fold sites enter Pn/Ps
    chrom_names = {cid: chrom for chrom, cid in chrom_ids.items()}
    pooled_counts = {}
    for key in pooled_keys:
        cid, pos = key >> POS_BITS, key & POS_MASK
        chrom = chrom_names[cid]
        arr = folded_pos.get(chrom)
        if not arr:
            continue
        i = bisect.bisect_left(arr, pos)
        if i >= len(arr) or arr[i] != pos:
            continue
        gid = folded_gene[chrom][i]
        entry = pooled_counts.get(gid)
        if entry is None:
            entry = pooled_counts[gid] = [0, 0]
        entry[0 if folded_cls[chrom][i] == 0 else 1] += 1
    write_counts(os.path.join(args.outdir, "pooled_polymorphism.tsv"),
                 pooled_counts, genes)

    chrom_names = {cid: chrom for chrom, cid in chrom_ids.items()}
    with open(os.path.join(args.outdir, "polymorphic_sites.tsv"), "w") as out:
        out.write("chrom\tpos\n")
        for key in sorted(pooled_keys,
                          key=lambda k: (chrom_names[k >> POS_BITS],
                                         k & POS_MASK)):
            out.write(f"{chrom_names[key >> POS_BITS]}\t{key & POS_MASK}\n")
    print(f"{len(pooled_keys)} polymorphic CDS sites (union of pools)")


if __name__ == "__main__":
    main()
