#!/usr/bin/env python3
"""Parse OrthoFinder 1:1 orthologs into an ingroup/outgroup pair table.

Reads <orthofinder>/Orthogroups/Orthogroups.tsv (species columns) and
keeps orthogroups with exactly one gene in both requested species whose
ID is a .t1 transcript. OrthoFinder proteome IDs carry a "<species>."
prefix (see orthofinder.smk orthofinder_unique_prefix), stripped here.

Output columns (tab-separated, header):
  orthogroup  ingroup_transcript  outgroup_transcript  ingroup_gene  outgroup_gene
"""

import argparse
import re
import sys


def strip_species(ids, species):
    prefix = species + "."
    return [i[len(prefix):] if i.startswith(prefix) else i for i in ids]


def split_cell(cell):
    # OrthoFinder separates multiple IDs in a cell with ", " (and ",")
    return [x.strip() for x in re.split(r",\s*", cell.strip()) if x.strip()]


def gene_of(transcript):
    return re.sub(r"\.t\d+$", "", transcript)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--orthofinder", required=True,
                    help="OrthoFinder results directory")
    ap.add_argument("--ingroup", required=True)
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    tsv = f"{args.orthofinder}/Orthogroups/Orthogroups.tsv"
    try:
        fh = open(tsv)
    except OSError:
        sys.exit(f"cannot open {tsv}")

    with fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            og_col = 0
            in_col = header.index(args.ingroup)
            out_col = header.index(args.outgroup)
        except ValueError:
            sys.exit(f"{tsv} lacks species columns {args.ingroup}/{args.outgroup}: {header}")

        rows = []
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(in_col, out_col):
                continue
            ing = strip_species(split_cell(fields[in_col]), args.ingroup)
            out = strip_species(split_cell(fields[out_col]), args.outgroup)
            if len(ing) != 1 or len(out) != 1:
                continue
            if not ing[0].endswith(".t1") or not out[0].endswith(".t1"):
                continue
            rows.append((fields[og_col], ing[0], out[0],
                         gene_of(ing[0]), gene_of(out[0])))

    with open(args.out, "w") as out:
        out.write("orthogroup\tingroup_transcript\toutgroup_transcript\t"
                  "ingroup_gene\toutgroup_gene\n")
        for row in rows:
            out.write("\t".join(row) + "\n")

    if not rows:
        sys.exit("no 1:1 .t1 ortholog pairs found")
    print(f"{len(rows)} 1:1 ortholog pairs -> {args.out}")


if __name__ == "__main__":
    main()
