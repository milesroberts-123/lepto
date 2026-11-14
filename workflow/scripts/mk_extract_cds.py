#!/usr/bin/env python3
"""Extract .t1 CDS sequences and a CDS-position -> genome map from a GFF3.

Only CDS features whose Parent transcript ends in the transcript suffix
(default ".t1") are kept. Segments are concatenated in transcript
orientation; the map file records, for every CDS position (1-based, in
transcript order), its scaffold and 1-based genomic position, so codon
alignments can be mapped back onto reference coordinates.

Outputs:
  --fna  FASTA, one record per transcript (ID = transcript)
  --map  TSV: transcript gene cdspos chrom pos
"""

import argparse
import re
import sys
from collections import defaultdict

import pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def gene_of(transcript):
    return re.sub(r"\.t\d+$", "", transcript)


def parse_cds(gff, suffix):
    # transcript -> list of (chrom, start, end, strand); GFF is 1-based
    segments = defaultdict(list)
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r"(?:^|;)Parent=([^;]+)", f[8])
            if not m:
                continue
            for parent in m.group(1).split(","):
                parent = parent.strip()
                if not parent.endswith(suffix):
                    continue
                segments[parent].append(
                    (f[0], int(f[3]), int(f[4]), f[6]))
    return segments


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--fna", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--transcript-suffix", default=".t1")
    ap.add_argument("--chrom-prefix", default="",
                    help="prefix added to scaffold names in the map "
                         "(e.g. 'hap1_' to match grenedalf/degenotate coords)")
    args = ap.parse_args()

    segments = parse_cds(args.gff, args.transcript_suffix)
    if not segments:
        sys.exit(f"no CDS features with Parent *{args.transcript_suffix} in {args.gff}")

    fasta = pysam.FastaFile(args.fasta)
    contigs = set(fasta.references)

    n_ok = n_bad = 0
    with open(args.fna, "w") as fna, open(args.map, "w") as mp:
        mp.write("transcript\tgene\tcdspos\tchrom\tpos\n")
        for transcript in sorted(segments):
            segs = segments[transcript]
            chrom = segs[0][0]
            strand = segs[0][3]
            if chrom not in contigs:
                n_bad += 1
                continue
            if strand == "+":
                segs.sort(key=lambda s: s[1])
                strand = "+"
            else:
                segs.sort(key=lambda s: s[1], reverse=True)
                strand = "-"

            seq_parts = []
            positions = []  # genomic positions in transcript order
            for _, start, end, _ in segs:
                seq_parts.append(fasta.fetch(chrom, start - 1, end))
                if strand == "+":
                    positions.extend(range(start, end + 1))
                else:
                    positions.extend(range(end, start - 1, -1))

            seq = "".join(seq_parts).upper()
            if strand == "-":
                seq = revcomp(seq)
            if not seq:
                n_bad += 1
                continue

            fna.write(f">{transcript}\n{seq}\n")
            gene = gene_of(transcript)
            for i, pos in enumerate(positions, start=1):
                mp.write(f"{transcript}\t{gene}\t{i}\t{args.chrom_prefix}{chrom}\t{pos}\n")
            n_ok += 1

    fasta.close()
    print(f"{args.species}: {n_ok} CDS extracted, {n_bad} skipped")
    if not n_ok:
        sys.exit("no CDS sequences extracted")


if __name__ == "__main__":
    main()
