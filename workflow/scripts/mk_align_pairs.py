#!/usr/bin/env python3
"""Codon-align one 1:1 ortholog pair (MAFFT proteins -> pal2nal codons).

Invoked once per orthogroup (the rule drives it with xargs -P). Writes
<outdir>/<orthogroup>.fna on success, <outdir>/<orthogroup>.failed with a
reason otherwise. Alignment headers are "<species>|<transcript>" so the
ingroup sequence is identifiable downstream.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mk_common import read_fasta  # noqa: E402


def load_pairs(path):
    pairs = {}
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                pairs[f[0]] = (f[1], f[2])
    return pairs


def fail(outdir, og, reason):
    with open(os.path.join(outdir, og + ".failed"), "w") as fh:
        fh.write(reason + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--only", required=True, help="orthogroup to align")
    ap.add_argument("--ingroup", required=True)
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--ingroup-fna", required=True)
    ap.add_argument("--outgroup-fna", required=True)
    ap.add_argument("--ingroup-faa", required=True)
    ap.add_argument("--outgroup-faa", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-aa", type=int, default=30)
    args = ap.parse_args()

    pairs = load_pairs(args.pairs)
    if args.only not in pairs:
        sys.exit(f"orthogroup {args.only} not in {args.pairs}")
    in_tx, out_tx = pairs[args.only]

    cds = {}
    cds.update(read_fasta(args.ingroup_fna))
    cds.update(read_fasta(args.outgroup_fna))
    prot = {}
    prot.update(read_fasta(args.ingroup_faa))
    prot.update(read_fasta(args.outgroup_faa))

    if in_tx not in cds or out_tx not in cds:
        return fail(args.outdir, args.only, "missing_cds")
    if in_tx not in prot or out_tx not in prot:
        return fail(args.outdir, args.only, "missing_protein")

    in_name = f"{args.ingroup}|{in_tx}"
    out_name = f"{args.outgroup}|{out_tx}"
    cds_records = {in_name: cds[in_tx], out_name: cds[out_tx]}
    prot_records = {in_name: prot[in_tx], out_name: prot[out_tx]}

    if min(len(p) for p in prot_records.values()) < args.min_aa:
        return fail(args.outdir, args.only, "short_protein")

    mafft = shutil.which("mafft")
    pal2nal = shutil.which("pal2nal.pl") or shutil.which("pal2nal")
    if mafft is None or pal2nal is None:
        sys.exit("mafft and pal2nal.pl must be on PATH (paml env)")

    with tempfile.TemporaryDirectory() as tmp:
        aa_in = os.path.join(tmp, "aa.fa")
        cds_in = os.path.join(tmp, "cds.fa")
        aa_aln = os.path.join(tmp, "aa.aln")
        with open(aa_in, "w") as fh:
            for name, seq in prot_records.items():
                fh.write(f">{name}\n{seq}\n")
        with open(cds_in, "w") as fh:
            for name, seq in cds_records.items():
                fh.write(f">{name}\n{seq}\n")

        with open(aa_aln, "w") as out:
            r = subprocess.run([mafft, "--quiet", "--auto", aa_in],
                               stdout=out, stderr=subprocess.DEVNULL)
        if r.returncode != 0 or os.path.getsize(aa_aln) == 0:
            return fail(args.outdir, args.only, "mafft_failed")

        r = subprocess.run([pal2nal, aa_aln, cds_in, "-output", "fasta"],
                           capture_output=True, text=True)
        if r.returncode != 0 or not r.stdout.strip():
            return fail(args.outdir, args.only, "pal2nal_failed")

    aln = read_fasta_text(r.stdout)
    if sorted(aln) != sorted(cds_records):
        return fail(args.outdir, args.only, "pal2nal_names")
    lens = {len(v) for v in aln.values()}
    if len(lens) != 1 or lens.pop() % 3 != 0:
        return fail(args.outdir, args.only, "frameshift")

    with open(os.path.join(args.outdir, args.only + ".fna"), "w") as out:
        for name in (in_name, out_name):
            out.write(f">{name}\n{aln[name]}\n")


def read_fasta_text(text):
    seqs = {}
    name = None
    chunks = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


if __name__ == "__main__":
    main()
