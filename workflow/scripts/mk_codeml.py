#!/usr/bin/env python3
"""Run pairwise codeml (runmode -2) on one codon alignment.

Writes <outdir>/<orthogroup>.tsv with one row: orthogroup gene dN dS
omega t status. Alignments with a stop codon or a non-multiple-of-three
length are reported with status=stop_codon/frameshift instead of being
run (codeml aborts on those).
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mk_common import (gene_of, read_fasta, split_alignment,  # noqa: E402
                       translate)

CTL = """\
seqfile = {seqfile}
treefile = {treefile}
outfile = {outfile}

noisy = 0
verbose = 0
runmode = -2
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
icode = 0
cleandata = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
"""

NUM = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"


def write_alignment(path, seqs):
    names = list(seqs)
    length = len(seqs[names[0]])
    with open(path, "w") as out:
        out.write(f" 2 {length}\n")
        for i, name in enumerate(names, start=1):
            # PAML strict PHYLIP: name padded/truncated to 10 characters
            out.write(f"{name:<10}{seqs[name]}\n")
    return names


def parse_mlc(path):
    values = {}
    with open(path, errors="replace") as fh:
        text = fh.read()
    for key, pattern in (
        ("omega", rf"omega \(dN/dS\)\s*=\s*({NUM})"),
        ("dN/dS", rf"dN/dS=\s*({NUM})"),
        ("dN", rf"(?<![\w/])dN\s*=\s*({NUM})"),
        ("dS", rf"(?<![\w/])dS\s*=\s*({NUM})"),
        ("t", rf"t=\s*({NUM})"),
    ):
        m = re.search(pattern, text)
        if m:
            values[key] = float(m.group(1))
    return values


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--ingroup", required=True)
    ap.add_argument("--outgroup", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    og = os.path.basename(args.alignment)[:-4]
    seqs = read_fasta(args.alignment)
    split = split_alignment(seqs, args.ingroup, args.outgroup)
    out_path = os.path.join(args.outdir, og + ".tsv")

    if split is None:
        with open(out_path, "w") as out:
            out.write(f"{og}\tNA\tNA\tNA\tNA\tNA\talignment_names\n")
        return
    ing, out_seq, tx = split
    gene = gene_of(tx)

    def emit(dn, ds, omega, t, status):
        with open(out_path, "w") as out:
            out.write(f"{og}\t{gene}\t{dn}\t{ds}\t{omega}\t{t}\t{status}\n")

    if len(ing) != len(out_seq) or len(ing) % 3 != 0:
        return emit("NA", "NA", "NA", "NA", "frameshift")
    if any(translate(ing[i:i + 3]) == "*" for i in range(0, len(ing), 3)) or \
            any(translate(out_seq[i:i + 3]) == "*" for i in range(0, len(out_seq), 3)):
        return emit("NA", "NA", "NA", "NA", "stop_codon")

    with tempfile.TemporaryDirectory() as tmp:
        seqfile = os.path.join(tmp, "aln.phy")
        treefile = os.path.join(tmp, "tree.nwk")
        outfile = os.path.join(tmp, "mlc")
        ctlfile = os.path.join(tmp, "codeml.ctl")
        names = write_alignment(seqfile, {"hap1": ing, "linanthus": out_seq})
        with open(treefile, "w") as out:
            out.write(f"({names[0]},{names[1]});\n")
        with open(ctlfile, "w") as out:
            out.write(CTL.format(seqfile=seqfile, treefile=treefile,
                                 outfile=outfile))
        r = subprocess.run(["codeml", ctlfile], cwd=tmp,
                           capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(outfile):
            return emit("NA", "NA", "NA", "NA", "codeml_failed")
        values = parse_mlc(outfile)

    if "dN" not in values or "dS" not in values:
        return emit("NA", "NA", "NA", "NA", "no_numbers")
    emit(f"{values['dN']:.6f}", f"{values['dS']:.6f}",
         f"{values.get('omega', values.get('dN/dS', float('nan'))):.6f}",
         f"{values.get('t', float('nan')):.6f}", "ok")


if __name__ == "__main__":
    main()
