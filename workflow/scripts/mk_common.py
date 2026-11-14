#!/usr/bin/env python3
"""Shared helpers for the McDonald-Kreitman scripts."""

import itertools
import re

# standard genetic code (NCBI table 1), TAA/TAG/TGA stops
_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
CODON_TABLE = {
    a + b + c: _AAS[i]
    for i, (a, b, c) in enumerate(itertools.product(_BASES, repeat=3))
}


def read_fasta(path, upper=True):
    seqs = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    if upper:
        seqs = {k: v.upper() for k, v in seqs.items()}
    return seqs


def gene_of(transcript):
    return re.sub(r"\.t\d+$", "", transcript)


def translate(codon):
    return CODON_TABLE.get(codon, "X")


def classify_codon_change(codon_from, codon_to):
    """Return (n_syn, n_nonsyn) averaged over minimum-replacement paths.

    A single-position difference is classified directly. Codons differing
    at several positions are connected by every shortest mutational path
    (degenotate / mkado convention); each path's synonymous and
    nonsynonymous steps are counted and averaged equally across paths.
    """
    if codon_from == codon_to:
        return (0.0, 0.0)
    diff = [i for i in range(3) if codon_from[i] != codon_to[i]]
    if len(diff) == 1:
        return (1.0, 0.0) if translate(codon_from) == translate(codon_to) else (0.0, 1.0)

    syn = nonsyn = 0.0
    paths = 0
    for order in itertools.permutations(diff):
        current = codon_from
        p_syn = p_non = 0.0
        ok = True
        for i in order:
            nxt = current[:i] + codon_to[i] + current[i + 1:]
            if translate(nxt) == "X":
                ok = False
                break
            if translate(current) == translate(nxt):
                p_syn += 1.0
            else:
                p_non += 1.0
            current = nxt
        if not ok:
            continue
        syn += p_syn
        nonsyn += p_non
        paths += 1
    if paths == 0:
        return (0.0, 0.0)
    return (syn / paths, nonsyn / paths)


def split_alignment(seqs, ingroup, outgroup):
    """Return (ingroup_seq, outgroup_seq, ingroup_transcript)."""
    ing = [(k, v) for k, v in seqs.items() if k.startswith(ingroup + "|")]
    out = [(k, v) for k, v in seqs.items() if k.startswith(outgroup + "|")]
    if len(ing) != 1 or len(out) != 1:
        return None
    return ing[0][1], out[0][1], ing[0][0].split("|", 1)[1]
