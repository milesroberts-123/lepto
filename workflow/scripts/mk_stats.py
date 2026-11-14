#!/usr/bin/env python3
"""Standard McDonald-Kreitman statistics per gene and per pool.

Joins the divergence table (Dn/Ds, fractional for multi-change codons)
with a polymorphism table (Pn/Ps, integers) and reports, per gene:

  alpha = 1 - (Ds * Pn) / (Dn * Ps)      (undefined -> NA)
  NI    = (Pn * Ds) / (Ps * Dn)          (neutrality index)
  DoS   = Pn/(Pn+Ps) - Dn/(Dn+Ds)
  fisher_p: two-sided Fisher exact test on [[Pn, Ps], [Dn, Ds]]

Because minimum-replacement-path averaging can leave fractional Dn/Ds,
the Fisher test uses rounded integer counts (reported alongside).
Benjamini-Hochberg q-values are computed within each pool.

Outputs (in --outdir): mk_stats_<ID>.tsv per pool, mk_stats_pooled.tsv.
"""

import argparse
import math
import os
import sys


def fisher_two_sided(a, b, c, d):
    """Fisher exact test p-value for a 2x2 table [[a, b], [c, d]]."""
    n1 = a + b
    n2 = c + d
    total = n1 + n2
    if total == 0:
        return float("nan")

    def hyper(i):
        # P(table with first-cell = i), fixed margins
        return (math.comb(n1, i) * math.comb(n2, a + c - i)) / math.comb(total, a + c)

    lo = max(0, n1 + a + c - total)
    hi = min(n1, a + c)
    observed = hyper(a)
    tol = observed * (1 + 1e-9)
    p = sum(prob for i in range(lo, hi + 1) if (prob := hyper(i)) <= tol)
    return min(p, 1.0)


def bh_adjust(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [float("nan")] * n
    running = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        running = min(running, pvals[i] * n / (rank + 1))
        q[i] = min(running, 1.0)
    return q


def read_divergence(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            rows.append({
                "gene": f[idx["gene"]],
                "orthogroup": f[idx["orthogroup"]],
                "dn": float(f[idx["dn"]]),
                "ds": float(f[idx["ds"]]),
            })
    return rows


def read_polymorphism(path):
    counts = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                counts[f[0]] = (int(f[1]), int(f[2]))
    return counts


def compute_table(divergence, poly):
    rows = []
    for rec in divergence:
        pn, ps = poly.get(rec["gene"], (0, 0))
        dn, ds = rec["dn"], rec["ds"]
        alpha = (1.0 - (ds * pn) / (dn * ps)) if dn > 0 and ps > 0 else float("nan")
        ni = (pn * ds) / (ps * dn) if ps > 0 and dn > 0 else float("nan")
        dos = (pn / (pn + ps) if pn + ps > 0 else float("nan"))
        dos -= (dn / (dn + ds) if dn + ds > 0 else float("nan"))
        dn_i, ds_i = int(round(dn)), int(round(ds))
        p = fisher_two_sided(pn, ps, dn_i, ds_i)
        rows.append({
            "gene": rec["gene"], "orthogroup": rec["orthogroup"],
            "dn": dn, "ds": ds, "pn": pn, "ps": ps,
            "dn_int": dn_i, "ds_int": ds_i,
            "alpha": alpha, "ni": ni, "dos": dos, "fisher_p": p,
        })
    finite = [i for i, r in enumerate(rows) if not math.isnan(r["fisher_p"])]
    qs = bh_adjust([rows[i]["fisher_p"] for i in finite])
    for i, q in zip(finite, qs):
        rows[i]["fisher_q"] = q
    for r in rows:
        r.setdefault("fisher_q", float("nan"))
    return rows


def fmt(x):
    return "NA" if isinstance(x, float) and math.isnan(x) else \
        (f"{x:.10g}" if isinstance(x, float) else str(x))


def write_table(path, rows):
    cols = ["gene", "orthogroup", "dn", "ds", "pn", "ps", "dn_int", "ds_int",
            "alpha", "ni", "dos", "fisher_p", "fisher_q"]
    with open(path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(fmt(r[c]) for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--divergence", required=True)
    ap.add_argument("--polymorphism", action="append", required=True,
                    help="ID=path of a <ID>_polymorphism.tsv table")
    ap.add_argument("--pooled", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    divergence = read_divergence(args.divergence)
    if not divergence:
        sys.exit("empty divergence table")

    os.makedirs(args.outdir, exist_ok=True)
    for spec in args.polymorphism:
        pool, sep, path = spec.partition("=")
        if not sep:
            sys.exit(f"bad --polymorphism value: {spec!r} (expected ID=path)")
        rows = compute_table(divergence, read_polymorphism(path))
        write_table(os.path.join(args.outdir, f"mk_stats_{pool}.tsv"), rows)
        print(f"pool {pool}: {len(rows)} genes")

    rows = compute_table(divergence, read_polymorphism(args.pooled))
    write_table(os.path.join(args.outdir, "mk_stats_pooled.tsv"), rows)
    print(f"pooled: {len(rows)} genes")


if __name__ == "__main__":
    main()
