#!/usr/bin/env bash
# Per-bin read depth from a pool's grenedalf frequency table.
#
# The frequency CSV has a row per dictionary site (it is written with
# --write-invariants), with REF_CNT/ALT_CNT/DEPTH in columns 5/6/7. One
# streaming pass over the file (tens of GB) accumulates depth sums into
# fixed-width bins, so the notebook can rebuild any window size that is a
# multiple of the bin width without touching the CSV again.
#
# Depth here is post-filter (the same map/base-quality thresholds that
# produced theta_pi), so it is the right covariate for asking whether
# CDS-dense windows are simply better covered.
#
# Output (tab, header): chrom, bin_start (1-based), depth_sum, n_sites
#
# Usage: grenedalf_window_depth.sh <workflow_dir> <pool_id> <out_tsv> [bin_bp]
set -euo pipefail

wf=${1:?workflow dir}
pool=${2:?pool id}
out=${3:?output tsv}
bin=${4:-1000}

csv="$wf/results/vg/hap2_vs_hap1/$pool/grenedalf_results_frequency.csv"
[ -s "$csv" ] || { echo "missing input: $csv" >&2; exit 1; }

mkdir -p "$(dirname "$out")"

awk -F, -v OFS='\t' -v bin="$bin" '
    BEGIN { print "chrom", "bin_start", "depth_sum", "n_sites" }
    $1 == "CHROM" { next }
    {
        b = int((($2 + 0) - 1) / bin) * bin + 1
        key = $1 SUBSEP b
        depth[key] += $7 + 0
        n[key]++
    }
    END {
        for (k in depth) {
            split(k, a, SUBSEP)
            print a[1], a[2], depth[k], n[k]
        }
    }
' "$csv" | sort -k1,1 -k2,2n > "$out"

echo "wrote $out ($(( $(wc -l < "$out") - 1 )) bins)" >&2
