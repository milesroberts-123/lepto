#!/usr/bin/env bash
# Join freqk per-pool allele frequencies to grenedalf per-site frequencies
# at the sites freqk indexed, for one pool.
#
# Both sides are position-sorted, so this is a streaming merge: the freqk
# keyed table (3.5 M rows at this annotation) is held in awk memory
# (~300-500 MB) and the grenedalf CSV (tens of GB, a row per dict site
# because it is written with --write-invariants) is read exactly once.
#
# Only rows whose REF/ALT agree between the two callers are written; the
# count of disagreeing rows is reported on stderr so a caller can say what
# fraction of sites the comparison covers.
#
# Usage: freqk_vs_grenedalf.sh <workflow_dir> <pool_id> <out_tsv>
set -euo pipefail

wf=${1:?workflow dir}
pool=${2:?pool id}
out=${3:?output tsv}

fk_dir="$wf/results/freqk/hap2_vs_hap1"
gren_csv="$wf/results/vg/hap2_vs_hap1/$pool/grenedalf_results_frequency.csv"
index="$fk_dir/filtered_index.txt"
calls="$fk_dir/${pool}_calls.txt"

for f in "$index" "$calls" "$gren_csv"; do
    [ -s "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done

mkdir -p "$(dirname "$out")"

# 1. freqk keyed table: chrom, pos (1-based), ref, alt, alt_af
#    index col2 = chrom, col3 = pos (0-based, one less than VCF POS),
#    col5 = REF|ALT; the call file is one line per index row, col2 = ALT af.
paste <(awk -F, '{print $2"\t"$3+2"\t"$5}' "$index") "$calls" |
    awk -F'\t' 'BEGIN{OFS="\t"}
        {split($3,a,"|"); split($4,c,"|"); print $1, $2, a[1], a[2], c[2]}' \
    > "$out.freqk"

# 2. stream grenedalf once, emit matching rows
awk -F, -v OFS='\t' -v out="$out" '
    NR == FNR { key = $1 SUBSEP $2; ref[key] = $3; alt[key] = $4; af[key] = $5; next }
    $1 == "CHROM" { next }
    {
        key = $1 SUBSEP $2
        if (!(key in af)) next
        if ($3 != ref[key] || $4 != alt[key]) { mismatch++; next }
        print $1, $2, ref[key], alt[key], af[key], $8, $6, $7
        n++
    }
    END {
        printf "freqk_vs_grenedalf: %d matched rows, %d REF/ALT mismatches\n",
               n, mismatch + 0 > "/dev/stderr"
    }
' "$out.freqk" "$gren_csv" > "$out.tmp"

# header + body
{
    printf 'chrom\tpos\tref\talt\tfreqk_af\tgren_af\tgren_alt_cnt\tgren_depth\n'
    cat "$out.tmp"
} > "$out"

rm -f "$out.freqk" "$out.tmp"
echo "wrote $out" >&2
