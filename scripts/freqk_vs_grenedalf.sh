#!/usr/bin/env bash
# Join freqk per-pool allele frequencies to grenedalf per-site frequencies
# at the sites freqk indexed, for one pool.
#
# Both sides are position-sorted, so this is a streaming merge: the freqk
# keyed table (3.5 M rows at this annotation) is held in awk memory
# (~300-500 MB) and the grenedalf CSV (tens of GB, a row per dict site
# because it is written with --write-invariants) is read exactly once.
#
# Allele-labelling caveat: grenedalf is given a sequence dictionary, not the
# genome, so it labels the two observed bases itself. Where the VCF ALT is
# the majority base the pair comes back swapped (REF=<VCF ALT>), and the
# join flips the frequency accordingly; sites whose bases are a third allele
# or N cannot be compared and are counted, per freqk-AF bin, on stderr.
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

# Only rows whose allele pair matches are written. grenedalf does not know
# the reference base (its input is a dict, not the genome), so it labels the
# two observed bases itself: at sites where the VCF ALT is the majority base
# grenedalf reports REF=<VCF ALT>, ALT=<VCF REF>, i.e. the pair is swapped.
# Accepting that orientation (and flipping the frequency: freq(VCF ALT) =
# 1 - FREQ) is what makes the comparison valid at high ALT frequency;
# without it every near-fixed site silently drops out and the scatter looks
# empty in the upper right. Rows whose two bases are not the VCF pair at all
# (third allele, N) stay unmatched and are counted by freqk AF bin.
asc="${out%.tsv}.ascertainment.tsv"
awk -v OFS='\t' -v asc="$asc" '
    BEGIN { FS = "\t" }
    NR == FNR {
        key = $1 SUBSEP $2; ref[key] = $3; alt[key] = $4; af[key] = $5;
        nkey++
        next
    }
    FNR == 1 { FS = "," }
    $1 == "CHROM" { next }
    {
        key = $1 SUBSEP $2
        if (!(key in af)) next
        seen[key] = 1
        if ($3 == ref[key] && $4 == alt[key]) {
            f = $8 + 0; altcnt = $6; swapped = 0
        } else if ($3 == alt[key] && $4 == ref[key]) {
            # gren REF is the VCF ALT here, so REF_CNT ($5) is the VCF ALT count
            f = 1 - ($8 + 0); altcnt = $5; swapped = 1
        } else {
            unmatched[key] = 1
            next
        }
        print $1, $2, ref[key], alt[key], af[key], f, altcnt, $7, swapped
        n++; nswap += swapped
    }
    END {
        for (k in af) if (!(k in seen)) unmatched[k] = 1
        nu = 0; for (k in unmatched) nu++
        printf "freqk_vs_grenedalf: %d matched rows (%d allele-swapped), %d of %d freqk sites unmatched\n",
               n + 0, nswap + 0, nu, nkey + 0 > "/dev/stderr"
        # ascertainment: distribution of freqk AF among matched vs unmatched
        split("0.0-0.2 0.2-0.4 0.4-0.6 0.6-0.8 0.8-1.0", lab, " ")
        for (k in af) {
            b = int(af[k] / 0.2); if (b > 4) b = 4
            if (k in unmatched) unm[b]++; else mat[b]++
        }
        printf "freqk_af_bin\tmatched\tunmatched\n" > asc
        for (b = 0; b <= 4; b++)
            printf "%s\t%d\t%d\n", lab[b+1], mat[b] + 0, unm[b] + 0 > asc
        close(asc)
    }
' "$out.freqk" "$gren_csv" > "$out.tmp"

# header + body
{
    printf 'chrom\tpos\tref\talt\tfreqk_af\tgren_af\tgren_alt_cnt\tgren_depth\tgren_swapped\n'
    cat "$out.tmp"
} > "$out"

rm -f "$out.freqk" "$out.tmp"
echo "wrote $out" >&2
