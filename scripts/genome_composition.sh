#!/usr/bin/env bash
# Assembly composition summary: sequence count, total length, N50, and base
# composition (uppercase / lowercase = soft-masked / N / other IUPAC), which
# is what decides whether BRAKER's --softmasking flag is meaningful.
#
# Usage: genome_composition.sh <fasta> [fai]
set -euo pipefail

fa=${1:?fasta}
fai=${2:-$fa.fai}

[ -s "$fa" ] || { echo "missing $fa" >&2; exit 1; }
[ -s "$fai" ] || { echo "missing $fai (run samtools faidx)" >&2; exit 1; }

awk -F'\t' -v OFS='\t' '{n++; len[$1]=$2; tot+=$2}
    END{
        asort(len, sorted)
        half = tot/2; cum = 0
        for (i = length(sorted); i >= 1; i--) {
            cum += sorted[i]
            if (cum >= half) { n50 = sorted[i]; l50 = length(sorted) - i + 1; break }
        }
        printf "sequences\t%d\ntotal_bp\t%d\nlongest_bp\t%d\nN50_bp\t%d\nL50\t%d\n", n, tot, sorted[length(sorted)], n50, l50
    }' "$fai"

perl -ne 'next if /^>/; chomp; $lo += ($_ =~ tr/acgt//); $up += ($_ =~ tr/ACGT//); $n += ($_ =~ tr/Nn//); $len += length($_); END{ $other = $len - $lo - $up - $n; printf "total_bases\t%d\nlowercase\t%d\t%.4f\nuppercase\t%d\t%.4f\nN\t%d\t%.4f\nother_iupac\t%d\t%.6f\n", $len, $lo, $lo/$len, $up, $up/$len, $n, $n/$len, $other, $other/$len }' "$fa"
