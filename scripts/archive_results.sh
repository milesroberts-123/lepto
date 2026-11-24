#!/usr/bin/env bash
# Bundle the keep-list of workflow/results subdirectories for cold storage
# (S3 Glacier Deep Archive) and write one sha256 per bundle.
#
# Usage: scripts/archive_results.sh [--force] [--verify]
#
#   --force   rebuild bundles that already exist
#   --verify  run sha256sum -c over every bundle after the run
#
# Output: workflow/bundles/<dir>.tar.zst and <dir>.tar.zst.sha256
#
# What is deliberately NOT bundled, and why:
#   kmate/kmerdb   >100 GB per pool, rebuilt in ~9 h from the reads
#   vg, ldx        graph/LD intermediates, deterministic from the assembly + VCFs
#   rna, hisat2    re-downloadable/re-mappable
#   intersect, union, fastp   k-mer set algebra and trimmed reads
#   results/backup.tar.gz     a previous tar of these same directories
#   repo-root results/        notebook by-date/by-analysis output, re-knittable
#
# Upload afterwards from a LOGIN node (compute nodes have no internet):
#   rclone copy workflow/bundles/ remote:lepto-archive/bundles/ \
#       --s3-storage-class DEEP_ARCHIVE --s3-chunk-size 64M \
#       --transfers 8 --checkers 16 --immutable --stats 1m --progress
set -euo pipefail

force=0
verify=0
for arg in "$@"; do
    case "$arg" in
        --force) force=1 ;;
        --verify) verify=1 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
results=$repo_root/workflow/results
bundles=$repo_root/workflow/bundles
threads=${ZSTD_THREADS:-4}

command -v zstd >/dev/null || { echo "zstd not found on PATH" >&2; exit 1; }
command -v sha256sum >/dev/null || { echo "sha256sum not found" >&2; exit 1; }
[ -d "$results" ] || { echo "no results directory: $results" >&2; exit 1; }

mkdir -p "$bundles"
# A .tmp can only be the debris of an interrupted run (bundles are renamed into
# place only after tar|zstd finishes), so it is always safe to clear.
rm -f "$bundles"/*.tar.zst.tmp

# zstd -9 for text (alignments, GFF, FASTA, TSV); -3 for content that is
# already compressed (BAM, VCF.gz, TSV.gz), where -9 only burns CPU.
text_dirs=(
    braker busco mk orthofinder panther degenotate msmc2 anchorwave
    mcscanx cds_windows specific dnds vigs tubemap bwa sourmash kmc grouped
)
packed_dirs=(freqk svim featurecounts)

make_bundle() {
    local name=$1 level=$2
    shift 2
    local out=$bundles/$name.tar.zst
    if [ ! -e "$results/$name" ]; then
        echo "skip (missing): results/$name"
        return 0
    fi
    if [ -e "$out" ] && [ "$force" -eq 0 ]; then
        echo "skip (exists): ${out#"$repo_root"/}  [use --force to rebuild]"
        return 0
    fi
    echo "==> $name (zstd -$level)"
    # Write to .tmp and rename so an interrupted run never leaves a bundle
    # that looks complete; the sha256 is only written alongside a finished file.
    tar -C "$results" -cf - "$@" | zstd -T"$threads" "-$level" >"$out.tmp"
    mv "$out.tmp" "$out"
    sha256sum "$out" >"$out.sha256"
}

for d in "${text_dirs[@]}"; do
    make_bundle "$d" 9 "$d"
done
for d in "${packed_dirs[@]}"; do
    make_bundle "$d" 3 "$d"
done
# kmate without the jellyfish DBs: those are the bulk of the directory and the
# only part that is cheap to regenerate relative to its size.
make_bundle kmate 3 --exclude='kmate/kmerdb' kmate

if [ "$verify" -eq 1 ]; then
    echo "==> verifying"
    (cd "$bundles" && sha256sum -c ./*.sha256)
fi

echo
echo "bundles in ${bundles#"$repo_root"/}:"
du -ch "$bundles"/*.tar.zst | tail -1
echo "checksums: cd ${bundles#"$repo_root"/} && sha256sum -c ./*.sha256"
