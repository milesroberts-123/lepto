# kMate panel building (`kmate.smk`)

`kmate build-index` (checkpoint `kmate_build_index`) turns the freqk-derived
`panel.vcf.gz` plus the hap1 backbone into the per-chromosome
`panel_<chrom>_kmers.tsv.gz` files that `kmate_build_kmer_pa` consumes.

## Node-local temp is mandatory

kmate 0.1.2 has **no `--tmp-dir`**. `JellyfishCounter` derives its work
directory as `<out>_tmp`, writes the jellyfish hash to `<out>_tmp/panel.jf`,
and then queries that hash through `dna_jellyfish.QueryMerFile`, i.e. an
`mmap` with random access. On Lustre every page fault is a network round
trip, so the query stage stalls; kMate's own `BUILDING_A_PANEL.md` §Step 1
warns it "can stall for hours". The symptom is a job that sits at 44 GB RSS
and 4 threads for its whole wall time with no output.

Because the prefix names both the hash and the TSVs, the rule redirects
`--out` itself:

```bash
mkdir -p ${SLURM_TMPDIR:-/tmp}/kmate_index
kmate build-index --out ${SLURM_TMPDIR:-/tmp}/kmate_index/panel ...
mkdir -p $(dirname {output})/index && \
    cp ${SLURM_TMPDIR:-/tmp}/kmate_index/panel_*_kmers.tsv.gz $(dirname {output})/index/
```

`SLURM_TMPDIR` is node-local NVMe on Savio, so the `mmap` reads stay on the
box. The declared `output:` directory is unchanged — the copy is what makes
the checkpoint's outputs exist where downstream rules expect them. Nothing
but the TSVs is written back to Lustre.

Note that `SLURM_TMPDIR` is per-node and per-job: a retry starts from an empty
tmpdir, and the copy-back happens only after a successful `build-index`, so a
timeout still leaves nothing behind and Snakemake's `--rerun-incomplete`
restarts the whole step.

## Do not split the panel by scaffold

K-mer uniqueness is a genome-wide property of the assembly: a k-mer that
looks unique within one scaffold may be repeated on another. Running
`build-index` per scaffold would change which k-mers survive the uniqueness
filter, so the outputs would not be equivalent to a single run. The panel is
processed whole.

## Memory on `savio4_htc`

Memory is derived by Slurm, not requested: `DefMemPerCPU=4000` MB, so a job
gets `cpus_per_task × 4 GB`. The profile gives this rule
`cpus_per_task: 16` → **64 GB**, and separately `set-threads: 8` for
`--jellyfish-threads`. Do not lower `cpus_per_task` to "match" the threads:
8 CPUs would mean a 32 GB allocation, below the observed peak.

Observed on the full `hap2_vs_hap1` panel (10,522,223 variants → 4,644,674
2-haplotype bubbles, segment FASTA ~2 Gbp, default `--jellyfish-hash 100M`):

| metric | value |
| --- | --- |
| peak RSS | 43.8 GB of 62.5 GB usable |
| jellyfish threads then | 4 |
| wall time then | exceeded the 24 h limit |

The hash size is left at the default on purpose: `-m 31 -s 100M` is ~1 byte
per distinct k-mer per entry and the observed peak fits comfortably inside
the 64 GB allocation, so raising it buys nothing and risks the allocation.
