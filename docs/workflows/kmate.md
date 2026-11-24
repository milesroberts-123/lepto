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

The hash size is left at the default (`-s 100M`) because it was out of scope
for this change, not because it is right-sized: a ~2 Gbp segment set has far
more than 100 M distinct k-mers, so jellyfish grows the table and rehashes —
which is why the peak (43.8 GB) is far above the 100 M initial allocation.

## K-mer DB sizing (`kmate_kmer_db`)

`kmate build-kmer-db` is a thin wrapper over

```
zcat -f <reads> | jellyfish count -m 31 -s <hash> -t <threads> -C -o <db>
```

and exposes none of the switches needed to make that tractable, so the rule
calls `jellyfish count` directly with the same binary and the same `-C`, which
keeps the DB readable by `kmate run --kmer-db` unchanged.

Two flags matter:

| flag | why |
| --- | --- |
| `--disk` | when the table fills, write a partial hash and keep going instead of doubling the table in RAM. Peak RSS stays near `-s`; the price is spill files plus a merge pass, both under `SLURM_TMPDIR`. |
| `--lower-count 2` | keeps singleton k-mers out of the written DB. This is an **output** filter — it shrinks the DB that `kmate run` mmaps, it does not reduce peak RSS. |

`-s` is the hash size in entries, ~4.6 B each (`8G` ≈ 37 GB). With `--disk` an
undersized table costs extra spill files, not a resize and not a wrong answer.

**Do not use `--bf-size`.** It enables jellyfish's one-pass Bloom filter, which
looks like the obvious fix and is not:

- the bloom filter must be sized to the *total* k-mer count, ~3×10¹¹ for a
  456 GB pool, i.e. hundreds of GB of filter — larger than the hash it replaces;
- at a smaller size it saturates, so nothing is filtered and the cost is pure
  overhead;
- and the one-pass method reports every surviving count **one lower than
  actual** (`count 2` k-mers come back as 1), which biases the k-mer weights
  that `kmate run` feeds to its EM.

Note also that `-c` is `--counter-len` (bits per counter, default 7), not a
minimum count: `-c 2` would silently cap every count at 3.

`cpus_per_task: 16` for this rule gives 64 GB of derived memory, comfortable
for a ~37 GB table plus spill buffers.
