# Grenedalf interval outputs and naming

Moved out of AGENTS.md on 2026-11-19 to keep per-turn context small. Content is verbatim.

- **Grenedalf interval filenames embed the window geometry** — `grenedalf_diversity_interval`/`grenedalf_diversity_interval_pooled`/`grenedalf_fst_interval` produce `grenedalf_results_interval_w{wlen}_s{sstep}_*.csv` or `grenedalf_results_pooled_w{wlen}_s{sstep}_diversity.csv` where `wlen`/`sstep` wildcards carry `grenedalf_window_width`/`grenedalf_window_interval_stride` straight from config (the rules' `--window-interval-*` flags and `--file-prefix` come from the wildcards, so filename and analysis parameters cannot drift; `wildcard_constraints: wlen/sstep=[0-9]+` in the Snakefile). Changing those config keys yields new filenames and a rebuild; genome-window outputs stay `grenedalf_results_genome_*` (never suffixed). Anything globbing interval outputs by name must use `interval_.*_(diversity|fst).csv` patterns (Rmd line ~1118 patched accordingly). Naming lesson (2026-10-15): the wildcards were first called `w`/`s`, which broke `snkc -n sourmash_all` — global `wildcard_constraints` apply name-wide across ALL subworkflows, and `{s}`/`{t}` are already taken by `sourmash_mutation_rate` (`mutation_rate_{s}_vs_{t}.txt`); the `[0-9]+` constraint made `s=hap1` unmatchable → MissingInputException. Before introducing a new wildcard name, grep every `.smk` + Snakefile for `{name}` collisions.

## Per-bin read depth (`grenedalf_depth_bins`)

Added 2026-11-20. One streaming pass per pool over
`grenedalf_results_frequency.csv` (which has a row per dictionary site because
of `--write-invariants`) produces

```
results/vg/{variant}_vs_{backbone}/{ID}/depth_bins_{bin}bp.tsv
```

with `chrom, bin_start, depth_sum, n_sites`. The bin width comes from
`config["grenedalf_depth_bin"]` and is **baked into the filename, not a
wildcard**: a global `{bin}` would collide with other rules' naming.

Purpose: the CDS-density / diversity correlation can be a mappability
artefact, so the notebook (`depth_vs_cds` chunk) rebuilds 50 kb / 100 kb /
500 kb mean depth from these bins and regresses π on CDS density with depth as
a covariate. Depth here is post-filter — the same base/map-quality thresholds
that produced θπ — which is what makes it the right confounder. Any window
size that is a multiple of the bin width can be rebuilt from the bins without
another CSV pass.

The work happens in `scripts/grenedalf_window_depth.sh <frequency_csv>
<out_tsv> [bin_bp]`; the rule just wires paths. `vg_diversity_all` expands it
for every pool, one 1-core job (~45 min) each, running in parallel.
