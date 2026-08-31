# AGENTS.md

## Project overview

Snakemake workflow for *Leptosiphon parviflorus* pool-seq analysis (k-mer counting, read mapping, graph genome construction, diversity statistics, functional annotation).

## Commands

```bash
# Dry-run a specific target rule
snakemake -s workflow/Snakefile --profile workflow/profiles/default -n

# Run the full default target (rule all)
snakemake -s workflow/Snakefile --profile workflow/profiles/default --use-conda

# Run a named target rule
snakemake -s workflow/Snakefile --profile workflow/profiles/default --use-conda vg_diversity_all
```

## Architecture

- **`config/config.yaml`** — all workflow parameters (k-mer sizes, BWA settings, vg/grenedalf thresholds, etc.)
- **`config/samples.tsv`** — sample metadata (sample_id, cram_path, group, sample_size)
- **`workflow/Snakefile`** — top-level Snakefile; parses config and samples, defines target rules, includes rule files
- **`workflow/rules/*.smk`** — 11 rule modules: `kmc.smk`, `panther.smk`, `orthofinder.smk`, `featurecounts.smk`, `degenotate.smk`, `vg_diversity.smk`, `dnds.smk`, `sourmash.smk`, `msmc2.smk`, `svim.smk`, `braker.smk`
- **`workflow/envs/*.yaml`** — conda environment specs for the tools used by rules (18 active envs; bcalm, metaspades, and bbtools were removed as unused)
- **`workflow/profiles/default/config.yaml`** — Slurm executor config for UC Berkeley Savio cluster

## Key constraints

- **This is not a Python package.** There is no `setup.py`, `pyproject.toml`, test suite, linting, or CI.
- **Runs only on Savio HPC.** The profile uses `executor: slurm` with `co_moilab` account and `savio4_htc` partition. Do not attempt to run the full workflow locally.
- **External data is not in the repo.** The `resources/` directory (reference genomes, annotations, CRAM files) is gitignored and must be provisioned separately. Paths in `config.yaml` use `../resources/` relative to `workflow/`.
- **`pool_sizes.csv` is auto-generated** at Snakefile parse time (lines 19-22 of `Snakefile`) from `samples.tsv`. Do not commit or manually edit it.
- **Conda environments are required.** Always use `--use-conda` when running snakemake. BRAKER rules additionally require `--use-singularity` (snakemake auto-pulls `braker_container` from Docker Hub into `.snakemake/singularity/`).
- **BRAKER singularity caveat:** Savio compute nodes have no internet. The first-time container pull must happen on a login node (either via a snakemake run from the login node, or `apptainer pull docker://teambraker/braker3:v3.1.1`); the cached image in `.snakemake/singularity/` is then readable by compute nodes. `sra_download` also needs internet, so run that rule locally with `--local-cores`.
- **Target rules** (entrypoints): `all`, `panther_all`, `orthofinder_all`, `featurecounts_all`, `degenotate_all`, `vg_diversity_all`, `dnds_all`, `sourmash_all`, `msmc2_all`, `svim_all`, `braker_all`.
- **The `.gitignore` is aggressive** — it excludes most bioinformatics file types (`.fastq`, `.bam`, `.fasta`, `.vcf`, `.tsv`, etc.). Be careful when adding new output types.

## Reverted approaches (do not re-add)

When a tool or subworkflow is reverted, record it here with the reason so it isn't re-added later, and commit the reversion with a `chore:` prefix referencing this section.

- **AnchorWave gene-anchored whole-genome alignment** (added 2026-08-31 in 2d75748, switched to prefixed assemblies in e9cedff, removed in [revert sha] as `chore: remove anchorwave subworkflow`) — **structurally incompatible with this project's hap1/hap2 assemblies; do not re-add.** All AnchorWave v1.3.1 anchor-generation modes (`genoAli`, `minimap2`-based `ali`, and `proali`/Quota-alignment) require the reference and query genomes to share contig names: anchor candidates are only kept when `refChr == queryChr` and the contig name exists in the reference GFF, the reference FASTA, and the query FASTA simultaneously (`src/service/TransferGffWithNucmerResult.cpp`; contigs missing from any of the three name spaces are dropped with "There is not enough anchors found on ..."). With zero surviving anchors the run aborts with "there is no match anchor found in the input sam file". Our hap1/hap2 scaffolds are numbered independently and ordered by length, not by homology — measured directly from the splice SAMs, 0 of 6,098 CDS with primary alignments in both haplotypes landed on same-numbered scaffolds (0.0% name correspondence), so no anchors can ever be constructed, with either raw (`scaffold_N`) or prefixed (`hapN_scaffold_N`) names. Do not re-add anchorwave unless it gains all-by-all scaffold-homology support. For synteny/breakpoint analysis use the existing whole-genome minimap2 alignment instead (`minimap2_asm_paf` in `vg_diversity.smk` → `results/vg/{variant}_vs_{backbone}.paf`). (Related earlier removals: bcalm, metaspades, bbtools conda envs were unused.)
