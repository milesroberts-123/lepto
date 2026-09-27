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
- **`workflow/rules/*.smk`** — 13 rule modules: `kmc.smk`, `panther.smk`, `orthofinder.smk`, `featurecounts.smk`, `degenotate.smk`, `vg_diversity.smk`, `dnds.smk`, `sourmash.smk`, `msmc2.smk`, `svim.smk`, `braker.smk`, `freqk.smk`, `anchorwave_pair.smk`
- **`workflow/envs/*.yaml`** — conda environment specs for the tools used by rules (20 active envs; bcalm, metaspades, and bbtools were removed as unused)
- **`workflow/profiles/default/config.yaml`** — Slurm executor config for UC Berkeley Savio cluster

## Key constraints

- **This is not a Python package.** There is no `setup.py`, `pyproject.toml`, test suite, linting, or CI.
- **Runs only on Savio HPC.** The profile uses `executor: slurm` with `co_moilab` account and `savio4_htc` partition. Do not attempt to run the full workflow locally.
- **External data is not in the repo.** The `resources/` directory (reference genomes, annotations, CRAM files) is gitignored and must be provisioned separately. Paths in `config.yaml` use `../resources/` relative to `workflow/`.
- **`pool_sizes.csv` is auto-generated** at Snakefile parse time (lines 19-22 of `Snakefile`) from `samples.tsv`. Do not commit or manually edit it.
- **Conda environments are required.** Always use `--use-conda` when running snakemake. BRAKER rules additionally require `--use-singularity` (snakemake auto-pulls `braker_container` from Docker Hub into `.snakemake/singularity/`). The `freqk` binary is NOT conda-managed: it must be placed manually at the path in `config["freqk_binary"]` (`/global/scratch/users/milesroberts/brandvain_lab_projects/lepto/bin/freqk`; grab the Linux release from github.com/milesroberts-123/freqk), like the `msmc2_binary` precedent.
- **BRAKER singularity caveat:** Savio compute nodes have no internet. The first-time container pull must happen on a login node (either via a snakemake run from the login node, or `apptainer pull docker://teambraker/braker3:v3.1.1`); the cached image in `.snakemake/singularity/` is then readable by compute nodes. `sra_download` and `odb_download` (OrthoDB Viridiplantae protein evidence, gated by `config["braker_use_protein"]`) also need internet, so run those rules on a login node with `--local-cores`.
- **Target rules** (entrypoints): `all`, `panther_all`, `orthofinder_all`, `featurecounts_all`, `degenotate_all`, `vg_diversity_all`, `dnds_all`, `sourmash_all`, `msmc2_all`, `svim_all`, `braker_all`, `freqk_all`, `anchorwave_pair_all`.
- **The `.gitignore` is aggressive** — it excludes most bioinformatics file types (`.fastq`, `.bam`, `.fasta`, `.vcf`, `.tsv`, etc.). Be careful when adding new output types.
- **Session transcripts must be archived before every push.** Run `scripts/export_session.sh` first: it exports the full transcript of the current opencode session as JSON (`opencode export <session-id>`), compresses it with `xz -9` into `transcripts/<session-id>.json.xz` (one file per session, overwritten on each push), and stages it. Commit the staged archive with a `chore:` prefix before pushing; if the script fails, do not push.

## Reverted approaches (do not re-add)

When a tool or subworkflow is reverted, record it here with the reason so it isn't re-added later, and commit the reversion with a `chore:` prefix referencing this section.

- **AnchorWave gene-anchored whole-genome alignment** (added 2026-08-31 in 2d75748, switched to prefixed assemblies in e9cedff, removed in 84fa7ec as `chore: remove anchorwave subworkflow`) — **structurally incompatible with this project's hap1/hap2 assemblies; do not re-add.** All AnchorWave v1.3.1 anchor-generation modes (`genoAli`, `minimap2`-based `ali`, and `proali`/Quota-alignment) require the reference and query genomes to share contig names: anchor candidates are only kept when `refChr == queryChr` and the contig name exists in the reference GFF, the reference FASTA, and the query FASTA simultaneously (`src/service/TransferGffWithNucmerResult.cpp`; contigs missing from any of the three name spaces are dropped with "There is not enough anchors found on ..."). With zero surviving anchors the run aborts with "there is no match anchor found in the input sam file". Our hap1/hap2 scaffolds are numbered independently and ordered by length, not by homology — measured directly from the splice SAMs, 0 of 6,098 CDS with primary alignments in both haplotypes landed on same-numbered scaffolds (0.0% name correspondence), so no anchors can ever be constructed, with either raw (`scaffold_N`) or prefixed (`hapN_scaffold_N`) names. Do not re-add anchorwave unless it gains all-by-all scaffold-homology support. For synteny/breakpoint analysis use the existing whole-genome minimap2 alignment instead (`minimap2_asm_paf` in `vg_diversity.smk` → `results/vg/{variant}_vs_{backbone}.paf`). (Related earlier removals: bcalm, metaspades, bbtools conda envs were unused.) **Sanctioned exception (added 2026-09-25, `anchorwave_pair_all` in `anchorwave_pair.smk`):** single scaffold-pair alignments where both scaffolds are explicitly extracted and renamed to a shared name (`anchorwave_pair_common_name`, e.g. `scaffold_pair`) before `genoAli`, satisfying the matching-name requirement by construction. The name mismatch that killed the whole-genome branch cannot occur; keep this pattern (extract → rename → gff2seq → splice SAMs → genoAli) if another pair is needed.
- **freqk grenedalf polymorphic-site exclusion** (added 2026-09-21 in 61499f5, removed 2026-09-22, re-added 2026-09-26 for W-group pools only) — **couples freqk to the vg mapping chain; keep it scoped to W-group samples.** The original `freqk_exclusion_bed` consumed all 8 per-pool `grenedalf_results_frequency.csv` tables to build a BED of sites polymorphic (0.02 < freq < 0.98) in any pool, which `freqk_filter_vcf` excluded with `bcftools view -T ^bed`; it was reverted because adding those CSVs as freqk inputs made any freqk-only target (`snkc -n freqk_all`) schedule the entire vg chain (74+ jobs: samtools_fastq, fastp, vg_giraffe, grenedalf, …) — consumed temp files force upstream rebuilds into the consumer's DAG (verified: `--rulegraph`/`--dag` showed grenedalf→freqk edges; `ancient()` on the CSVs and on `results/fastp/{ID}.fastq` did NOT break the cascade — snakemake 9.16.2 only honors ancient when all such inputs both are ancient AND exist, so missing temps drag consumers regardless; dag.py:1609-1612). The current sanctioned design (2026-09-26) re-adds `freqk_exclusion_bed` but consumes only the 4 W-group (`group == "W"`, `w_sample_ids` in Snakefile) frequency tables, with thresholds as config keys (`freqk_min_alt_count: 1`, `freqk_min_depth: 10`, `freqk_min_polymorphic_freq: 0.02`, `freqk_max_polymorphic_freq: 0.98`; awk `$1 != "CHROM" && $6 > c && $7 >= d && $8 > f && $8 < x`), unioned into a 0-based BED (`sort -u -k1,1 -k2,2n`); `freqk_filter_vcf` pipes `bcftools norm -m -any | bcftools view -T ^bed | bcftools view -i 'QUAL>=…' -v snps -m2 -M2`. Run `vg_diversity_all` first and keep its products; freqk-only targets still schedule the vg chain when the CSVs are missing — do not widen the CSV set (back to all 8 pools or other groups) without re-accepting that tradeoff.
