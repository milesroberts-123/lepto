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
- **`workflow/rules/*.smk`** — 20 rule modules: `kmc.smk`, `panther.smk`, `orthofinder.smk`, `featurecounts.smk`, `degenotate.smk`, `vg_diversity.smk`, `dnds.smk`, `sourmash.smk`, `msmc2.smk`, `svim.smk`, `braker.smk`, `freqk.smk`, `anchorwave_pair.smk`, `vigs.smk`, `mcscanx.smk`, `kmate.smk`, `cds_windows.smk`, `ldx.smk`, `busco.smk`, `tubemap.smk`
- **`workflow/envs/*.yaml`** — conda environment specs for the tools used by rules (26 active envs; bcalm, metaspades, and bbtools were removed as unused)
- **`workflow/profiles/default/config.yaml`** — Slurm executor config for UC Berkeley Savio cluster

## Key constraints

- **This is not a Python package.** There is no `setup.py`, `pyproject.toml`, test suite, linting, or CI.
- **Runs only on Savio HPC.** The profile uses `executor: slurm` with `co_moilab` account and `savio4_htc` partition. Do not attempt to run the full workflow locally. On htc partitions memory must NOT be requested explicitly: Slurm derives job memory from cores at `DefMemPerCPU=4000` MB/core (verified via `scontrol show partition savio4_htc`; 56-core/503 GiB nodes). Size `cpus_per_task` for memory-hungry rules (e.g. `kmate_kmer_db`/`kmate_run` 12 cores -> 48 GB covers the 32G jellyfish hash; with 8 cores the ~32 GiB table OOM-killed the step in 25 s — 2026-10-22). The only sanctioned memory-like key is `vg_autoindex`'s `max_mem: "32G"`, which vg autoindex consumes itself (`-M`), not a `--mem` request.
- **External data is not in the repo.** The `resources/` directory (reference genomes, annotations, CRAM files) is gitignored and must be provisioned separately. Paths in `config.yaml` use `../resources/` relative to `workflow/`.
- **`pool_sizes.csv` is auto-generated** at Snakefile parse time (lines 19-22 of `Snakefile`) from `samples.tsv`. Do not commit or manually edit it.
- **Conda environments are required.** Always use `--use-conda` when running snakemake. BRAKER rules additionally require `--use-singularity` (snakemake auto-pulls `braker_container` from Docker Hub into `.snakemake/singularity/`). The `freqk` binary is NOT conda-managed: it must be placed manually at the path in `config["freqk_binary"]` (`/global/scratch/users/milesroberts/brandvain_lab_projects/lepto/bin/freqk`; grab the Linux release from github.com/milesroberts-123/freqk), like the `msmc2_binary` precedent.
- **BRAKER singularity caveat:** Savio compute nodes have no internet. The first-time container pull must happen on a login node (either via a snakemake run from the login node, or `apptainer pull docker://teambraker/braker3:v3.1.1`); the cached image in `.snakemake/singularity/` is then readable by compute nodes. `sra_download` and `odb_download` (OrthoDB Viridiplantae protein evidence, gated by `config["braker_use_protein"]`) also need internet, so run those rules on a login node with `--local-cores`. Same for `busco_lineage_download` (see `docs/workflows/busco.md`).
- **Target rules** (entrypoints): `all`, `panther_all`, `orthofinder_all`, `featurecounts_all`, `degenotate_all`, `vg_diversity_all`, `dnds_all`, `sourmash_all`, `msmc2_all`, `svim_all`, `braker_all`, `freqk_all`, `anchorwave_pair_all`, `vigs_all`, `mcscanx_all`, `kmate_all`, `cds_windows_all`, `ldx_all`, `busco_all`, `vg_tubemap_all`, `mk_all`.
- **The `.gitignore` is aggressive** — it excludes most bioinformatics file types (`.fastq`, `.bam`, `.fasta`, `.vcf`, `.tsv`, etc.). Be careful when adding new output types.
- **Session transcripts must be archived before every push.** Run `scripts/export_session.sh` first: it exports the full transcript of the current opencode session as JSON (`opencode export <session-id>`), compresses it with `xz -9` into `transcripts/<session-id>.json.xz` (one file per session, overwritten on each push), and stages it. Commit the staged archive with a `chore:` prefix before pushing; if the script fails, do not push.

## Where the details live

Do not re-derive these from code archaeology. Read the matching file before editing, and update it when behaviour changes.

| Area | File |
| --- | --- |
| BUSCO chain, lineage download | `docs/workflows/busco.md` |
| Sequence tube map (`tubemap.smk`) | `docs/workflows/tubemap.md` |
| BRAKER3 annotation, AUGUSTUS config rebuild, masking checks | `docs/workflows/braker.md` |
| McDonald-Kreitman chain, `mk_polymorphism` memory design | `docs/workflows/mk.md` |
| Grenedalf interval outputs, filename geometry | `docs/workflows/grenedalf.md` |
| Approaches that were rejected — why, and what not to re-add | `docs/reverted-approaches.md` |
| Full pre-2026-11-19 AGENTS.md, for anything missing | `docs/legacy/AGENTS-2026-11-19.md` |

Rule: if a change alters a workflow's behaviour, add the detail to the matching doc and keep `AGENTS.md` to commands and hard constraints. Never move operational history back into this file.

## Manuscript submodule

- **`manuscript/` is a git submodule** pointing at `github.com/milesroberts-123/manuscript-leptosiphon-parviflorus` (LaTeX source, synced to Overleaf through GitHub). The analysis repo records an exact manuscript commit, so any analysis commit names the text it belongs to.
- **Commit order for cross-repo changes:** commit and push inside `manuscript/` first, then `git add manuscript` in this repo and commit the updated pointer with a `chore:` prefix. Never `git add -A` while the submodule is dirty — it stages a pointer to a commit nobody else can fetch.
- **Clone with submodules:** `git clone --recurse-submodules <url>`, or `git submodule update --init` in an existing clone. `git worktree add` does not initialise submodules; run `git submodule update --init` in every new lane before touching `manuscript/`.
- **Do not edit the manuscript's git history or push it from here** unless the task is explicitly about the manuscript. Prose edits belong in the manuscript repo; analysis changes here should only bump the pointer when the figure or text contract changed.
- **Detached HEAD is normal** inside `manuscript/`: check out the branch you mean (`git -C manuscript switch main`) before editing, or commit on a new branch and update the pointer to it.
