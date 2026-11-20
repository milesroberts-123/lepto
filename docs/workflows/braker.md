# BRAKER3 annotation (`braker.smk`)

Target: `braker_all`. Chain: OrthoDB/SRA download (login node, needs
internet) → fastp trimming → HISAT2 `--dta` → `braker.pl` in the
`teambraker/braker3:v3.1.1` container → featureCounts expression filter.

## Running it

The rule declares `container: braker_container`, so **the run needs
singularity/apptainer in addition to conda**:

```bash
snkc braker_all --sdm conda singularity    # or: --sdm conda apptainer
```

The `snkc` alias only enables conda, so without the extra flag snakemake
prints `Singularity containers: ignored`, the image is skipped, and the
rule dies immediately with `AUGUSTUS_CONFIG_PATH: unbound variable`
(the rule now emits an explicit error for this case). The first pull needs
internet, so cache the image from a login node:

```bash
apptainer pull docker://teambraker/braker3:v3.1.1   # into .snakemake/singularity/
```

## AUGUSTUS config is rebuilt on every run

`augustus_config/` under `results/braker/{ref}/` is **not** a cached
artifact. The rule deletes it and re-copies `$AUGUSTUS_CONFIG_PATH` from
the container before calling `braker.pl`; `--species={braker_species}`
then creates the species parameter set and BRAKER's etraining /
`optimize_augustus.pl` writes the trained parameters into it.

Why: BRAKER writes training results into that directory, so a config left
over from an earlier attempt makes the annotation depend on run history.
The failure that motivated this (2026-11-19): an `--skipOptimize` run
left an untrained `linanthus_parryae` parameter set in place, and a later
run with optimization enabled would have "optimized" on top of it instead
of starting from the container default.

Never add an `if [ ! -d ... ]` guard here; a resumed run must retrain.

## Genome masking

`--softmasking` is passed unconditionally. It is only meaningful when the
assembly actually carries lowercase repeat masking; verify before blaming
parameters for poor BUSCO scores:

```bash
grep -v '^>' genome.fasta | tr -d 'ACGTNacgtn\n' | wc -c   # other IUPAC codes
grep -v '^>' genome.fasta | tr -cd 'a-z' | wc -c           # lowercase (soft-masked)
grep -v '^>' genome.fasta | tr -cd 'N' | wc -c             # hard-masked
```

## Evidence

- RNA-seq: `rna_sra_accessions` (5 runs), HISAT2 `--dta`, no dedup or
  multimap filtering. Check `braker.log` for the RNA evidence summary and
  the AUGUSTUS training accuracy line ("accuracy after initial training").
- Protein: OrthoDB v12 Viridiplantae (`braker_odb_url`), toggled by
  `braker_use_protein`; ETP mode runs `gmetp.pl`, which dominates runtime
  (44 h at 32 cores for *L. parryae*).

## Reading a failed or weak run

`braker.log` in `results/braker/{ref}/braker/` records the exact
`braker.pl` call, the ETP and etraining commands, and the AUGUSTUS
training accuracy. A high "error rate of missing stop codon" in the
training set, or an accuracy near 0.64, points at the training gene set
rather than at the AUGUSTUS search itself.
