# BUSCO annotation completeness

Moved out of AGENTS.md on 2026-11-19 to keep per-turn context small. Content is verbatim.

- **BUSCO annotation completeness** (`busco.smk`, added 2026-10-21) — protein-mode BUSCO (`busco=6.1.0` + seqkit in `envs/busco.yaml`; 6.1.0 because the user's manual env empirically ran it with the version-suffixed lineage name) on every BRAKER3 annotation (wildcards over `braker_refs` = `braker_genomes` keys; just `linanthus` today). Chain: `busco_lineage_download` (`busco --download` of `config["busco_lineage"]` = `eudicotyledons_odb12.2` into `config["busco_downloads_dir"]`; login node like `odb_download`, confirm spelling with `busco --list-datasets`) → `busco_proteins_select` (seqkit `.t1` filter of `braker.aa` → `results/busco/{ref}/busco_proteins_select_{iso}.faa`; representative-isoform convention, TSEBRA t1 ≈ the 39,769 genes) → `busco_proteins` (`-m proteins -c {threads} --offline`, output `results/busco/{ref}/busco_proteins/short_summary.specific.{lineage}.busco_proteins.txt` — the only snakemake-tracked output; lineage name sits in both the out-path and the summary filename per the grenedalf naming lesson). `busco_isoform: "none"` switches the selection pattern to `.` (all transcripts; D inflated by 2,814 extra isoforms). Target: `busco_all`. Manual BUSCO 6.1.0 pre-runs exist untracked under `workflow/busco_{proteins,longest_proteins,downloads}/` — longest-isoform-on-expressed proteome scored C:71.9%[S:63.6%,D:8.2%],F:1.7%,M:26.4%,n:1990 (eudicotyledons_odb12.2) — the canonical t1-on-full-annotation numbers may differ.

## Genome mode (`busco_genome`)

Added 2026-11-20. `busco_genome` runs the same lineage and BUSCO version in
`-m genome` on the assembly BRAKER annotated (`braker_genomes[ref]`), giving
the ceiling an annotation is judged against:

```
results/busco/{ref}/busco_genome/short_summary.specific.{lineage}.busco_genome.txt
```

It needs only the assembly and the downloaded lineage, so it runs without the
BRAKER chain — useful while the annotation is still being tuned. The profile
allocates 32 cores / 600 min (genome mode is the memory-heavy BUSCO mode).

**Reading it against the protein run.** The two summary percentages are not
strictly comparable: genome mode sees introns and uses homology search while
protein mode scores spliced models, and their Complete/Fragmented definitions
are not nested. Compare per gene instead — both runs write `full_table.tsv`,
joined on `# Busco id`:

- Complete in genome **and** Missing/Fragmented in protein = conserved genes
  the assembly carries but the annotation lost: the recoverable gap.
- Missing in genome mode = outside any annotation's reach; subtract those
  before judging BRAKER.

Clustering the lost genes by masked-region overlap, scaffold-end proximity and
single-exon status is what distinguishes masking, assembly fragmentation and
model prediction as the cause.
