# BUSCO protein-mode annotation completeness (odb12 lineage) on the
# BRAKER3 annotations. Chain: lineage dataset download -> one-isoform
# protein selection -> busco protein run per braker_genomes key.
# Protein mode only (transcriptome/genome modes not run). One
# representative transcript per gene (.t1) caps BUSCO duplication from
# alternative isoforms (repo-wide mcscanx/orthofinder convention).
# Output filenames embed the lineage (grenedalf naming lesson):
# switching busco_lineage produces fresh paths rather than mixing
# datasets in one directory. The lineage download needs internet: run
# busco_lineage_download once on a login node (odb_download precedent;
# Savio compute nodes have none).

busco_lineage = config["busco_lineage"]
busco_isoform = config["busco_isoform"]
busco_downloads_dir = config["busco_downloads_dir"]

busco_proteins_select_faa = (
    f"results/busco/{{ref}}/busco_proteins_select_{busco_isoform}.faa"
)
busco_summary_txt = (
    f"results/busco/{{ref}}/busco_proteins/"
    f"short_summary.specific.{busco_lineage}.busco_proteins.txt"
)
# genome mode: the ceiling an annotation is judged against (same lineage,
# same BUSCO version; the mode suffix keeps the two runs from colliding)
busco_genome_summary_txt = (
    f"results/busco/{{ref}}/busco_genome/"
    f"short_summary.specific.{busco_lineage}.busco_genome.txt"
)


rule busco_lineage_download:
    # Needs internet: run on a login node.
    output:
        directory(f"results/busco/downloads/lineages/{busco_lineage}")
    params:
        lineage=busco_lineage,
        downloads_dir=busco_downloads_dir
    conda: "../envs/busco.yaml"
    shell:
        """
        busco --download {params.lineage} --download_path {params.downloads_dir}
        """


rule busco_proteins_select:
    input:
        faa="results/braker/{ref}/braker.aa"
    output:
        faa=busco_proteins_select_faa
    params:
        # \.t1$ keeps TSEBRA's primary transcript per gene; "none" keeps
        # every transcript (pattern "." matches every bare-ID header).
        pattern=r"\.t1$" if busco_isoform == "t1" else "."
    conda: "../envs/busco.yaml"
    shell:
        """
        seqkit grep -n -r -p '{params.pattern}' {input.faa} > {output.faa}
        """


rule busco_proteins:
    input:
        proteins=busco_proteins_select_faa,
        lineage=directory(f"results/busco/downloads/lineages/{busco_lineage}")
    output:
        summary=busco_summary_txt
    params:
        lineage=busco_lineage,
        downloads_dir=busco_downloads_dir
    conda: "../envs/busco.yaml"
    shell:
        """
        mkdir -p results/busco/{wildcards.ref}
        busco -i {input.proteins} \
            -l {params.lineage} -m proteins -c {threads} \
            --offline --download_path {params.downloads_dir} \
            --out_path results/busco/{wildcards.ref} \
            -o busco_proteins --force
        """


rule busco_genome:
    # Genome-mode BUSCO on the same assembly BRAKER annotated: the ceiling
    # for annotation completeness. full_table.tsv from this run and from
    # busco_proteins can be joined on "# Busco id" to see which conserved
    # genes the assembly has but the annotation lost.
    input:
        genome=lambda wildcards: braker_genomes[wildcards.ref],
        lineage=directory(f"results/busco/downloads/lineages/{busco_lineage}")
    output:
        summary=busco_genome_summary_txt
    params:
        lineage=busco_lineage,
        downloads_dir=busco_downloads_dir
    conda: "../envs/busco.yaml"
    shell:
        """
        mkdir -p results/busco/{wildcards.ref}
        busco -i {input.genome} \
            -l {params.lineage} -m genome -c {threads} \
            --offline --download_path {params.downloads_dir} \
            --out_path results/busco/{wildcards.ref} \
            -o busco_genome --force
        """