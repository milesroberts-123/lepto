# MCScanX collinearity between the hap1 and hap2 assemblies.
# Chain: featureCounts-expressed proteomes/GFFs -> t1-filtered, species-
# prefixed sequences and gene-level 4-column GFF -> combined all-vs-all
# BLASTP -> MCScanX. BLAST hit names equal GFF gene names by construction
# (hapN.gM), since both sides come from the same expressed-gene lists.


rule mcscanx_take_t1_prefix_faa:
    # Keep one transcript per gene (t1, same representative rule as
    # orthofinder.smk), strip the .tN suffix, prefix with the species so
    # gene names are unique across haplotypes.
    input:
        lambda wildcards: config["mcscanx_input_fastas"][wildcards.species]
    output:
        "results/mcscanx/prefixed/{species}.faa"
    conda: "../envs/mcscanx.yaml"
    shell:
        """
        mkdir -p results/mcscanx/prefixed
        seqkit grep -n -r -p "\\.t1$" {input} \\
            | seqkit replace -p '\\.t1$' -r '' \\
            | seqkit replace -p '^' -r '{wildcards.species}.' > {output}
        """


rule mcscanx_gff_to_4col:
    # Reduce the BRAKER GFF3 to MCScanX's 4-column format (gene, chrom,
    # start, end), keeping only gene-level rows and prefixing both the
    # gene ID and the scaffold name. Scaffold numbering collides between
    # haplotypes, so the chromosome prefix is mandatory.
    input:
        lambda wildcards: config["mcscanx_input_gffs"][wildcards.species]
    output:
        "results/mcscanx/prefixed/{species}.gff"
    shell:
        """
        mkdir -p results/mcscanx/prefixed
        awk -F'\\t' -v OFS='\\t' -v sp='{wildcards.species}' \
            '$3 == "gene" {{
                match($9, /ID=([^;]+)/, id)
                print sp"."id[1], sp"."$1, $4, $5
            }}' {input} > {output}
        """


rule mcscanx_concat_faa:
    input:
        expand("results/mcscanx/prefixed/{species}.faa",
               species=list(config["mcscanx_input_fastas"].keys()))
    output:
        "results/mcscanx/leptosiphon_all.faa"
    shell:
        "cat {input} > {output}"


rule mcscanx_concat_gffs:
    input:
        expand("results/mcscanx/prefixed/{species}.gff",
               species=list(config["mcscanx_input_gffs"].keys()))
    output:
        "results/mcscanx/leptosiphon_all.gff"
    shell:
        "cat {input} > {output}"


rule mcscanx_makeblastdb:
    input:
        "results/mcscanx/leptosiphon_all.faa"
    output:
        # representative makeblastdb product; .pin/.psq/.pog etc. land
        # beside it
        "results/mcscanx/db/leptosiphon.phr"
    conda: "../envs/blast.yaml"
    shell:
        """
        mkdir -p results/mcscanx/db
        rm -f results/mcscanx/db/leptosiphon.*
        makeblastdb -in {input} -dbtype prot -out results/mcscanx/db/leptosiphon
        """


rule mcscanx_blastp:
    # All-vs-all BLASTP on the combined proteome: detects hap1-hap2
    # collinearity as well as within-haplotype synteny.
    input:
        faa="results/mcscanx/leptosiphon_all.faa",
        db="results/mcscanx/db/leptosiphon.phr"
    output:
        "results/mcscanx/leptosiphon_all.hom"
    conda: "../envs/blast.yaml"
    params:
        evalue=config["mcscanx_blast_evalue"],
        max_targets=config["mcscanx_blast_max_targets"]
    shell:
        """
        blastp -query {input.faa} \
            -db results/mcscanx/db/leptosiphon \
            -outfmt 6 \
            -evalue {params.evalue} \
            -max_target_seqs {params.max_targets} \
            -num_threads {threads} \
            -out {output}
        """


rule mcscanx_run:
    # MCScanX wants the .hom and .gff files sharing a basename in its
    # working directory; it writes .collinearity, .tandem and .html/
    # next to its input prefix.
    input:
        hom="results/mcscanx/leptosiphon_all.hom",
        gff="results/mcscanx/leptosiphon_all.gff"
    output:
        collinearity="results/mcscanx/collinearity/leptosiphon_all.collinearity",
        tandem="results/mcscanx/collinearity/leptosiphon_all.tandem",
        html=directory("results/mcscanx/collinearity/leptosiphon_all.html")
    conda: "../envs/mcscanx.yaml"
    params:
        workdir="results/mcscanx/collinearity",
        prefix="leptosiphon_all"
    shell:
        """
        mkdir -p {params.workdir}
        cp {input.hom} {params.workdir}/{params.prefix}.hom
        cp {input.gff} {params.workdir}/{params.prefix}.gff
        cd {params.workdir} && MCScanX {params.prefix}
        """