rule featurecounts_run:
    input:
        gff=lambda wildcards: config["featurecounts_annotations"][wildcards.species],
        bams=lambda wildcards: config["featurecounts_bams"][wildcards.species]
    output:
        "results/featurecounts/{species}.counts"
    conda: "../envs/featurecounts.yaml"
    params:
        gff_feature=config["featurecounts_gff_feature"],
        paired="-p" if config["featurecounts_paired"] else ""
    shell:
        """
        mkdir -p results/featurecounts
        featureCounts -a {input.gff} \
            -o {output} \
            {params.paired} \
            -g {params.gff_feature} \
            {input.bams}
        """


rule featurecounts_filter_expressed:
    input:
        "results/featurecounts/{species}.counts"
    output:
        "results/featurecounts/{species}_expressed_genes.txt"
    shell:
        """
        grep -v -P "\t0\t0\t0\t0\t0\t0$" {input} | cut -f 1 | grep "^g" > {output}
        """


rule filter_expressed_fasta:
    input:
        genes="results/featurecounts/{species}_expressed_genes.txt",
        fasta=lambda wildcards: config["panther_input_fastas"][wildcards.species]
    output:
        "results/featurecounts/{species}_expressed.faa"
    conda: "../envs/featurecounts.yaml"
    shell:
        """
        sed 's/^/^/; s/$/\\.t1\\b/' {input.genes} > results/featurecounts/{wildcards.species}_patterns.txt
        seqkit grep -n -r -f results/featurecounts/{wildcards.species}_patterns.txt {input.fasta} > {output}
        """
