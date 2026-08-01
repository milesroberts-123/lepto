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
            -T {threads} \
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

rule featurecounts_filter_gff:
    input:
        glist="results/featurecounts/{species}_expressed_genes.txt",
        gff=lambda wildcards: config["featurecounts_annotations"][wildcards.species]
    output:
        "results/featurecounts/{species}_expressed.gff3"
    shell:
        """
        awk -F'\t' 'NR==FNR{{genes[$1]=1; next}}
            /^#/{{print; next}}
            {{
                match($9, /ID=([^;]+)/, id); match($9, /Parent=([^;]+)/, parent)
                if (id[1] in genes || parent[1] in genes) print
            }}' {input.glist} {input.gff} > {output}
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
        seqkit grep -n -r -f {input.genes} {input.fasta} > {output}
        """
