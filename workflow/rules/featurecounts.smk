# featureCounts quantification and filtering of expressed genes for each haplotype.

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


# Expression filtering of BRAKER gene models: count reads per transcript
# over the braker.gff3 annotations with the same RNA-seq BAMs BRAKER used,
# then drop gene models with zero counts and their offspring from the GFF3
# and the protein FASTA. Adapts the featurecounts subworkflow above;
# transcript IDs are also emitted at gene level (stripped .tN suffix) so
# the top-level gene features survive the GFF3 filter.

rule braker_featurecounts_run:
    input:
        gff="results/braker/{ref}/braker.gff3",
        bams=lambda wildcards: expand(
            "results/hisat2/{ref}/{acc}.bam",
            ref=wildcards.ref,
            acc=config["rna_sra_accessions"],
        )
    output:
        "results/braker/{ref}/braker.counts"
    conda: "../envs/featurecounts.yaml"
    params:
        gff_feature=config["featurecounts_gff_feature"],
        paired="-p" if config["featurecounts_paired"] else ""
    shell:
        """
        mkdir -p $(dirname {output})
        featureCounts -a {input.gff} \
            -o {output} \
            -T {threads} \
            {params.paired} \
            -g {params.gff_feature} \
            {input.bams}
        """


rule braker_featurecounts_filter_expressed:
    input:
        "results/braker/{ref}/braker.counts"
    output:
        "results/braker/{ref}/braker_expressed_genes.txt"
    conda: "../envs/featurecounts.yaml"
    params:
        # last column is length; counts are the remaining columns
        n_bams=lambda wildcards: len(config["rna_sra_accessions"])
    shell:
        """
        awk -F'\t' 'NR > 2 && /^g/ {{
            sum = 0
            for (i = NF - {params.n_bams} + 1; i <= NF; i++) sum += $i
            if (sum > 0) {{
                print $1
                sub(/\\.t[0-9]+$/, "", $1)
                print $1
            }}
        }}' {input} | sort -u > {output}
        """


rule braker_featurecounts_filter_gff:
    input:
        glist="results/braker/{ref}/braker_expressed_genes.txt",
        gff="results/braker/{ref}/braker.gff3"
    output:
        "results/braker/{ref}/braker_expressed.gff3"
    conda: "../envs/featurecounts.yaml"
    shell:
        """
        awk -F'\\t' 'NR==FNR{{genes[$1]=1; next}}
            /^#/{{print; next}}
            {{
                match($9, /ID=([^;]+)/, id); match($9, /Parent=([^;]+)/, parent)
                if (id[1] in genes || parent[1] in genes) print
            }}' {input.glist} {input.gff} > {output}
        """


rule braker_filter_expressed_faa:
    input:
        genes="results/braker/{ref}/braker_expressed_genes.txt",
        faa="results/braker/{ref}/braker.aa"
    output:
        "results/braker/{ref}/braker_expressed.faa"
    conda: "../envs/featurecounts.yaml"
    shell:
        """
        # anchor patterns to full header match so g1 does not match g10.t1
        sed 's/.*/^&$/' {input.genes} > {input.genes}.patterns
        seqkit grep -n -r -f {input.genes}.patterns {input.faa} > {output}
        rm {input.genes}.patterns
        """
