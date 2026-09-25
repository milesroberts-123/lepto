# AnchorWave alignment of a single scaffold pair: extract the backbone
# scaffold (and its annotation) and the variant scaffold from their
# assemblies, rename both to a shared name so AnchorWave's matching-name
# anchor requirement is satisfied, splice-align backbone-derived CDS to
# both scaffolds, then genoAli builds the pair alignment.
# One-off pair defined in config["anchorwave_pair"]; fixed paths, no
# wildcards.


rule anchorwave_pair_extract_backbone:
    input:
        fasta=config["anchorwave_pair"]["backbone"]["fasta"],
        gff=config["anchorwave_pair"]["backbone"]["gff"]
    params:
        scaffold=config["anchorwave_pair"]["backbone"]["scaffold"],
        common_name=config["anchorwave_pair_common_name"]
    output:
        fasta="results/anchorwave/pair/backbone.fasta",
        gff="results/anchorwave/pair/backbone.gff3"
    conda: "../envs/featurecounts.yaml"
    shell:
        """
        seqkit grep -p {params.scaffold} {input.fasta} \
            | sed 's/^>.*/>{params.common_name}/' > {output.fasta}
        awk -v OFS='\\t' -v scaf="{params.scaffold}" -v name="{params.common_name}" \
            '$1 == scaf {{ $1 = name; print }}' {input.gff} > {output.gff}
        """


rule anchorwave_pair_extract_variant:
    input:
        fasta=config["anchorwave_pair"]["variant"]["fasta"]
    params:
        scaffold=config["anchorwave_pair"]["variant"]["scaffold"],
        common_name=config["anchorwave_pair_common_name"]
    output:
        "results/anchorwave/pair/variant.fasta"
    conda: "../envs/featurecounts.yaml"
    shell:
        """
        seqkit grep -p {params.scaffold} {input.fasta} \
            | sed 's/^>.*/>{params.common_name}/' > {output}
        """


rule anchorwave_pair_gff2seq:
    input:
        gff="results/anchorwave/pair/backbone.gff3",
        ref="results/anchorwave/pair/backbone.fasta"
    output:
        "results/anchorwave/pair/cds.fa"
    conda: "../envs/anchorwave.yaml"
    shell:
        "anchorwave gff2seq -i {input.gff} -r {input.ref} -o {output}"


rule anchorwave_pair_splice_backbone:
    input:
        fasta="results/anchorwave/pair/backbone.fasta",
        cds="results/anchorwave/pair/cds.fa"
    output:
        "results/anchorwave/pair/backbone_cds.sam"
    threads: 8
    conda: "../envs/minimap2.yaml"
    params:
        k=config["anchorwave_minimap_k"],
        p=config["anchorwave_minimap_p"],
        n=config["anchorwave_minimap_n"]
    shell:
        "minimap2 -x splice -t {threads} -k {params.k} -a -p {params.p} -N {params.n} {input.fasta} {input.cds} > {output}"


rule anchorwave_pair_splice_variant:
    input:
        fasta="results/anchorwave/pair/variant.fasta",
        cds="results/anchorwave/pair/cds.fa"
    output:
        "results/anchorwave/pair/variant_cds.sam"
    threads: 8
    conda: "../envs/minimap2.yaml"
    params:
        k=config["anchorwave_minimap_k"],
        p=config["anchorwave_minimap_p"],
        n=config["anchorwave_minimap_n"]
    shell:
        "minimap2 -x splice -t {threads} -k {params.k} -a -p {params.p} -N {params.n} {input.fasta} {input.cds} > {output}"


rule anchorwave_pair_genoali:
    input:
        gff="results/anchorwave/pair/backbone.gff3",
        cds="results/anchorwave/pair/cds.fa",
        ref="results/anchorwave/pair/backbone.fasta",
        variant="results/anchorwave/pair/variant.fasta",
        ref_sam="results/anchorwave/pair/backbone_cds.sam",
        variant_sam="results/anchorwave/pair/variant_cds.sam"
    output:
        anchors="results/anchorwave/pair/anchors",
        maf="results/anchorwave/pair/maf",
        fmaf="results/anchorwave/pair/f.maf"
    threads: 8
    log:
        "results/anchorwave/pair/genoAli.log"
    conda: "../envs/anchorwave.yaml"
    shell:
        "anchorwave genoAli -t {threads} -i {input.gff} -as {input.cds} -r {input.ref} -a {input.variant_sam} -ar {input.ref_sam} -s {input.variant} -n {output.anchors} -o {output.maf} -f {output.fmaf} > {log} 2>&1"