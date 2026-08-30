# AnchorWave gene-anchored whole-genome alignment: extract CDS from the
# backbone's expressed-gene GFF3 (featurecounts branch), splice-align the
# CDS to both haplotypes with minimap2, then AnchorWave genoAli constructs
# the whole-genome alignment anchored on conserved coding sequence.


rule anchorwave_gff2seq:
    input:
        gff="results/featurecounts/{ref}_expressed.gff3",
        ref_fasta="../resources/reference/{ref}.fasta"
    output:
        "results/anchorwave/{ref}/cds.fa"
    conda: "../envs/anchorwave.yaml"
    shell:
        "anchorwave gff2seq -i {input.gff} -r {input.ref_fasta} -o {output}"


rule minimap2_splice_cds:
    input:
        genome_fasta="../resources/reference/{genome}.fasta",
        cds="results/anchorwave/{backbone}/cds.fa"
    output:
        "results/anchorwave/{variant}_vs_{backbone}/{genome}_cds.sam"
    conda: "../envs/minimap2.yaml"
    params:
        k=config["anchorwave_minimap_k"],
        p=config["anchorwave_minimap_p"],
        n=config["anchorwave_minimap_n"]
    shell:
        "minimap2 -x splice -t {threads} -k {params.k} -a -p {params.p} -N {params.n} {input.genome_fasta} {input.cds} > {output}"


rule anchorwave_genoali:
    input:
        gff="results/featurecounts/{backbone}_expressed.gff3",
        cds="results/anchorwave/{backbone}/cds.fa",
        ref_fasta="../resources/reference/{backbone}.fasta",
        variant_fasta="../resources/reference/{variant}.fasta",
        ref_sam="results/anchorwave/{variant}_vs_{backbone}/{backbone}_cds.sam",
        variant_sam="results/anchorwave/{variant}_vs_{backbone}/{variant}_cds.sam"
    output:
        anchors="results/anchorwave/{variant}_vs_{backbone}/{variant}.anchors",
        maf="results/anchorwave/{variant}_vs_{backbone}/{variant}.maf",
        fmaf="results/anchorwave/{variant}_vs_{backbone}/{variant}.f.maf"
    log:
        "results/anchorwave/{variant}_vs_{backbone}/{variant}.genoAli.log"
    conda: "../envs/anchorwave.yaml"
    shell:
        "anchorwave genoAli -i {input.gff} -as {input.cds} -r {input.ref_fasta} -a {input.variant_sam} -ar {input.ref_sam} -s {input.variant_fasta} -n {output.anchors} -o {output.maf} -f {output.fmaf} > {log} 2>&1"