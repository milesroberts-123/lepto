# MSMC2: pairwise sequentially Markovian coalescent analysis.
# Chain: prefixed backbone fasta -> faidx -> contig list (checkpoint)
#        -> per-chromosome generate_multihetsep -> msmc2 run.
# Input VCF is the variant-vs-backbone VCF from the vg branch;
# callable mask is the degenotate CDS bed for the backbone.


rule msmc2_faidx:
    input:
        "results/vg/{backbone}_prefixed.fasta"
    output:
        "results/vg/{backbone}_prefixed.fasta.fai"
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools faidx {input}"


rule msmc2_mask:
    input:
        "results/degenotate/{backbone}/degeneracy-cds-sites.bed"
    output:
        "results/msmc2/{backbone}_cds.mask.bed.gz"
    conda: "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p results/msmc2
        bgzip -c {input} > {output}
        """


checkpoint msmc2_contigs:
    input:
        "results/vg/{backbone}_prefixed.fasta.fai"
    output:
        "results/msmc2/{backbone}_contigs.txt"
    shell:
        """
        mkdir -p results/msmc2
        cut -f1 {input} > {output}
        """


rule msmc2_filter_snps:
    input:
        vcfgz="results/vg/{variant}_vs_{backbone}.vcf.gz",
        tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi"
    output:
        temp("results/msmc2/{variant}_vs_{backbone}.snps.vcf.gz")
    conda: "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p results/msmc2
        bcftools view -v snps -m2 -M2 {input.vcfgz} -Oz -o {output}
        """


rule msmc2_het_vcf:
    input:
        "results/msmc2/{variant}_vs_{backbone}.snps.vcf.gz"
    output:
        vcfgz=temp("results/msmc2/{variant}_vs_{backbone}.het.vcf.gz"),
        tbi=temp("results/msmc2/{variant}_vs_{backbone}.het.vcf.gz.tbi")
    conda: "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p results/msmc2
        zcat {input} | sed 's:1/1:0/1:g' | bgzip -c > {output.vcfgz}
        tabix -p vcf {output.vcfgz}
        """


rule msmc2_subset_vcf:
    group: "msmc2_subset"
    input:
        vcfgz="results/msmc2/{variant}_vs_{backbone}.het.vcf.gz",
        tbi="results/msmc2/{variant}_vs_{backbone}.het.vcf.gz.tbi"
    output:
        temp("results/msmc2/{variant}_vs_{backbone}/{chrom}.vcf.gz")
    conda: "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p results/msmc2/{wildcards.variant}_vs_{wildcards.backbone}
        bcftools view -r {wildcards.chrom} {input.vcfgz} -Oz -o {output}
        """


rule msmc2_subset_mask:
    group: "msmc2_subset"
    input:
        "results/msmc2/{backbone}_cds.mask.bed.gz"
    output:
        temp("results/msmc2/{backbone}/{chrom}.cds.mask.bed.gz")
    conda: "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p results/msmc2/{wildcards.backbone}
        zcat {input} | awk -v c={wildcards.chrom} '$1==c' | bgzip -c > {output}
        """


rule msmc2_generate:
    group: "msmc2_subset"
    input:
        vcf="results/msmc2/{variant}_vs_{backbone}/{chrom}.vcf.gz",
        mask="results/msmc2/{backbone}/{chrom}.cds.mask.bed.gz"
    output:
        "results/msmc2/{variant}_vs_{backbone}/{chrom}.multihetsep.txt"
    conda: "../envs/msmc2.yaml"
    params:
        script=config["msmc2_generate_multihetsep"]
    shell:
        """
        mkdir -p results/msmc2/{wildcards.variant}_vs_{wildcards.backbone}
        python {params.script} --chr {wildcards.chrom} --mask {input.mask} {input.vcf} > {output}
        """


def msmc2_inputs(wildcards):
    """Return per-chromosome multihetsep files listed in the contig checkpoint."""
    contigs = checkpoints.msmc2_contigs.get(backbone=wildcards.backbone).output[0]
    with open(contigs) as f:
        chroms = [line.strip() for line in f if line.strip()]
    return expand(
        "results/msmc2/{variant}_vs_{backbone}/{chrom}.multihetsep.txt",
        variant=wildcards.variant,
        backbone=wildcards.backbone,
        chrom=chroms,
    )


rule msmc2_run:
    input:
        msmc2_inputs
    output:
        "results/msmc2/{variant}_vs_{backbone}/msmc2.final.txt"
    params:
        binary=config["msmc2_binary"],
        p=config["msmc2_p"],
        prefix=lambda wildcards: f"results/msmc2/{wildcards.variant}_vs_{wildcards.backbone}/msmc2"
    shell:
        """
        {params.binary} -I 0-1 -t {threads} -p {params.p} -o {params.prefix} {input}
        """
