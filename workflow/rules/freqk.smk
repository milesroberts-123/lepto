# Estimate allele frequencies of known variants in pool-seq reads from
# k-mer counts. Variants come from the vg pipeline's paftools-call VCF on
# the hap1 backbone (run vg_diversity_all first: this branch expects the
# VCF, its tabix index, and the prefixed backbone fasta to already
# exist); reads are the same fastp-cleaned pooled fastqs fed to vg
# giraffe. Before indexing, the VCF is normalized and subset to
# biallelic SNPs. The freqk binary is provisioned outside the repo and
# referenced by the freqk_binary config key (see AGENTS.md).

rule freqk_filter_vcf:
    input:
        "results/vg/{variant}_vs_{backbone}.vcf.gz",
        tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi"
    output:
        vcfgz=temp("results/freqk/{variant}_vs_{backbone}/norm.vcf.gz"),
        tbi=temp("results/freqk/{variant}_vs_{backbone}/norm.vcf.gz.tbi")
    conda: "../envs/bcftools.yaml"
    params:
        freqk_min_qual=config["freqk_min_qual"]
    shell:
        """
        mkdir -p results/freqk/{wildcards.variant}_vs_{wildcards.backbone}
        bcftools norm -m -any {input} \\
            | bcftools view -i 'QUAL>={params.freqk_min_qual}' -v snps -m2 -M2 -Oz -o {output.vcfgz}
        tabix -p vcf {output.vcfgz}
        """

rule freqk_index:
    input:
        fasta="results/vg/{backbone}_prefixed.fasta",
        fai="results/vg/{backbone}_prefixed.fasta.fai",
        #vcf="results/vg/{variant}_vs_{backbone}.vcf.gz",
        #tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi",
        vcf="results/freqk/{variant}_vs_{backbone}/norm.vcf.gz",
        tbi="results/freqk/{variant}_vs_{backbone}/norm.vcf.gz.tbi"
    output:
        temp("results/freqk/{variant}_vs_{backbone}/index.txt")
    benchmark:
        "benchmarks/freqk_index/freqk_index_{variant}_vs_{backbone}.bench"
    params:
        binary=config["freqk_binary"],
        k=config["freqk_k"]
    shell:
        "{params.binary} index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output}"

rule freqk_var_dedup:
    input:
        "results/freqk/{variant}_vs_{backbone}/index.txt"
    output:
        temp("results/freqk/{variant}_vs_{backbone}/var_index.txt")
    benchmark:
        "benchmarks/freqk_index/freqk_var_dedup_{variant}_vs_{backbone}.bench"
    params:
        binary=config["freqk_binary"]
    shell:
        "{params.binary} var-dedup --index {input} --output {output}"

rule freqk_ref_dedup:
    input:
        index="results/freqk/{variant}_vs_{backbone}/var_index.txt",
        fasta="results/vg/{backbone}_prefixed.fasta",
        fai="results/vg/{backbone}_prefixed.fasta.fai",
        #vcf="results/vg/{variant}_vs_{backbone}.vcf.gz",
        #tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi",
        vcf="results/freqk/{variant}_vs_{backbone}/norm.vcf.gz",
        tbi="results/freqk/{variant}_vs_{backbone}/norm.vcf.gz.tbi"
    output:
        "results/freqk/{variant}_vs_{backbone}/ref_index.txt"
    benchmark:
        "benchmarks/freqk_index/freqk_ref_dedup_{variant}_vs_{backbone}.bench"
    params:
        binary=config["freqk_binary"]
    shell:
        "{params.binary} ref-dedup --index {input.index} --fasta {input.fasta} --vcf {input.vcf} --output {output}"

rule freqk_count:
    input:
        reads="results/fastp/{ID}.fastq",
        index="results/freqk/{variant}_vs_{backbone}/ref_index.txt"
    output:
        counts="results/freqk/{variant}_vs_{backbone}/{ID}_counts.txt",
        freqs="results/freqk/{variant}_vs_{backbone}/{ID}_freqs.txt"
    benchmark:
        "benchmarks/freqk_count/freqk_count_{variant}_vs_{backbone}_{ID}.bench"
    params:
        binary=config["freqk_binary"]
    shell:
        "{params.binary} count -vvvv --nthreads {threads} --index {input.index} --reads {input.reads} --freq-output {output.freqs} --count-output {output.counts}"

rule freqk_call:
    input:
        freqs="results/freqk/{variant}_vs_{backbone}/{ID}_freqs.txt",
        index="results/freqk/{variant}_vs_{backbone}/ref_index.txt"
    output:
        "results/freqk/{variant}_vs_{backbone}/{ID}_calls.txt"
    benchmark:
        "benchmarks/freqk_call/freqk_call_{variant}_vs_{backbone}_{ID}.bench"
    params:
        binary=config["freqk_binary"]
    shell:
        "{params.binary} call --index {input.index} -c {input.freqs} --output {output}"
