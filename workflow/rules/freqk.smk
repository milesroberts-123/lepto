# Estimate allele frequencies of known variants in pool-seq reads from
# k-mer counts. Variants come from the vg pipeline's paftools-call VCF on
# the hap1 backbone (run vg_diversity_all first: this branch expects the
# VCF, its tabix index, and the prefixed backbone fasta to already
# exist); reads are the same fastp-cleaned pooled fastqs fed to vg
# giraffe. Before indexing, the VCF is normalized and subset to
# biallelic SNPs, excluding sites polymorphic (0.02 < alt freq < 0.98)
# in any W-group pool per the per-pool grenedalf frequency tables. The
# freqk binary is provisioned outside the repo and referenced by the
# freqk_binary config key (see AGENTS.md).

rule freqk_exclusion_sites:
    # Unsorted, un-deduplicated stream of sites polymorphic in any
    # W-group pool, from the per-pool grenedalf frequency tables (column
    # 6 alt count, 7 depth, 8 alt freq; 0-based BED). Consuming these
    # tables couples freqk to the vg mapping chain (see AGENTS.md
    # reverted-approaches entry): freqk-only targets will schedule the
    # vg chain if the CSVs are missing.
    input:
        expand("results/vg/{variant}_vs_{backbone}/{ID}/grenedalf_results_frequency.csv",
               variant=[config["vg_ref_variant"]],
               backbone=[config["vg_ref_backbone"]],
               ID=w_sample_ids)
    output:
        temp("results/freqk/{variant}_vs_{backbone}/polymorphic_sites.unsorted.bed")
    params:
        min_alt_count=config["freqk_min_alt_count"],
        min_depth=config["freqk_min_depth"],
        min_freq=config["freqk_min_polymorphic_freq"],
        max_freq=config["freqk_max_polymorphic_freq"]
    benchmark:
        "benchmarks/freqk_exclusion_bed/freqk_exclusion_sites_{variant}_vs_{backbone}.bench"
    shell:
        """
        mkdir -p results/freqk/{wildcards.variant}_vs_{wildcards.backbone}
        for table in {input}; do
            awk -F, -v c={params.min_alt_count} -v d={params.min_depth} \\
                -v f={params.min_freq} -v x={params.max_freq} \\
                '$1 != "CHROM" && $6 > c && $7 >= d && $8 > f && $8 < x \\
                 {{print $1"\\t"$2-1"\\t"$2}}' "$table"
        done > {output}
        """

rule freqk_exclusion_bed:
    # Sort and merge overlapping/adjacent polymorphic sites into the
    # final 0-based BED consumed by freqk_filter_vcf.
    input:
        "results/freqk/{variant}_vs_{backbone}/polymorphic_sites.unsorted.bed"
    output:
        "results/freqk/{variant}_vs_{backbone}/polymorphic_exclusion.bed"
    conda: "../envs/bedtools.yaml"
    benchmark:
        "benchmarks/freqk_exclusion_bed/freqk_exclusion_bed_{variant}_vs_{backbone}.bench"
    shell:
        """
        sort -k1,1 -k2,2n {input} | bedtools merge -i - > {output}
        """


rule freqk_filter_vcf:
    input:
        "results/vg/{variant}_vs_{backbone}.vcf.gz",
        tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi",
        excl="results/freqk/{variant}_vs_{backbone}/polymorphic_exclusion.bed"
    output:
        vcfgz=temp("results/freqk/{variant}_vs_{backbone}/norm.vcf.gz"),
        tbi=temp("results/freqk/{variant}_vs_{backbone}/norm.vcf.gz.tbi")
    conda: "../envs/bcftools.yaml"
    params:
        freqk_min_qual=config["freqk_min_qual"]
    shell:
        """
        mkdir -p results/freqk/{wildcards.variant}_vs_{wildcards.backbone}
        bcftools norm -m -any {input[0]} \\
            | bcftools view -T ^{input.excl} \\
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
