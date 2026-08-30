# Graph-genome construction and diversity analysis: build a pangenome
# graph from two haplotypes, map pool-seq reads with vg giraffe, and
# compute diversity/Fst with grenedalf over degenerate-site masks.


def count_vcf_contigs(vcfgz):
    """Return the number of contigs in a tabix-indexed VCF."""
    import subprocess
    result = subprocess.run(["tabix", "-l", vcfgz], capture_output=True, text=True)
    return len(result.stdout.strip().split("\n"))


rule add_prefix:
    input:
        lambda wildcards: config["bwa_refs"][wildcards.reference]
    output:
        "results/vg/{reference}_prefixed.fasta"
    conda: "../envs/featurecounts.yaml"
    shell:
        """
        seqkit replace -p '^' -r '{wildcards.reference}_' {input} > {output}
        """

rule create_sequence_dict:
    input:
        "results/vg/{backbone}_prefixed.fasta"
    output:
        "results/vg/{backbone}_prefixed.dict"
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools dict {input} > {output}"


rule minimap2_asm_paf:
    input:
        backbone="results/vg/{backbone}_prefixed.fasta",
        variant="results/vg/{variant}_prefixed.fasta"
    output:
        "results/vg/{variant}_vs_{backbone}.paf"
    conda: "../envs/minimap2.yaml"
    shell:
        "minimap2 -t {threads} -c --eqx -x asm10 -f 0.02 --cs {input.backbone} {input.variant} > {output}"


rule paftools_call_vcf:
    input:
        paf="results/vg/{variant}_vs_{backbone}.paf",
        backbone="results/vg/{backbone}_prefixed.fasta"
    output:
        "results/vg/{variant}_vs_{backbone}.vcf"
    conda: "../envs/minimap2.yaml"
    params:
        L=config["paftools_call_L"],
        q=config["paftools_call_q"]
    shell:
        "sort -k6,6 -k8,8n {input.paf} | paftools.js call -f {input.backbone} -L {params.L} -q {params.q} -s {wildcards.variant} - > {output}"


rule bgzip_tabix_vcf:
    input:
        "results/vg/{variant}_vs_{backbone}.vcf"
    output:
        vcfgz="results/vg/{variant}_vs_{backbone}.vcf.gz",
        tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi"
    conda: "../envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output.vcfgz} && tabix -p vcf {output.vcfgz}"


rule vg_autoindex:
    input:
        backbone="results/vg/{backbone}_prefixed.fasta",
        vcfgz="results/vg/{variant}_vs_{backbone}.vcf.gz",
        tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi"
    output:
        gbz="results/vg/{variant}_{backbone}_giraffe.giraffe.gbz",
        min="results/vg/{variant}_{backbone}_giraffe.shortread.withzip.min",
        dist="results/vg/{variant}_{backbone}_giraffe.dist",
        zip="results/vg/{variant}_{backbone}_giraffe.shortread.zipcodes"
    threads: lambda wildcards, input: count_vcf_contigs(input.vcfgz)
    conda: "../envs/vg.yaml"
    params:
        prefix=lambda wildcards: f"results/vg/{wildcards.variant}_{wildcards.backbone}_giraffe"
    shell:
        """
        # Split VCF by contig and use local scratch for temp files to speed up indexing.
        # Passing one VCF per contig lets vg autoindex parallelize across contigs,
        # and --tmp-dir on local storage avoids network filesystem overhead.
        LOCAL_TMP=/tmp/$USER/vgtmp
        mkdir -p "$LOCAL_TMP/vcf_split"
        trap "rm -rf $LOCAL_TMP" EXIT

        tabix -l {input.vcfgz} | while read contig; do
            tabix -h {input.vcfgz} "$contig" | bgzip > "$LOCAL_TMP/vcf_split/$contig.vcf.gz"
            tabix -p vcf "$LOCAL_TMP/vcf_split/$contig.vcf.gz"
        done

        vg autoindex --workflow sr-giraffe \
            -r {input.backbone} \
            $(for f in "$LOCAL_TMP/vcf_split"/*.vcf.gz; do echo -v "$f"; done) \
            -p {params.prefix} \
            -t {threads} \
            -M {resources.max_mem} \
            --tmp-dir "$LOCAL_TMP"

        # make indices read only so that vg giraffe isn't slow
        chmod 444 {params.prefix}.giraffe.gbz {params.prefix}.shortread.withzip.min {params.prefix}.dist {params.prefix}.shortread.zipcodes
        """


rule vg_giraffe:
    input:
        gbz="results/vg/{variant}_{backbone}_giraffe.giraffe.gbz",
        min="results/vg/{variant}_{backbone}_giraffe.shortread.withzip.min",
        dist="results/vg/{variant}_{backbone}_giraffe.dist",
        zip="results/vg/{variant}_{backbone}_giraffe.shortread.zipcodes",
        fastq="results/fastp/{ID}.fastq"
    output:
        temp("results/vg/{variant}_vs_{backbone}/{ID}/mapped.gam")
    conda: "../envs/vg.yaml"
    shell:
        "vg giraffe -Z {input.gbz} -m {input.min} -d {input.dist} -z {input.zip} -f {input.fastq} -b fast -t {threads} -p -o gam > {output}"


rule vg_surject:
    input:
        gam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.gam",
        gbz="results/vg/{variant}_{backbone}_giraffe.giraffe.gbz"
    output:
        temp("results/vg/{variant}_vs_{backbone}/{ID}/mapped.bam")
    conda: "../envs/vg.yaml"
    shell:
        "vg surject -x {input.gbz} -b -t {threads} {input.gam} > {output}"


rule samtools_sort_vg:
    input:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.bam"
    output:
        bam=temp("results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam"),
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools sort -m 2G -@ {threads} -o {output.bam} {input}"

rule samtools_index_vg:
    input:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam"
    output:
        bai=temp("results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam.bai")
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools index {input}"

rule picard_reorder_sam:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam.bai",
        dict="results/vg/{backbone}_prefixed.dict"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam"
    conda: "../envs/picard.yaml"
    shell:
        "picard ReorderSam INPUT={input.bam} OUTPUT={output} SEQUENCE_DICTIONARY={input.dict} ALLOW_CONTIG_LENGTH_DISCORDANCE=false ALLOW_INCOMPLETE_DICT_CONCORDANCE=false"


rule samtools_index_reordered:
    input:
        "results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai"
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools index {input}"

rule samtools_idxstats:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/idxstats.txt"
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools idxstats {input.bam} > {output}"

rule samtools_mpileup:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai",
        ref="results/vg/{backbone}_prefixed.fasta"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/{ID}.mpileup"
    conda: "../envs/bcftools.yaml"
    params:
        Q=config["samtools_mpileup_Q"],
        q=config["samtools_mpileup_q"]
    shell:
        "samtools mpileup -r hap1_scaffold_25 -f {input.ref} -Q {params.Q} -q {params.q} {input.bam} > {output}"

rule grenedalf_diversity_interval:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai",
        bed="results/degenotate/{backbone}/degeneracy-{site}-sites.bed",
        dict="results/vg/{backbone}_prefixed.dict",
        pool_sizes="pool_sizes.csv"
    output:
        div="results/vg/{variant}_vs_{backbone}/{ID}/diversity/{site}/grenedalf_results_interval_diversity.csv",
    conda: "../envs/grenedalf.yaml"
    params:
        window_width=config["grenedalf_window_width"],
        filter_sample_min_read_depth=config["grenedalf_filter_sample_min_read_depth"],
        window_average_policy=config["grenedalf_window_average_policy"],
        filter_sample_min_count=config["grenedalf_filter_sample_min_count"],
        min_map_qual=config["grenedalf_min_map_qual"],
        min_base_qual=config["grenedalf_min_base_qual"],
        subsample_max_read_depth=config["grenedalf_subsample_max_read_depth"]
    shell:
        """
        grenedalf diversity --window-type interval \\
            --window-interval-width {params.window_width} \\
            --pool-sizes {input.pool_sizes} \\
            --filter-sample-min-read-depth {params.filter_sample_min_read_depth} \\
            --filter-sample-min-count {params.filter_sample_min_count} \\
            --subsample-max-read-depth {params.subsample_max_read_depth} \\
            --window-average-policy {params.window_average_policy} \\
            --filter-region-bed {input.bed} \\
            --reference-genome-dict {input.dict} \\
            --out-dir $(dirname {output.div}) \\
            --file-prefix grenedalf_results_interval_ \\
            --sam-min-map-qual {params.min_map_qual} \\
            --sam-min-base-qual {params.min_base_qual} \\
            --allow-file-overwriting \\
            --threads {threads} \\
            --sam-path {input.bam}
        """

rule grenedalf_diversity_genome:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai",
        bed="results/degenotate/{backbone}/degeneracy-{site}-sites.bed",
        dict="results/vg/{backbone}_prefixed.dict",
        pool_sizes="pool_sizes.csv"
    output:
        div="results/vg/{variant}_vs_{backbone}/{ID}/diversity/{site}/grenedalf_results_genome_diversity.csv",
    conda: "../envs/grenedalf.yaml"
    params:
        window_width=config["grenedalf_window_width"],
        filter_sample_min_read_depth=config["grenedalf_filter_sample_min_read_depth"],
        window_average_policy=config["grenedalf_window_average_policy"],
        filter_sample_min_count=config["grenedalf_filter_sample_min_count"],
        min_map_qual=config["grenedalf_min_map_qual"],
        min_base_qual=config["grenedalf_min_base_qual"],
        subsample_max_read_depth=config["grenedalf_subsample_max_read_depth"]
    shell:
        """
        grenedalf diversity --window-type genome \\
            --pool-sizes {input.pool_sizes} \\
            --filter-sample-min-read-depth {params.filter_sample_min_read_depth} \\
            --filter-total-snp-min-count {params.filter_sample_min_count} \\
            --window-average-policy {params.window_average_policy} \\
            --filter-region-bed {input.bed} \\
            --reference-genome-dict {input.dict} \\
            --out-dir $(dirname {output.div}) \\
            --file-prefix grenedalf_results_genome_ \\
            --sam-min-map-qual {params.min_map_qual} \\
            --sam-min-base-qual {params.min_base_qual} \\
            --allow-file-overwriting \\
            --threads {threads} \\
            --sam-path {input.bam}
        """


rule grenedalf_frequency:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai",
        #bed="results/degenotate/{backbone}/degeneracy-{site}-sites.bed",
        dict="results/vg/{backbone}_prefixed.dict"
    output:
        fst="results/vg/{variant}_vs_{backbone}/{ID}/grenedalf_results_frequency.csv",
    conda: "../envs/grenedalf.yaml"
    params:
        window_width=config["grenedalf_window_width"],
        window_average_policy=config["grenedalf_window_average_policy"],
        filter_sample_min_count=config["grenedalf_filter_sample_min_count"],
        min_map_qual=config["grenedalf_min_map_qual"],
        min_base_qual=config["grenedalf_min_base_qual"]
    shell:
        """
        grenedalf frequency \\
            --write-sample-counts \\
            --write-sample-read-depth \\
            --write-sample-alt-freq \\
            --write-invariants \\
            --reference-genome-dict {input.dict} \\
            --out-dir $(dirname {output.fst}) \\
            --file-prefix grenedalf_results_ \\
            --sam-min-map-qual {params.min_map_qual} \\
            --sam-min-base-qual {params.min_base_qual} \\
            --allow-file-overwriting \\
            --threads {threads} \\
            --sam-path {input.bam}
        """


rule grenedalf_frequency_filter:
    input: "results/vg/{variant}_vs_{backbone}/{ID}/grenedalf_results_frequency.csv"
    output: "results/vg/{variant}_vs_{backbone}/{ID}/grenedalf_results_frequency_filtered.csv"
    shell: 
        """
        awk -F, '$1 == "hap1_scaffold_25" && $2 > 10960000 && $2 < 12050000 && $8 > 0.02 && $7 >= 10 && $6 > 1' {input} > {output}
        """ 


rule grenedalf_fst_interval:
    input:
        bam=expand("results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam", variant=[config["vg_ref_variant"]], backbone=[config["vg_ref_backbone"]], ID=sample_ids),
        bai=expand("results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai", variant=[config["vg_ref_variant"]], backbone=[config["vg_ref_backbone"]], ID=sample_ids),
        bed="results/degenotate/{backbone}/degeneracy-{site}-sites.bed",
        dict="results/vg/{backbone}_prefixed.dict",
        pool_sizes="pool_sizes.csv"
    output:
        fst="results/vg/{variant}_vs_{backbone}/fst/{site}/grenedalf_results_interval_fst.csv",
    conda: "../envs/grenedalf.yaml"
    params:
        window_width=config["grenedalf_window_width"],
        filter_sample_min_read_depth=config["grenedalf_filter_sample_min_read_depth"],
        window_average_policy=config["grenedalf_window_average_policy"],
        filter_total_snp_min_count=config["grenedalf_filter_sample_min_count"],
        min_map_qual=config["grenedalf_min_map_qual"],
        min_base_qual=config["grenedalf_min_base_qual"],
        method=config["grenedalf_fst_method"]
    shell:
        """
        grenedalf fst --window-type interval \\
            --window-interval-width {params.window_width} \\
            --pool-sizes {input.pool_sizes} \\
            --filter-region-bed {input.bed} \\
            --filter-sample-min-read-depth {params.filter_sample_min_read_depth} \\
            --filter-total-snp-min-count {params.filter_total_snp_min_count} \\
            --window-average-policy {params.window_average_policy} \\
            --reference-genome-dict {input.dict} \\
            --out-dir results/vg/{wildcards.variant}_vs_{wildcards.backbone}/fst/{wildcards.site}/ \\
            --file-prefix grenedalf_results_interval_ \\
            --sam-min-map-qual {params.min_map_qual} \\
            --sam-min-base-qual {params.min_base_qual} \\
            --method {params.method} \\
            --allow-file-overwriting \\
            --write-pi-tables \\
            --threads {threads} \\
            --sam-path {input.bam}
        """

rule grenedalf_fst_genome:
    input:
        bam=expand("results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam", variant=[config["vg_ref_variant"]], backbone=[config["vg_ref_backbone"]], ID=sample_ids),
        bai=expand("results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai", variant=[config["vg_ref_variant"]], backbone=[config["vg_ref_backbone"]], ID=sample_ids),
        dict="results/vg/{backbone}_prefixed.dict",
        bed="results/degenotate/{backbone}/degeneracy-{site}-sites.bed",
        pool_sizes="pool_sizes.csv"
    output:
        fst="results/vg/{variant}_vs_{backbone}/fst/{site}/grenedalf_results_genome_fst.csv",
    conda: "../envs/grenedalf.yaml"
    params:
        filter_sample_min_read_depth=config["grenedalf_filter_sample_min_read_depth"],
        window_average_policy=config["grenedalf_window_average_policy"],
        filter_total_snp_min_count=config["grenedalf_filter_sample_min_count"],
        min_map_qual=config["grenedalf_min_map_qual"],
        min_base_qual=config["grenedalf_min_base_qual"],
        method=config["grenedalf_fst_method"]
    shell:
        """
        grenedalf fst --window-type genome \\
            --pool-sizes {input.pool_sizes} \\
            --filter-sample-min-read-depth {params.filter_sample_min_read_depth} \\
            --filter-total-snp-min-count {params.filter_total_snp_min_count} \\
            --window-average-policy {params.window_average_policy} \\
            --filter-region-bed {input.bed} \\
            --reference-genome-dict {input.dict} \\
            --out-dir results/vg/{wildcards.variant}_vs_{wildcards.backbone}/fst/{wildcards.site}/ \\
            --file-prefix grenedalf_results_genome_ \\
            --sam-min-map-qual {params.min_map_qual} \\
            --sam-min-base-qual {params.min_base_qual} \\
            --method {params.method} \\
            --allow-file-overwriting \\
            --write-pi-tables \\
            --threads {threads} \\
            --sam-path {input.bam}
        """

