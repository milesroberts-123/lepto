vg_backbone = config["vg_ref_backbone"]
vg_variant = config["vg_ref_variant"]


def count_vcf_contigs(vcfgz):
    import subprocess
    result = subprocess.run(["tabix", "-l", vcfgz], capture_output=True, text=True)
    return len(result.stdout.strip().split("\n"))


rule add_prefix:
    input:
        lambda wildcards: config["bwa_refs"][wildcards.reference]
    output:
        "results/vg/{reference}_prefixed.fasta"
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
        "minimap2 -t {threads} -cx asm10 -f 0.02 --cs {input.backbone} {input.variant} > {output}"


rule paftools_call_vcf:
    input:
        paf="results/vg/{variant}_vs_{backbone}.paf",
        backbone="results/vg/{backbone}_prefixed.fasta"
    output:
        "results/vg/{variant}_vs_{backbone}.vcf"
    conda: "../envs/minimap2.yaml"
    shell:
        "sort -k6,6 -k8,8n {input.paf} | paftools.js call -f {input.backbone} -s {wildcards.variant} - > {output}"


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
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.gam"
    conda: "../envs/vg.yaml"
    shell:
        "vg giraffe -Z {input.gbz} -m {input.min} -d {input.dist} -z {input.zip} -f {input.fastq} -b fast -t {threads} -p -o gam > {output}"


rule vg_surject:
    input:
        gam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.gam",
        gbz="results/vg/{variant}_{backbone}_giraffe.giraffe.gbz"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.bam"
    conda: "../envs/vg.yaml"
    shell:
        "vg surject -x {input.gbz} -b -t {threads} {input.gam} > {output}"


rule samtools_sort_index_vg:
    input:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.bam"
    output:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam.bai"
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools sort -m 2G -@ {threads} -o {output.bam} {input} && samtools index {output.bam}"


rule picard_reorder_sam:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam",
        dict="results/vg/{backbone}_prefixed.dict"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam"
    conda: "../envs/picard.yaml"
    shell:
        "picard ReorderSam INPUT={input.bam} OUTPUT={output} SEQUENCE_DICTIONARY={input.dict} ALLOW_CONTIG_LENGTH_DISCORDANCE=false ALLOW_INCOMPLETE_DICT_CONCORDANCE=false"


rule samtools_index_reordered:
    input:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam.bai"
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools index {input}"


rule grenedalf_diversity:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam.bai",
        bed="results/degenotate/{backbone}/degeneracy-{site}-sites.bed",
        dict="results/vg/{backbone}_prefixed.dict"
    output:
        directory("results/vg/{variant}_vs_{backbone}/{ID}/diversity/{site}")
    conda: "../envs/grenedalf.yaml"
    params:
        window_width=config["grenedalf_window_width"],
        pool_sizes=config["grenedalf_pool_sizes"],
        filter_min_count=config["grenedalf_filter_min_count"],
        filter_min_read_depth=config["grenedalf_filter_min_read_depth"],
        window_average_policy=config["grenedalf_window_average_policy"],
        filter_sample_min_count=config["grenedalf_filter_sample_min_count"]
    shell:
        """
        grenedalf diversity --window-type interval \\
            --window-interval-width {params.window_width} \\
            --pool-sizes {params.pool_sizes} \\
            --filter-sample-min-read-depth {params.filter_min_read_depth} \\
            --filter-sample-min-count {params.filter_sample_min_count} \\
            --window-average-policy {params.window_average_policy} \\
            --filter-mask-total-bed {input.bed} \\
            --filter-mask-total-bed-invert \\
            --reference-genome-dict {input.dict} \\
            --out-dir {output}/ \\
            --file-prefix pi_windows_ \\
            --sam-path {input.bam}
        """


rule grenedalf_frequency:
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam.bai",
        bed="results/degenotate/{backbone}/degeneracy-{site}-sites.bed",
        dict="results/vg/{backbone}_prefixed.dict"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/frequency/{site}/grenedalf_results_frequency.txt"
    conda: "../envs/grenedalf.yaml"
    params:
        window_width=config["grenedalf_window_width"],
        pool_sizes=config["grenedalf_pool_sizes"],
        filter_min_count=config["grenedalf_filter_min_count"],
        filter_min_read_depth=config["grenedalf_filter_min_read_depth"],
        window_average_policy=config["grenedalf_window_average_policy"],
        filter_sample_min_count=config["grenedalf_filter_sample_min_count"]
    shell:
        """
        grenedalf frequency --window-type interval \\
            --window-interval-width {params.window_width} \\
            --pool-sizes {params.pool_sizes} \\
            --filter-sample-min-read-depth {params.filter_min_read_depth} \\
            --filter-sample-min-count {params.filter_sample_min_count} \\
            --window-average-policy {params.window_average_policy} \\
            --filter-mask-total-bed {input.bed} \\
            --filter-mask-total-bed-invert \\
            --reference-genome-dict {input.dict} \\
            --out-dir results/vg/{wildcards.variant}_vs_{wildcards.backbone}/{wildcards.ID}/frequency/{wildcards.site} \\
            --file-prefix grenedalf_results_ \\
            --sam-path {input.bam}
        """


rule grenedalf_fst:
    input:
        bam=expand("results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam", variant=[vg_variant], backbone=[vg_backbone], ID=sample_ids),
        bai=expand("results/vg/{variant}_vs_{backbone}/{ID}/mapped.reordered.bam.bai", variant=[vg_variant], backbone=[vg_backbone], ID=sample_ids),
        bed="results/degenotate/{backbone}/degeneracy-{site}-sites.bed",
        dict="results/vg/{backbone}_prefixed.dict"
    output:
        "results/vg/{variant}_vs_{backbone}/fst/{site}/grenedalf_results_fst.txt"
    conda: "../envs/grenedalf.yaml"
    params:
        window_width=config["grenedalf_window_width"],
        pool_sizes=config["grenedalf_pool_sizes"],
        filter_min_count=config["grenedalf_filter_min_count"],
        filter_min_read_depth=config["grenedalf_filter_min_read_depth"],
        window_average_policy=config["grenedalf_window_average_policy"],
        filter_sample_min_count=config["grenedalf_filter_sample_min_count"],
        min_map_qual=config["grenedalf_min_map_qual"],
        min_base_qual=config["grenedalf_min_base_qual"]
    shell:
        """
        grenedalf fst --window-type interval \\
            --window-interval-width {params.window_width} \\
            --pool-sizes {params.pool_sizes} \\
            --filter-sample-min-read-depth {params.filter_min_read_depth} \\
            --filter-sample-min-count {params.filter_sample_min_count} \\
            --window-average-policy {params.window_average_policy} \\
            --filter-mask-total-bed {input.bed} \\
            --filter-mask-total-bed-invert \\
            --reference-genome-dict {input.dict} \\
            --out-dir results/vg/{wildcards.variant}_vs_{wildcards.backbone}/fst/{wildcards.site}/ \\
            --file-prefix grenedalf_results_ \\
            --sam-min-map-qual {params.min_map_qual} \\
            --sam-min-base-qual {params.min_base_qual} \\
            --sam-path {input.bam}
        """

rule vg_diversity_all:
    input:
        expand("results/vg/{variant}_vs_{backbone}/{ID}/frequency/{site}/grenedalf_results_frequency.txt",
            variant=[vg_variant],
            backbone=[vg_backbone],
            ID=sample_ids,
            site = ["cds", "four", "zero"]),
        expand("results/vg/{variant}_vs_{backbone}/fst/{site}/grenedalf_results_fst.txt",
            variant=[vg_variant],
            backbone=[vg_backbone],
            site = ["cds", "four", "zero"]),
        expand("results/vg/{variant}_vs_{backbone}/{ID}/diversity/{site}",
            variant=[vg_variant],
            backbone=[vg_backbone],
            ID=sample_ids,
            site = ["cds", "four", "zero"]),

