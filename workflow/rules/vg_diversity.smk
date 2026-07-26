vg_backbone = config["vg_ref_backbone"]
vg_variant = config["vg_ref_variant"]


rule minimap2_asm_paf:
    input:
        backbone=lambda wildcards: config["bwa_refs"][wildcards.backbone],
        variant=lambda wildcards: config["bwa_refs"][wildcards.variant]
    output:
        "results/vg/{variant}_vs_{backbone}.paf"
    conda: "../envs/minimap2.yaml"
    shell:
        "minimap2 -t {threads} -cx asm10 -f 0.02 --cs {input.backbone} {input.variant} > {output}"


rule paftools_call_vcf:
    input:
        paf="results/vg/{variant}_vs_{backbone}.paf",
        backbone=lambda wildcards: config["bwa_refs"][wildcards.backbone]
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
        backbone=lambda wildcards: config["bwa_refs"][wildcards.backbone],
        vcfgz="results/vg/{variant}_vs_{backbone}.vcf.gz",
        tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi"
    output:
        gbz="results/vg/{variant}_{backbone}_giraffe.giraffe.gbz",
        min="results/vg/{variant}_{backbone}_giraffe.min",
        dist="results/vg/{variant}_{backbone}_giraffe.dist"
    conda: "../envs/vg.yaml"
    params:
        prefix=lambda wildcards: f"results/vg/{wildcards.variant}_{wildcards.backbone}_giraffe"
    shell:
        "vg autoindex --workflow sr-giraffe -r {input.backbone} -v {input.vcfgz} -p {params.prefix} -t {threads}"


rule vg_giraffe:
    input:
        gbz="results/vg/{variant}_{backbone}_giraffe.giraffe.gbz",
        min="results/vg/{variant}_{backbone}_giraffe.min",
        dist="results/vg/{variant}_{backbone}_giraffe.dist",
        fastq="results/fastp/{ID}.fastq"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.gam"
    conda: "../envs/vg.yaml"
    shell:
        "vg giraffe -Z {input.gbz} -m {input.min} -d {input.dist} -f {input.fastq} -p > {output}"


rule vg_surject:
    input:
        gam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.gam",
        gbz="results/vg/{variant}_{backbone}_giraffe.giraffe.gbz"
    output:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.bam"
    conda: "../envs/vg.yaml"
    params:
        prefixes=lambda wildcards: " ".join(f"-p {p}" for p in config["vg_surject_prefixes"])
    shell:
        "vg surject -x {input.gbz} {params.prefixes} -b -t {threads} {input.gam} > {output}"


rule samtools_sort_index_vg:
    input:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.bam"
    output:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam.bai"
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output.bam} {input} && samtools index {output.bam}"


rule grenedalf_diversity:
    input:
        "results/vg/{variant}_vs_{backbone}/{ID}/mapped.sorted.bam"
    output:
        directory("results/vg/{variant}_vs_{backbone}/{ID}/diversity")
    conda: "../envs/grenedalf.yaml"
    params:
        window_width=config["grenedalf_window_width"],
        pool_sizes=config["grenedalf_pool_sizes"],
        filter_min_count=config["grenedalf_filter_min_count"],
        filter_min_read_depth=config["grenedalf_filter_min_read_depth"],
        window_average_policy=config["grenedalf_window_average_policy"]
    shell:
        "grenedalf diversity --window-type sliding --window-width {params.window_width} --pool-sizes {params.pool_sizes} --filter-min-count {params.filter_min_count} --filter-min-read-depth {params.filter_min_read_depth} --window-average-policy {params.window_average_policy} --file-prefix {output}/ {input}"


rule vg_diversity_all:
    input:
        expand("results/vg/{variant}_vs_{backbone}/{ID}/diversity",
               variant=[vg_variant], backbone=[vg_backbone], ID=sample_ids)
