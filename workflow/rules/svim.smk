# SV calling with svim-asm: align variant assembly to backbone with
# minimap2, sort/index, then call structural variants in haploid mode.

svim_backbone = config["vg_ref_backbone"]
svim_variant = config["vg_ref_variant"]


rule minimap2_svim_asm:
    input:
        ref="results/vg/{backbone}_prefixed.fasta",
        qry="results/vg/{variant}_prefixed.fasta"
    output:
        "results/svim/{variant}_vs_{backbone}.sam"
    conda: "../envs/minimap2.yaml"
    params:
        x=config["svim_minimap_x"],
        r2k=config["svim_minimap_r2k"]
    shell:
        "minimap2 -a -x {params.x} --cs -r2k -t {threads} {input.ref} {input.qry} > {output}"


rule samtools_sort_svim:
    input:
        "results/svim/{variant}_vs_{backbone}.sam"
    output:
        temp("results/svim/{variant}_vs_{backbone}.sorted.bam")
    conda: "../envs/minimap2.yaml"
    params:
        m=config["svim_samtools_sort_m"]
    shell:
        "samtools sort -m{params.m} -@ {threads} -o {output} {input}"


rule samtools_index_svim:
    input:
        "results/svim/{variant}_vs_{backbone}.sorted.bam"
    output:
        "results/svim/{variant}_vs_{backbone}.sorted.bam.bai"
    conda: "../envs/minimap2.yaml"
    shell:
        "samtools index {input}"


rule svim_asm_run:
    input:
        bam="results/svim/{variant}_vs_{backbone}.sorted.bam",
        bai="results/svim/{variant}_vs_{backbone}.sorted.bam.bai",
        ref="results/vg/{backbone}_prefixed.fasta"
    output:
        "results/svim/{variant}_vs_{backbone}/svim-asm.vcf"
    conda: "../envs/svim.yaml"
    shell:
        "svim-asm haploid results/svim/{wildcards.variant}_vs_{wildcards.backbone} {input.bam} {input.ref}"


rule svim_all:
    input:
        expand("results/svim/{variant}_vs_{backbone}/svim-asm.vcf",
               variant=[svim_variant],
               backbone=[svim_backbone])
