# SV calling with svim-asm: chunk both assemblies into fixed-size windows
# (to keep contigs within SAM length limits), align variant assembly to
# backbone with minimap2, sort/index, then call structural variants in
# haploid mode.


rule seqkit_sliding:
    input:
        "results/vg/{reference}_prefixed.fasta"
    output:
        "results/svim/{reference}_chunked.fasta"
    conda: "../envs/minimap2.yaml"
    params:
        chunk=config["svim_chunk_size"]
    shell:
        """
        mkdir -p results/svim
        seqkit sliding -W {params.chunk} -s {params.chunk} -g -o {output} {input}
        """


rule minimap2_svim_asm:
    input:
        ref="results/svim/{backbone}_chunked.fasta",
        qry="results/svim/{variant}_chunked.fasta"
    output:
        "results/svim/{variant}_vs_{backbone}.sam"
    conda: "../envs/minimap2.yaml"
    params:
        x=config["svim_minimap_x"],
    shell:
        "minimap2 -a -f 0.02 -r2k -x {params.x} --cs -t {threads} {input.ref} {input.qry} > {output}"


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
        ref="results/svim/{backbone}_chunked.fasta"
    output:
        "results/svim/{variant}_vs_{backbone}/variants.vcf"
    conda: "../envs/svim.yaml"
    params:
        gap_tolerance=config["svim_reference_gap_tolerance"],
        max_sv_size=config["svim_max_sv_size"]
    shell:
        "svim-asm haploid --types DEL,INS,INV --reference_gap_tolerance {params.gap_tolerance} --max_sv_size {params.max_sv_size} results/svim/{wildcards.variant}_vs_{wildcards.backbone} {input.bam} {input.ref}"
