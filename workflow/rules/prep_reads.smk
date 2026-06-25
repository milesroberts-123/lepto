rule samtools_fastq:
    input:
        cram="{cram}"
    output:
        r1=temp("{config[output_dirs]["fastq"]}/{{ID}}_R1.fq"),
        r2=temp("{config[output_dirs]["fastq"]}/{{ID}}_R2.fq")
    conda: "../envs/bcftools.yaml"
    params:
        threads=4
    shell:
        """
        samtools fastq -1 {output.r1} -2 {output.r2} -n {input.cram} -@ {params.threads}
        """
