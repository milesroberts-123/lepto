rule samtools_fastq:
    input:
    output:
    conda: "../envs/bcftools.yaml"
    shell:
        "samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}"

rule fastp:
    input:

    output:
        duread1=temp("fastp_results/dedup_unpaired_R1_{ID}.fastq.gz"),
        duread2=temp("fastp_results/dedup_unpaired_R2_{ID}.fastq.gz"),
        pread1=temp("fastp_results/trimmed_paired_R1_{ID}.fastq.gz"),
        pread2=temp("fastp_results/trimmed_paired_R2_{ID}.fastq.gz"),
        uread1=temp("fastp_results/trimmed_unpaired_R1_{ID}.fastq.gz"),
        uread2=temp("fastp_results/trimmed_unpaired_R2_{ID}.fastq.gz"),
        jsonR1R2=temp("fastp_results/{ID}_R1R2.json"),
        jsonU1=temp("fastp_results/{ID}_U1.json"),
        jsonU2=temp("fastp_results/{ID}_U2.json"),
    conda:
        "../envs/fastp.yaml"
    log:
        "logs/fastp/{ID}.log",
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
        nBaseLimit=config["nBaseLimit"]
    shell:
        """
        # remove duplicates, do read correction, drop low quality reads
        fastp --thread {threads} --n_base_limit {params.nBaseLimit} -u {params.unqualLimit} -q {params.qualThresh} -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} --json {output.jsonR1R2} --dedup --correction -i {input.read1} -I {input.read2} -o {output.dpread1} -O {output.dpread2} --unpaired1 {output.duread1} --unpaired2 {output.duread2}
        """

