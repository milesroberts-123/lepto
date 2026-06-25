rule samtools_fastq:
    input:
        cram="../config/cram/{ID}.cram"
    output:
        r1=temp("results/samtools/{{ID}}_R1.fq"),
        r2=temp("results/samtools/{{ID}}_R2.fq")
    conda: "../envs/bcftools.yaml"
    shell:
        """
        samtools fastq -1 {output.r1} -2 {output.r2} -n {input.cram} -@ {params.threads}
        """

# ---- read trimming ----
rule fastp:
    input:
        r1=temp("results/samtools/{{ID}}_R1.fq"),
        r2=temp("results/samtools/{{ID}}_R2.fq")
    output:
        pread1=temp("results/fastp/paired_R1_{ID}.fastq.gz"),
        pread2=temp("results/fastp/paired_R2_{ID}.fastq.gz"),
        uread1=temp("results/fastp/unpaired_R1_{ID}.fastq.gz"),
        uread2=temp("results/fastp/unpaired_R2_{ID}.fastq.gz"),
        jsonR1R2=temp("fastp_results/{ID}_R1R2.json"),
        jsonU1=temp("fastp_results/{ID}_U1.json"),
        jsonU2=temp("fastp_results/{ID}_U2.json"),
    conda:
        "../envs/fastp.yaml"
    params:
        unqualLimit=config["unqualLimit"],
        k=config["k"],
        qualThresh=config["qualThresh"],
        windowLength=config["windowLength"],
        nBaseLimit=config["nBaseLimit"]
    shell:
        """
        # remove duplicates, do read correction, drop low quality reads
        fastp --thread {threads} --n_base_limit {params.nBaseLimit} -u {params.unqualLimit} -q {params.qualThresh} --dedup --correction -l {params.k} --cut_tail --cut_tail_window_size {params.windowLength} --cut_tail_mean_quality {params.qualThresh} -i {input.r1} -I {input.r2} -o {output.dpread1} -O {output.dpread2} --unpaired1 {output.duread1} --unpaired2 {output.duread2} --json {output.jsonR1R2}
        """

# ---- per-sample k-mer counting ----
rule kmc_count:
    input:
        dpread1=temp("results/fastp/dedup_paired_R1_{{ID}}.fastq.gz"),
        dpread2=temp("results/fastp/dedup_paired_R2_{{ID}}.fastq.gz"),
        uread1=temp("results/fastp/unpaired_R1_{{ID}}.fastq.gz"),
        uread2=temp("results/fastp/unpaired_R2_{{ID}}.fastq.gz")
    output:
        db=temp("results/kmc/{{ID}}/kmc_db")
    conda: "../envs/kmc.yaml"
    params:
        mincount=config["mincount"],
        maxcount=config.get("maxcount", "auto"),
        k=config["k"]
    shell:
        """
        local_tmp=tmp_kmc_{wildcards.ID}
        rm -rf "$local_tmp"
        mkdir -p "$local_tmp"

        ls -d {input} > {log}.kmc_input

        if [ "{params.maxcount}" = "auto" ]; then
            kmc -m15 -t{threads} -ci{params.mincount} -k{params.k} \\
                @{output.db}.inlist \\
                {output.db} "$local_tmp"
        else
            kmc -m15 -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} \\
                @{output.db}.inlist \\
                {output.db} "$local_tmp"
        fi

        # remove tmp directory
        rm -rf "$local_tmp"
        """

# ---- Union k-mers within a group ----
rule kmc_intersect_group:
    input:
        dbs=expand("results/kmc/{{ID}}/kmc_db", ID=lambda wildcards: samples_by_group[wildcards.group])
    output:
        pre="results/grouped/{{group}}.kmc_pre",
        suf="results/grouped/{{group}}.kmc_suf"
    conda: "../envs/kmc.yaml"
    params:
        depth=config.get("kmc_intersect_depth", 100)
    log:
        "logs/kmc_group/{group}.log"
    shell:
        """
        ls {input.dbs} > {log}.inlist

        kmc_tools complex @{log}.inlist union -kci1 -kcs{params.depth} {output.pre}
        """

# ---- Subtract non-target group k-mers ----
rule kmc_subtract:
    input:
        target_db="results/grouped/{{group}}.kmc_pre",
        other_dbs=expand("results/grouped/{{group}}.kmc_pre",
                        group=lambda wildcards: [g for g in groups if g != wildcards.group])
    output:
        pre="results/grouped/{{group}}_specific.kmc_pre",
        suf="results/grouped/{{group}}_specific.kmc_suf"
    conda: "../envs/kmc.yaml"
    shell:
        """
        echo {input.target_db} kci1 > {log}.complex_input
        printf '%s kci1\n' {input.other_dbs} >> {log}.complex_input

        kmc_tools complex @{log}.complex_input subtract {output.pre}
        """

# ---- Filter reads containing group-specific k-mers (per sample) ----
rule kmc_filter_reads:
    input:
        kmc_db="results/grouped/{{group}}_specific.kmc_pre",
        fastq="results/fastp/{{read_type}}/{{read_name}}.fastq.gz"
    output:
        filtered="results/filtered/{{group}}/{{read_type}}_{{ID}}_filtered.fastq"
    conda: "../envs/kmc.yaml"
    params:
        min_support=config["filter_min_kmer_support"]
    shell:
        """
        kmc_tools extract -db {input.kmc_db} -ci {params.min_support} \\
            {input.fastq} {output.filtered}
        """

# ---- Combine filtered reads for a group and assemble ----
rule combine_group_reads:
    input:
        filtered=temp("results/filtered/{{group}}/dp_{{ID}}_filtered.fastq",
                     ID=lambda wildcards: samples_by_group[wildcards.group]),
        unfiltered=temp("results/filtered/{{group}}/u_{{ID}}_filtered.fastq",
                       ID=lambda wildcards: samples_by_group[wildcards.group])
    output:
        all_reads="{results/combined/{{group}}_combined.fastq.gz"
    shell:
        """
        cat {input.filtered} {input.unfiltered} | gzip > {output.all_reads}
        """

# ---- MetaSPAdes assembly of group-specific reads ----
rule metaspades:
    input:
        reads="results/combined/{{group}}_combined.fastq.gz"
    output:
        assembly="results/assembly/{{group}}_assembly.fasta"
    conda: "../envs/metaspades.yaml"
    params:
        kmers=config["metaspades"]["kmers"],
    shell:
        """
        mkdir -p {output.assembly}

        metaspades.py \\
            --meta \\
            -m {params.mem} \\
            -t {params.threads} \\
            -k {params.kmers} \\
            -1 {input.reads} \\
            -o {output.assembly}
        """
