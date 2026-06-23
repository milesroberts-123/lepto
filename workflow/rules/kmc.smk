rule kmc:
    input:
        pread1=temp("fastp_results/trimmed_paired_R1_{ID}.fastq.gz"),
        pread2=temp("fastp_results/trimmed_paired_R2_{ID}.fastq.gz"),
        uread1=temp("fastp_results/trimmed_unpaired_R1_{ID}.fastq.gz"),
        uread2=temp("fastp_results/trimmed_unpaired_R2_{ID}.fastq.gz"),
    output:
        tmp_pre=temp("results/kmc/{ID}.kmc_pre"),
        tmp_suf=temp("results/kmc/{ID}.kmc_suf")
        tmp_list=temp("results/kmc/input_files_{ID}.txt")        
    conda:
        "../envs/kmc.yaml"
    params:
        mincount=config["mincount"],
        maxcount=config["maxcount"],
        k=config["k"],
    shell:
        """
        # create directory
        if [ -d "tmp_kmc_{wildcards.ID}" ]; then
            rm -r tmp_kmc_{wildcards.ID}
        fi

        mkdir tmp_kmc_{wildcards.ID}

        # list of input files
        ls -d {input} > {output.tmp_list}

        # count k-mers
        kmc -m15 -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} @{output.tmp_list} kmc_db_{wildcards.ID} tmp_kmc_{wildcards.ID}

        # delete tmp directories
        rm -r tmp_kmc_{wildcards.ID}
        """


rule kmc_intersect_group:
    input:
        pre=expand("results/kmc/{ID}.kmc_pre", ID=lookup())
    output:
        pre="results/grouped/{group}.kmc_pre"
        suf="results/grouped/{group}.kmc_suf"
    shell:
        """
        kmc_tools complex {input}
        """

rule kmc_subtract:

