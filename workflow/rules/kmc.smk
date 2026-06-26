rule samtools_fastq:
    input:
        cram=lambda wildcards: cram_paths[wildcards.ID]
    output:
        temp("results/samtools/{ID}.fq")
    conda: "../envs/bcftools.yaml"
    shell:
        """
        samtools fastq -0 {output} -n {input.cram} -@ {threads}
        """

rule fastp:
    input:
        "results/samtools/{ID}.fq",
    output:
        fastq=temp("results/fastp/{ID}.fastq"),
        json="results/fastp/{ID}.json"
    conda:
        "../envs/fastp.yaml"
    params:
        unqual_limit=config["fastp_unqual_limit"],
        min_len=config["fastp_min_len"],
        qual_thresh=config["fastp_qual_thresh"],
        window_length=config["fastp_window_length"],
        n_base_limit=config["fastp_n_base_limit"]
    shell:
        """
        fastp --thread {threads} \
            --n_base_limit {params.n_base_limit} \
            -u {params.unqual_limit} \
            -q {params.qual_thresh} \
            --dedup --correction \
            -l {params.min_len} \
            --cut_tail \
            --cut_tail_window_size {params.window_length} \
            --cut_tail_mean_quality {params.qual_thresh} \
            --json {output.json} \
            -i {input} \
            -o {output.fastq}
        """

rule kmc_count:
    input:
        "results/fastp/{ID}.fastq"
    output:
        pre=temp("results/kmc/{ID}/kmc_db.kmc_pre"),
        suf=temp("results/kmc/{ID}/kmc_db.kmc_suf")
    conda: "../envs/kmc.yaml"
    params:
        mincount=config["kmc_mincount"],
        maxcount=config["kmc_maxcount"],
        k=config["kmc_k"]
    shell:
        """
        local_tmp=tmp_kmc_{wildcards.ID}
        rm -rf "$local_tmp"
        mkdir -p "$local_tmp"
        mkdir -p results/kmc/{wildcards.ID}

        if [ "{params.maxcount}" = "auto" ]; then
            kmc -m28 -sm -t{threads} -ci{params.mincount} -k{params.k} \
                {input} \
                results/kmc/{wildcards.ID}/kmc_db "$local_tmp"
        else
            kmc -m28 -sm -t{threads} -ci{params.mincount} -cs{params.maxcount} -k{params.k} \
                {input} \
                results/kmc/{wildcards.ID}/kmc_db "$local_tmp"
        fi

        rm -rf "$local_tmp"
        """

rule kmc_intersect_group:
    input:
        dbs=expand(["results/kmc/{ID}/kmc_db.kmc_pre", "results/kmc/{ID}/kmc_db.kmc_suf"],
                   ID=lambda wildcards: samples_by_group[wildcards.group])
    output:
        pre=temp("results/grouped/{group}.kmc_pre"),
        suf=temp("results/grouped/{group}.kmc_suf"),
        complex=temp("results/grouped/{group}.complex")
    conda: "../envs/kmc.yaml"
    params:
        depth=config.get("kmc_intersect_depth", 100)
    shell:
        """
        mkdir -p results/grouped

        {{
            echo "INPUT:"
            printf '%s\\n' {input.dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{print "set" NR " = " $0 " -ci1"}}'
            echo "OUTPUT:"
            printf "results/grouped/{wildcards.group} = "
            printf '%s\\n' {input.dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{printf "%sset%d", (NR>1?" + ":""), NR}} END{{print ""}}'
            echo "OUTPUT_PARAMS:"
            echo "-cs{params.depth}"
        }} > {output.complex}

        kmc_tools -t{threads} complex {output.complex}
        """

rule kmc_subtract:
    input:
        target_db=["results/grouped/{group}.kmc_pre", "results/grouped/{group}.kmc_suf"],
        other_dbs=expand(["results/grouped/{other}.kmc_pre", "results/grouped/{other}.kmc_suf"],
                        other=lambda wildcards: [g for g in groups if g != wildcards.group])
    output:
        pre=temp("results/specific/{group}_specific.kmc_pre"),
        suf=temp("results/specific/{group}_specific.kmc_suf"),
        complex=temp("results/specific/{group}_specific.complex")
    conda: "../envs/kmc.yaml"
    shell:
        """
        mkdir -p results/specific

        if [ -z "{input.other_dbs}" ]; then
            cp results/grouped/{wildcards.group}.kmc_pre {output.pre}
            cp results/grouped/{wildcards.group}.kmc_suf {output.suf}
            touch {output.complex}
        else
            {{
                echo "INPUT:"
                target_prefix=$(echo results/grouped/{wildcards.group}.kmc_pre | sed 's/\\.kmc_pre$//')
                echo "target = $target_prefix -ci1"
                printf '%s\\n' {input.other_dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{print "set" NR " = " $0 " -ci1"}}'
                echo "OUTPUT:"
                printf "results/specific/{wildcards.group}_specific = target"
                printf '%s\\n' {input.other_dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{printf " - set%d", NR}} END{{print ""}}'
            }} > {output.complex}

            kmc_tools -t{threads} complex {output.complex}
        fi
        """

rule kmc_dump_kmers:
    input:
        kmc_db=["results/specific/{group}_specific.kmc_pre", "results/specific/{group}_specific.kmc_suf"]
    output:
        dump="results/specific/{group}.dump",
        fasta="results/specific/{group}_specific.fasta"
    conda: "../envs/kmc.yaml"
    shell:
        """
        db_prefix=$(echo {input.kmc_db[0]} | sed 's/\\.kmc_pre$//')
        kmc_tools -t{threads} transform "$db_prefix" dump {output.dump}
        cat {output.dump} | awk '{{print ">kmer_" NR "\\n" $1}}' > {output.fasta}
        """

rule metaspades:
    input:
        kmers="results/specific/{group}_specific.fasta"
    output:
        assembly="results/assembly/{group}_assembly.fasta"
    conda: "../envs/metaspades.yaml"
    params:
        kmers=config["metaspades_kmers"],
        mem=config["metaspades_mem"],
        tmp_dir=config["metaspades_tmp_dir"]
    shell:
        """
        outdir=$(dirname {output.assembly})/metaspades_{wildcards.group}
        rm -rf "$outdir"
        mkdir -p "$outdir"

        spades.py \
            --meta \
            -m {params.mem} \
            -t {threads} \
            -k {params.kmers} \
            --trusted-contigs {input.kmers} \
            --tmp-dir {params.tmp_dir} \
            -o "$outdir"

        cp "$outdir"/scaffolds.fasta {output.assembly}
        rm -rf "$outdir"
        """
