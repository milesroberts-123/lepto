# k-mer counting, group k-mer sets, and read mapping to references.
# Chain: CRAM -> fastq -> fastp -> kmc -> group intersect/union/subtract
#        -> filtered kmers -> BWA mem -> bam/bed -> minimap2 ref alignment.

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
        pre=temp("results/intersect/{group}.kmc_pre"),
        suf=temp("results/intersect/{group}.kmc_suf"),
        complex=temp("results/intersect/{group}.complex")
    conda: "../envs/kmc.yaml"
    params:
        depth=config.get("kmc_intersect_depth", 1000)
    shell:
        """
        mkdir -p results/grouped

        {{
            echo "INPUT:"
            printf '%s\\n' {input.dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{print "set" NR " = " $0 " -ci1"}}'
            echo "OUTPUT:"
            printf "results/intersect/{wildcards.group} = "
            printf '%s\\n' {input.dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{printf "%sset%d", (NR>1?" * ":""), NR}} END{{print ""}}'
            echo "OUTPUT_PARAMS:"
            echo "-cs{params.depth}"
        }} > {output.complex}

        kmc_tools -t{threads} complex {output.complex}
        """

rule kmc_union_group:
    input:
        dbs=expand(["results/kmc/{ID}/kmc_db.kmc_pre", "results/kmc/{ID}/kmc_db.kmc_suf"],
                   ID=lambda wildcards: samples_by_group[wildcards.group])
    output:
        pre=temp("results/union/{group}.kmc_pre"),
        suf=temp("results/union/{group}.kmc_suf"),
        complex=temp("results/union/{group}.complex")
    conda: "../envs/kmc.yaml"
    params:
        depth=config.get("kmc_intersect_depth", 1000)
    shell:
        """
        mkdir -p results/grouped

        {{
            echo "INPUT:"
            printf '%s\\n' {input.dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{print "set" NR " = " $0 " -ci1"}}'
            echo "OUTPUT:"
            printf "results/union/{wildcards.group} = "
            printf '%s\\n' {input.dbs} | grep '\\.kmc_pre$' | sed 's/\\.kmc_pre$//' | awk '{{printf "%sset%d", (NR>1?" + ":""), NR}} END{{print ""}}'
            echo "OUTPUT_PARAMS:"
            echo "-cs{params.depth}"
        }} > {output.complex}

        kmc_tools -t{threads} complex {output.complex}
        """


rule kmc_subtract:
    input:
        target_db=["results/intersect/{group}.kmc_pre", "results/intersect/{group}.kmc_suf"],
        other_dbs=expand(["results/union/{other}.kmc_pre", "results/union/{other}.kmc_suf"],
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
            cp results/intersect/{wildcards.group}.kmc_pre {output.pre}
            cp results/intersect/{wildcards.group}.kmc_suf {output.suf}
            touch {output.complex}
        else
            {{
                echo "INPUT:"
                target_prefix=$(echo results/intersect/{wildcards.group}.kmc_pre | sed 's/\\.kmc_pre$//')
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
        temp("results/specific/{group}.dump")
    conda: "../envs/kmc.yaml"
    shell:
        """
        db_prefix=$(echo {input.kmc_db[0]} | sed 's/\\.kmc_pre$//')
        kmc_tools -t{threads} transform "$db_prefix" dump {output}
        """

rule filter_dump:
    input:
        "results/specific/{group}.dump"
    output:
        temp("results/specific/{group}_40.dump")
    conda: "../envs/kmc.yaml"
    params:
        threshold=config["bwa_count_threshold"]
    shell:
        "awk '($2 >= {params.threshold})' {input} > {output}"

rule dump_to_fasta:
    input:
        "results/specific/{group}_40.dump"
    output:
        temp("results/specific/{group}_40.fasta")
    conda: "../envs/kmc.yaml"
    shell:
        "awk '{{print \">kmer_\" NR \"\\n\" $1}}' {input} > {output}"

rule bwa_index:
    input:
        "../resources/reference/{ref_name}.fasta"
    output:
        amb="../resources/reference/{ref_name}.fasta.amb",
        ann="../resources/reference/{ref_name}.fasta.ann",
        bwt="../resources/reference/{ref_name}.fasta.bwt",
        pac="../resources/reference/{ref_name}.fasta.pac",
        sa="../resources/reference/{ref_name}.fasta.sa"
    conda: "../envs/bwa.yaml"
    shell:
        "bwa index {input}"

rule bwa_mem:
    input:
        fasta="results/specific/{group}_40.fasta",
        ref="../resources/reference/{ref_name}.fasta",
        index=rules.bwa_index.output
    output:
        "results/bwa/{ref_name}/{group}.bam"
    conda: "../envs/bwa.yaml"
    params:
        k=config["bwa_k"],
        A=config["bwa_A"],
        B=config["bwa_B"],
        O=config["bwa_O"],
        E=config["bwa_E"],
        L=config["bwa_L"],
        T=config["bwa_T"]
    shell:
        "bwa mem -t {threads} -k {params.k} -A {params.A} -B {params.B} -O {params.O} -E {params.E} -L {params.L} -T {params.T} -a {input.ref} {input.fasta} | samtools sort -@ {threads} -o {output}"

rule samtools_index:
    input:
        "results/bwa/{ref_name}/{group}.bam"
    output:
        "results/bwa/{ref_name}/{group}.bam.bai"
    conda: "../envs/bwa.yaml"
    shell:
        "samtools index {input}"

rule bam_to_bed:
    input:
        "results/bwa/{ref_name}/{group}.bam"
    output:
        "results/bwa/{ref_name}/{group}.bed"
    conda: "../envs/bwa.yaml"
    params:
        merge_distance=config["bwa_merge_distance"]
    shell:
        "samtools view -b -F 4 {input} | bedtools bamtobed -i - | sort -k1,1 -k2,2n | bedtools merge -d {params.merge_distance} -i - > {output}"

rule minimap2_refs:
    input:
        query="../resources/reference/{query}.fasta",
        target="../resources/reference/{target}.fasta"
    output:
        "results/minimap2/{query}_vs_{target}.paf"
    conda: "../envs/minimap2.yaml"
    params:
        f=config["minimap2_f"]
    shell:
        "minimap2 -t {threads} -f {params.f} {input.target} {input.query} > {output}"

