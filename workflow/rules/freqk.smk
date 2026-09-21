rule freqk_index:
    input:
        vcf="slim_results/norm_{SID}_{PID}.vcf.gz",
        fasta="ancestral_genome_results/{SID}.fasta",
        fai="ancestral_genome_results/{SID}.fasta.fai",
    output:
        index=temp("results/{SID}_{PID}.txt"),
    benchmark:
        "benchmarks/freqk_index/{SID}_{PID}.bench"
    group: "freqk_index"
    params:
        k=lookup(query="ID == '{SID}'", within=parameters, cols="k"),
    shell:
        """
        # index panel of variants
        ./scripts/freqk index --fasta {input.fasta} --vcf {input.vcf} -k {params.k} --output {output.index}
        """

rule freqk_var_dedup:
    input:
        "freqk_indices/{SID}_{PID}.txt",
    output:
        temp("freqk_var_dedup/{SID}_{PID}.txt"),
    benchmark:
        "benchmarks/freqk_var_dedup/{SID}_{PID}.bench"
    group: "freqk_index"
    shell:
        """
        ./scripts/freqk var-dedup --index {input} --output {output}
        """

rule freqk_ref_dedup:
    input:
        index="freqk_var_dedup/{SID}_{PID}.txt",
        vcf="slim_results/norm_{SID}_{PID}.vcf.gz",
        fasta="ancestral_genome_results/{SID}.fasta",
        fai="ancestral_genome_results/{SID}.fasta.fai",
    output:
        "freqk_ref_dedup/{SID}_{PID}.txt",
    benchmark:
        "benchmarks/freqk_ref_dedup/{SID}_{PID}.bench"
    group: "freqk_index"
    shell:
        """
        ./scripts/freqk ref-dedup --index {input.index} --fasta {input.fasta} --vcf {input.vcf} --output {output}
        """

rule freqk_count:
    input:
        reads="all_{SID}_{PID}.fastq",
        index="freqk_ref_dedup/{SID}_{PID}.txt",
    output:
        counts="freqk_results/{SID}_{PID}_counts.txt",
        freqs="freqk_results/{SID}_{PID}_freqs.txt",
    group: "freqk_count"
    benchmark:
        "benchmarks/freqk_count/{SID}_{PID}.bench"
    shell:
        """
        ./scripts/freqk count --nthreads {threads} --index {input.index} --reads {input.reads} --freq-output {output.freqs} --count-output {output.counts}
        """

rule freqk_call:
    input:
        counts="freqk_results/{SID}_{PID}_freqs.txt",
        index="freqk_ref_dedup/{SID}_{PID}.txt",
    output:
        "freqk_results/{SID}_{PID}_calls.txt",
    group: "freqk_count"
    priority: 1000
    benchmark:
        "benchmarks/freqk_call/{SID}_{PID}.bench"
    shell:
        "./scripts/freqk call --index {input.index} -c {input.counts} --output {output}"
