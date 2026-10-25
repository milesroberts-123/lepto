# VIGS target design: map the targets' k-mers back to the whole
# transcriptome (config["vigs_transcriptome_fasta"], guidelines doc:
# reject any target containing a >=19-21 nt contiguous match with a
# non-target transcript, either orientation). Three rules: isolate the
# chosen target records, count their distinct k-mers with KMC, then
# map those k-mers back with bwa mem and stop at the coordinate-sorted,
# indexed BAM. Downstream analysis (which k-mers only hit the chosen
# targets, stretch assembly) happens outside the workflow.


rule vigs_target:
    # Isolate the chosen target records by exact record ID (seqkit grep
    # matches whole IDs, not substrings). Fails loudly if the target
    # list is empty, has duplicates, or contains IDs absent from the
    # transcriptome fasta.
    input:
        config["vigs_transcriptome_fasta"]
    output:
        target=temp("results/vigs/target.fasta")
    log:
        "logs/vigs/target.log"
    conda: "../envs/vigs.yaml"
    params:
        # Space-joined so a multi-target list stays on one shell line;
        # printf '%s\n' expands it to one ID per line again.
        ids=lambda wc: " ".join(config["vigs_target_transcripts"])
    shell:
        """
        mkdir -p results/vigs logs/vigs
        n=$(printf '%s\\n' {params.ids} | grep -c .)
        if [ "$n" -eq 0 ]; then
            echo "vigs_target_transcripts is empty: set record IDs in config.yaml" >&2
            exit 1
        fi
        if printf '%s\\n' {params.ids} | sort | uniq -d | grep .; then
            echo "vigs_target_transcripts contains duplicate IDs" >&2
            exit 1
        fi
        printf '%s\\n' {params.ids} > results/vigs/.target_ids.txt
        seqkit grep -f results/vigs/.target_ids.txt {input} > {output.target} 2> {log}
        rm results/vigs/.target_ids.txt
        ng=$(seqkit seq -n -i {output.target} | wc -l)
        if [ "$ng" -ne "$n" ]; then
            echo "targets missing from the transcriptome fasta:" >&2
            printf '%s\\n' {params.ids} | sort > results/vigs/.ids_sorted.txt
            seqkit seq -n -i {output.target} | sort > results/vigs/.targets_sorted.txt
            comm -23 results/vigs/.ids_sorted.txt results/vigs/.targets_sorted.txt >&2
            rm results/vigs/.ids_sorted.txt results/vigs/.targets_sorted.txt
            exit 1
        fi
        """


rule vigs_kmers:
    # Count the target's k-mers with KMC (fasta input needs -fm; option
    # flags must precede the file arguments, KMC silently ignores them
    # otherwise) and dump the distinct k-mers as a FASTA whose read
    # names are the k-mer sequences themselves - bwa mem's QNAME then
    # keys every alignment back to its k-mer.
    input:
        "results/vigs/target.fasta"
    output:
        kmers=temp("results/vigs/kmers.fasta")
    log:
        "logs/vigs/kmers.log"
    conda: "../envs/kmc.yaml"
    params:
        k=config["vigs_k"]
    shell:
        """
        mkdir -p results/vigs logs/vigs
        rm -rf results/vigs/.kmc_tmp && mkdir -p results/vigs/.kmc_tmp
        kmc -k{params.k} -m8 -sm -t{threads} -ci1 -fm \\
            {input} results/vigs/.kmers_db results/vigs/.kmc_tmp > {log} 2>&1
        kmc_tools -t{threads} transform results/vigs/.kmers_db dump \\
            results/vigs/.kmers_dump.txt >> {log} 2>&1
        awk -F'\\t' '{{print ">" $1 "\\n" $1}}' results/vigs/.kmers_dump.txt > {output.kmers}
        rm -f results/vigs/.kmers_db.kmc_pre results/vigs/.kmers_db.kmc_suf \\
            results/vigs/.kmers_dump.txt
        rm -rf results/vigs/.kmc_tmp
        """


rule vigs_locate:
    input:
        trans=config["vigs_transcriptome_fasta"],
        kmers="results/vigs/kmers.fasta"
    output:
        "results/vigs/match_list.txt"
    conda: "../envs/vigs.yaml"
    shell:
        "seqkit locate --threads {threads} -f {input.kmers} {input.trans} > {output}"

rule vigs_map:
    # bwa mem over the whole transcriptome; -T k makes every reported
    # alignment of a k-nt query a full-length exact match (any mismatch
    # scores below the threshold) and -a reports every locus as its own
    # record, reverse-strand hits included. samtools sort + index make
    # a standard coordinate-sorted, indexed BAM for downstream analysis;
    # each k-mer's records are reunifiable because the QNAME is the
    # k-mer sequence itself. The index prefix keeps bwa files out of
    # resources/ and is removed after the run.
    input:
        ref=config["vigs_transcriptome_fasta"],
        target="results/vigs/target.fasta",
        kmers="results/vigs/kmers.fasta"
    output:
        bam="results/vigs/kmers_mapped.bam",
        bai="results/vigs/kmers_mapped.bam.bai"
    log:
        "logs/vigs/map.log"
    conda: "../envs/vigs.yaml"
    params:
        k=config["vigs_k"]
    shell:
        """
        mkdir -p results/vigs logs/vigs
        bwa index -p results/vigs/.idx {input.ref} > {log} 2>&1
        bwa mem -a -T {params.k} -k {params.k} -t {threads} \\
            results/vigs/.idx {input.kmers} 2>> {log} \\
            | samtools sort -o {output.bam}
        samtools index {output.bam}
        rm results/vigs/.idx.amb results/vigs/.idx.ann \\
            results/vigs/.idx.bwt results/vigs/.idx.pac results/vigs/.idx.sa
        """
