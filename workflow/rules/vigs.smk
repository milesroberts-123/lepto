# VIGS target design: find regions of target transcripts with no exact
# >=vigs_k nt match to any other transcript in the transcriptome
# (config["vigs_transcriptome_fasta"], guidelines doc: reject any target
# containing a >=19-21 nt contiguous match with a non-target
# transcript). Logic: split the transcriptome into target and
# non-target records, count non-target k-mers with KMC (both strands),
# subtract them from the sliding-window k-mers of the targets, then
# assemble the surviving runs into 200-300 bp insert candidate blocks
# passing GC 40-55% and no homopolymer runs > vigs_max_homopolymer.


rule vigs_split:
    # Split the transcriptome into target and non-target records by
    # exact record ID (seqkit grep matches whole IDs, not substrings).
    # Fails loudly if the target list is empty, has duplicates, contains
    # IDs absent from the fasta, or covers every record (empty
    # non-target set).
    input:
        config["vigs_transcriptome_fasta"]
    output:
        target=temp("results/vigs/target.fasta"),
        nontarget=temp("results/vigs/nontarget.fasta")
    log:
        "logs/vigs/split.log"
    threads: 2
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
        seqkit grep -f results/vigs/.target_ids.txt -v {input} > {output.nontarget} 2>> {log}
        rm results/vigs/.target_ids.txt
        nt=$(seqkit seq -n -i {output.nontarget} | wc -l)
        ng=$(seqkit seq -n -i {output.target} | wc -l)
        if [ "$nt" -eq 0 ]; then
            echo "all transcriptome records are targets: nothing to compare against" >&2
            exit 1
        fi
        if [ "$ng" -ne "$n" ]; then
            echo "targets missing from the transcriptome fasta:" >&2
            printf '%s\\n' {params.ids} | sort > results/vigs/.ids_sorted.txt
            seqkit seq -n -i {output.target} | sort > results/vigs/.targets_sorted.txt
            comm -23 results/vigs/.ids_sorted.txt results/vigs/.targets_sorted.txt >&2
            rm results/vigs/.ids_sorted.txt results/vigs/.targets_sorted.txt
            exit 1
        fi
        """


rule vigs_nontarget_rc:
    # Reverse complement of the non-target records: the KMC k-mer set
    # must cover both strands.
    input:
        "results/vigs/nontarget.fasta"
    output:
        temp("results/vigs/nontarget_rc.fasta")
    log:
        "logs/vigs/nontarget_rc.log"
    threads: 2
    conda: "../envs/vigs.yaml"
    shell:
        "seqkit seq -r {input} > {output} 2> {log}"


rule vigs_nontarget_kmers:
    # KMC set of every non-target k-mer, forward and reverse strand,
    # dumped as one k-mer per line.
    input:
        fwd="results/vigs/nontarget.fasta",
        rev="results/vigs/nontarget_rc.fasta"
    output:
        "results/vigs/nontarget_kmers.txt"
    log:
        "logs/vigs/nontarget_kmers.log"
    threads: 8
    resources:
        mem_mb=12000
    conda: "../envs/kmc.yaml"
    params:
        k=config["vigs_k"]
    shell:
        """
        mkdir -p results/vigs logs/vigs
        printf '%s\\n' {input.fwd} {input.rev} > results/vigs/.kmc_inputs.txt
        local_tmp=.snakemake/tmp_kmc_vigs
        mkdir -p $local_tmp
        # KMC option flags must precede the @input-list argument; flags
        # after the working directory are silently ignored (verified on
        # kmc 3.2.4: -k19 after @file ran as k=25).
        kmc -k{params.k} -m8 -sm -t{threads} -ci1 -fm \\
            @results/vigs/.kmc_inputs.txt results/vigs/.nontarget_db $local_tmp \\
            > {log} 2>&1
        kmc_tools -t{threads} transform results/vigs/.nontarget_db dump \\
            results/vigs/.nontarget_dump.txt >> {log} 2>&1
        cut -f1 results/vigs/.nontarget_dump.txt > {output}
        rm -f results/vigs/.kmc_inputs.txt results/vigs/.nontarget_db.kmc_pre \\
            results/vigs/.nontarget_db.kmc_suf results/vigs/.nontarget_dump.txt
        """


rule vigs_target_seq:
    # Target sequences as id<TAB>seq, for GC/homopolymer filters and
    # block extraction.
    input:
        "results/vigs/target.fasta"
    output:
        temp("results/vigs/target_seq.tsv")
    log:
        "logs/vigs/target_seq.log"
    threads: 2
    conda: "../envs/vigs.yaml"
    shell:
        "mkdir -p results/vigs logs/vigs; seqkit fx2tab -i {input} > {output} 2> {log}"


rule vigs_target_kmers:
    # Sliding k-mers of every target record (fx2tab adds a trailing
    # empty column after the sequence), filtered to ACGT k-mers that
    # occur in neither orientation in the non-target set. Columns:
    # target, 1-based start, end, kmer.
    input:
        kmers_src="results/vigs/target.fasta",
        set="results/vigs/nontarget_kmers.txt"
    output:
        temp("results/vigs/target_kmers.tsv")
    log:
        "logs/vigs/target_kmers.log"
    threads: 4
    conda: "../envs/vigs.yaml"
    params:
        k=config["vigs_k"]
    shell:
        """
        mkdir -p results/vigs logs/vigs
        seqkit sliding -W {params.k} -s 1 {input.kmers_src} \\
            | seqkit fx2tab -i \\
            | awk -F'\\t' -v OFS='\\t' -v k={params.k} \\
                -v setfile={input.set} '
                BEGIN {{
                    while ((getline s < setfile) > 0) set[s] = 1
                }}
                {{
                    id = $1; sub(/_sliding:.*/, "", id)
                    pos = $1; sub(/^.*_sliding:/, "", pos)
                    sub(/-.*/, "", pos)
                    seq = $2
                    if (seq !~ /^[ACGT]+$/) next
                    rc = ""
                    for (i = length(seq); i >= 1; i--) {{
                        c = substr(seq, i, 1)
                        rc = rc ((c == "A") ? "T" : (c == "C") ? "G" : (c == "G") ? "C" : "A")
                    }}
                    if (seq in set || rc in set) next
                    end = pos + k - 1
                    print id, pos, end, seq
                }}' > {output} 2> {log}
        """


rule vigs_assemble:
    # Cluster surviving k-mers into maximal runs of consecutive
    # positions: a run of L k-mers at positions p..p+L-1 is the unique
    # stretch [p, p+L+k-2] (every k-mer fully inside it survives, and
    # the flanking k-mers fail only because they hit a non-target
    # transcript). Each stretch yields one recommended insert block:
    # the whole stretch if 200-300 bp (vigs_insert_min/max), else the
    # earliest longest <=vigs_insert_max window >=vigs_insert_min
    # passing GC and homopolymer filters. No block is chosen for
    # stretches < vigs_insert_min.
    input:
        seq="results/vigs/target_seq.tsv",
        kmers="results/vigs/target_kmers.tsv"
    output:
        blocks="results/vigs/vigs_targets.tsv",
        fasta="results/vigs/vigs_target_blocks.fasta",
        unitigs="results/vigs/vigs_unitigs.tsv"
    log:
        "logs/vigs/assemble.log"
    threads: 1
    conda: "../envs/vigs.yaml"
    params:
        k=config["vigs_k"],
        imin=config["vigs_insert_min"],
        imax=config["vigs_insert_max"],
        mingc=config["vigs_min_gc"],
        maxgc=config["vigs_max_gc"],
        maxhp=config["vigs_max_homopolymer"]
    shell:
        """
        mkdir -p results/vigs logs/vigs
        awk -F'\\t' -v OFS='\\t' \\
            -v k={params.k} -v imin={params.imin} -v imax={params.imax} \\
            -v mingc={params.mingc} -v maxgc={params.maxgc} -v maxhp={params.maxhp} \\
            -v seqfile={input.seq} '
            function maxrun(s,   n, i, m) {{
                m = 1; n = 1
                for (i = 2; i <= length(s); i++) {{
                    n = (substr(s, i, 1) == substr(s, i - 1, 1)) ? n + 1 : 1
                    if (n > m) m = n
                }}
                return m
            }}
            function gcf(s,   L, n) {{
                L = length(s)
                n = gsub(/[GC]/, "", s)
                return n / L
            }}
            function pass(s) {{
                return gcf(s) >= mingc && gcf(s) <= maxgc && maxrun(s) <= maxhp
            }}
            function emit_block(   L, seq, gc, W, j, cand, blk, bs, be) {{
                L = run_end - run_start + 1
                seq = substr(seqs[run_id], run_start, L)
                gc = gcf(seq)
                print run_id, run_start, run_end, L, sprintf("%.3f", gc) >> unitigs_path
                blk = ""; bs = 0; be = 0
                if (L >= imin) {{
                    for (W = imax; W >= imin && blk == ""; W--) {{
                        if (W > L) continue
                        for (j = run_start; j <= run_start + L - W && blk == ""; j++) {{
                            cand = substr(seq, j - run_start + 1, W)
                            if (pass(cand)) {{ blk = cand; bs = j; be = j + W - 1 }}
                        }}
                    }}
                }}
                if (blk != "") {{
                    print run_id, bs, be, length(blk), sprintf("%.3f", gcf(blk)) >> blocks_path
                    printf ">%s_p%d-%d_len%d_gc%.3f\\n%s\\n", run_id, bs, be, length(blk), gcf(blk), blk >> fasta_path
                }}
            }}
            BEGIN {{
                while ((getline line < seqfile) > 0) {{
                    split(line, a, "\\t"); seqs[a[1]] = a[2]
                }}
                unitigs_path = "{output.unitigs}"
                blocks_path = "{output.blocks}"
                fasta_path = "{output.fasta}"
                print "target\\tstretch_start\\tstretch_end\\tlength\\tgc" >> unitigs_path
                print "target\\tblock_start\\tblock_end\\tlength\\tgc" >> blocks_path
                run_id = ""; run_start = 0; run_end = -1
            }}
            {{
                if ($1 != run_id || $2 != prev_start + 1) {{
                    if (run_id != "" && run_end >= run_start) emit_block()
                    run_id = $1; run_start = $2
                }}
                run_end = $3; prev_start = $2
            }}
            END {{
                if (run_id != "" && run_end >= run_start) emit_block()
            }}' 2> {log}
        """