# VIGS target design: find stretches of consecutive target k-mers that
# occur nowhere else in the transcriptome (config["vigs_transcriptome_fasta"],
# guidelines doc: reject any target containing a >=19-21 nt contiguous
# match with a non-target transcript, either orientation). Three rules:
# isolate the chosen target records, count their distinct k-mers with
# KMC, then map those k-mers back to the whole transcriptome with
# bwa mem and keep the stretches whose every k-mer aligns only to
# chosen target records. bwa checks both strands and, with -T k on
# k-nt queries, an alignment scores >= k only as a full-length exact
# match; -a reports every locus. Stretches of consecutive surviving
# k-mers land in results/vigs/vigs_stretches.tsv.


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


rule vigs_map:
    # bwa mem over the whole transcriptome with -a (-T k: on k-nt
    # queries any mismatch drops the score below threshold, so only
    # full-length exact alignments are reported; -a reports every
    # locus, both strands). A k-mer survives when its complete set of
    # hit loci lies inside the chosen target records. Surviving k-mers
    # are located back onto the target records (either orientation) and
    # compressed into maximal runs of consecutive k-mer positions -
    # one tsv row per stretch.
    input:
        ref=config["vigs_transcriptome_fasta"],
        target="results/vigs/target.fasta",
        kmers="results/vigs/kmers.fasta"
    output:
        stretches="results/vigs/vigs_stretches.tsv"
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
            results/vigs/.idx {input.kmers} > results/vigs/.kmers.sam 2>> {log}
        # header first, rows appended by the awk's END block
        printf 'target\\tstretch_start\\tstretch_end\\tlength\\tgc\\n' > {output.stretches}
        awk -F'\\t' -v OFS='\\t' -v k={params.k} -v out_path={output.stretches} '
            function revcomp(s,   i, c, o) {{
                o = ""
                for (i = length(s); i >= 1; i--) {{
                    c = substr(s, i, 1)
                    o = o ((c == "A") ? "T" : (c == "C") ? "G" : (c == "G") ? "C" : (c == "T") ? "A" : c)
                }}
                return o
            }}
            function gcf(s,   L, n) {{
                L = length(s)
                n = gsub(/[GC]/, "", s)
                return n / L
            }}
            # file 1: the chosen target fasta (record IDs + sequences)
            FNR == NR {{
                if ($0 ~ /^>/) {{
                    id = substr($0, 2)
                    ids[id] = 1
                    nt++; order[nt] = id
                }} else {{
                    seqs[id] = seqs[id] $0
                }}
                next
            }}
            # file 2: the bwa mem -a SAM (qname = the k-mer sequence)
            {{
                if ($1 ~ /^@/) next
                if ($3 == "*" || $3 == "") next
                seen[$1] = 1
                hits[$1 SUBSEP $3] = 1
            }}
            END {{
                if (!nt) {{
                    print "target fasta has no records" > "/dev/stderr"
                    exit 1
                }}
                # a k-mer survives when no reported locus falls outside
                # the chosen target records
                for (key in hits) {{
                    split(key, p, SUBSEP)
                    if (p[2] in ids) ontarget[p[1]] = 1
                    else off[p[1]] = 1
                }}
                for (q in seen) {{
                    if (ontarget[q] && !off[q]) kept[q] = 1
                }}
                # per target, ascending pass: a window counts when it
                # hosts a kept k-mer in either orientation; gaps split
                # the stretches
                for (i = 1; i <= nt; i++) {{
                    t = order[i]; s = seqs[t]
                    L = length(s)
                    if (L < k) continue
                    run = 0; rs = 0; re = 0
                    for (ps = 1; ps <= L - k + 1; ps++) {{
                        x = substr(s, ps, k)
                        if ((x in kept) || (revcomp(x) in kept)) {{
                            if (!run) {{ rs = ps; run = 1 }}
                            re = ps + k - 1
                        }} else if (run) {{
                            seq2 = substr(s, rs, re - rs + 1)
                            print t, rs, re, re - rs + 1, sprintf("%.3f", gcf(seq2)) >> out_path
                            run = 0
                        }}
                    }}
                    if (run) {{
                        seq2 = substr(s, rs, re - rs + 1)
                        print t, rs, re, re - rs + 1, sprintf("%.3f", gcf(seq2)) >> out_path
                    }}
                }}
            }}' {input.target} {input.sam} 2> {log}
        rm results/vigs/.kmers.sam results/vigs/.idx.amb results/vigs/.idx.ann \\
            results/vigs/.idx.bwt results/vigs/.idx.pac results/vigs/.idx.sa
        """