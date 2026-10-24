# VIGS target design: find regions of target transcripts with no exact
# >=vigs_k nt match to any other transcript in the transcriptome
# (config["vigs_transcriptome_fasta"], guidelines doc: reject any target
# containing a >=19-21 nt contiguous match with a non-target
# transcript). Logic: split the transcriptome into target and
# non-target records; for EACH target build a Bifrost coloured de
# Bruijn graph over just that transcript plus the non-target records
# with unitig-caller --call (canonical k-mers cover both strands; the
# per-target graph prevents cross-transcript graph glue from stealing
# another target's attribution); keep the unitigs present only in the
# target colour; probe their endpoint k-mers against the non-target
# set with unitig-caller --simple (the tool's head/tail colour
# intersection can over-claim shared endpoint k-mers, so any endpoint
# k-mer found in the non-target transcripts costs its 19 bp); locate
# the trimmed unitigs back onto the target transcript, then cut the
# surviving stretches into 200-300 bp insert candidate blocks passing
# GC 40-55% and no homopolymer runs > vigs_max_homopolymer.


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


rule vigs_target_slice:
    # One target record (id + fasta copy) on its own; {idx} indexes
    # config["vigs_target_transcripts"]. Each index gets its own graph
    # in vigs_unitig_caller so unitigs cannot glue across target
    # transcripts. The rtab colour column for the slice is basenamed
    # target_{idx}.
    input:
        "results/vigs/target.fasta"
    output:
        fasta=temp("results/vigs/target_{idx}.fasta"),
        seq=temp("results/vigs/target_{idx}.seq.tsv")
    log:
        "logs/vigs/slice_{idx}.log"
    conda: "../envs/vigs.yaml"
    params:
        name=lambda wc: config["vigs_target_transcripts"][int(wc.idx) - 1],
        ntargets=lambda wc: len(config["vigs_target_transcripts"])
    shell:
        """
        mkdir -p results/vigs logs/vigs
        if [ "{wildcards.idx}" -lt 1 ] || [ "{wildcards.idx}" -gt "{params.ntargets}" ]; then
            echo "target index {wildcards.idx} out of range 1..{params.ntargets}" >&2
            exit 1
        fi
        seqkit grep -p '{params.name}' {input} > {output.fasta} 2> {log}
        if [ ! -s {output.fasta} ]; then
            echo "target transcript '{params.name}' not found in the transcriptome fasta" >&2
            exit 1
        fi
        recs=$(seqkit seq -n -i {output.fasta} | wc -l)
        if [ "$recs" -ne 1 ]; then
            echo "record ID '{params.name}' matches $recs records (duplicate IDs in the transcriptome fasta)" >&2
            exit 1
        fi
        seqkit fx2tab -i {output.fasta} > {output.seq} 2>> {log}
        """


rule vigs_unitig_caller:
    # unitig-caller --call builds a Bifrost compacted coloured de Bruijn
    # graph over this target's slice plus the non-target records
    # (canonical k-mers: both strands handled, no separate RC copy) and
    # reports every unitig with its colours. Colours are per input
    # FILE: the .rtab columns are the ref-file basenames minus
    # extension in refs.txt order (target_{idx}, nontarget); rows are
    # unitig sequence and 1/0 per colour. The head/tail colour
    # intersection keeps the claim target-only only when both path ends
    # agree; interior k-mers of a simple path cannot be shared (a
    # shared k-mer's off-path edges break the path at it), but shared
    # END k-mers get intersected away - the probe rule re-checks them.
    input:
        slice="results/vigs/target_{idx}.fasta",
        nontarget="results/vigs/nontarget.fasta"
    output:
        rtab="results/vigs/unitig_caller_{idx}.rtab"
    log:
        "logs/vigs/unitig_caller_{idx}.log"
    conda: "../envs/unitig_caller.yaml"
    params:
        k=config["vigs_k"]
    shell:
        """
        mkdir -p results/vigs logs/vigs
        printf '%s\\n' {input.slice} {input.nontarget} > results/vigs/.refs_{wildcards.idx}.txt
        unitig-caller --call --refs results/vigs/.refs_{wildcards.idx}.txt \\
            --kmer {params.k} --rtab --threads {threads} \\
            --out results/vigs/.caller_{wildcards.idx} > {log} 2>&1
        mv results/vigs/.caller_{wildcards.idx}.rtab {output.rtab}
        rm results/vigs/.refs_{wildcards.idx}.txt
        """


rule vigs_probe:
    # unitig-caller's --simple mode (FM-index exact matching, both
    # strands) over the non-target fasta, queried with both endpoint
    # k-mers of every target-only unitig from the rtab. --simple omits
    # queries found nowhere, so the output is a SET of non-target-shared
    # endpoint k-mers; the map rule trims 19 bp behind any member.
    input:
        rtab="results/vigs/unitig_caller_{idx}.rtab",
        nontarget="results/vigs/nontarget.fasta"
    output:
        probe="results/vigs/probe_{idx}.pyseer"
    log:
        "logs/vigs/probe_{idx}.log"
    conda: "../envs/unitig_caller.yaml"
    params:
        k=config["vigs_k"]
    shell:
        """
        mkdir -p results/vigs logs/vigs
        awk -F'\\t' -v OFS='\\t' -v k={params.k} '
            NR == 1 {{
                for (c = 2; c <= NF; c++) {{
                    if ($c == "nontarget") ncol = c
                    else if ($c != "Unitig_sequence") tcol = c
                }}
                if (!tcol || !ncol) {{
                    print "rtab header lacks expected colour columns" > "/dev/stderr"
                    exit 1
                }}
                print "probe"
                next
            }}
            {{
                u = $1
                if (u !~ /^[ACGT]+$/) next
                if ($tcol != 1 || $ncol != 0) next
                print substr(u, 1, k)
                print substr(u, length(u) - k + 1, k)
            }}' {input.rtab} > results/vigs/.probe_{wildcards.idx}.txt
        printf '%s\\n' {input.nontarget} > results/vigs/.probe_refs_{wildcards.idx}.txt
        unitig-caller --simple --refs results/vigs/.probe_refs_{wildcards.idx}.txt \\
            --unitigs results/vigs/.probe_{wildcards.idx}.txt --no-save-idx \\
            --threads {threads} --out results/vigs/.probe_{wildcards.idx} > {log} 2>&1
        mv results/vigs/.probe_{wildcards.idx}.pyseer {output.probe}
        rm results/vigs/.probe_{wildcards.idx}.txt results/vigs/.probe_refs_{wildcards.idx}.txt
        """


rule vigs_unitigs_map:
    # Trim target-only unitigs whose endpoint k-mer is non-target-shared
    # (both the head and the tail k-mer are checked against the probe
    # set, 19 bp each), locate each trimmed unitig - in whichever
    # orientation occurs - onto its target transcript (all
    # occurrences), and emit one unitig row per occurrence plus insert
    # candidate blocks: the whole occurrence if 200-300 bp
    # (vigs_insert_min/max), else the earliest longest
    # <=vigs_insert_max window >=vigs_insert_min passing GC and
    # homopolymer filters. Occurrences < vigs_insert_min yield no
    # block. The per-target graph should make every target-only unitig
    # map somewhere; a residual warning on stderr keeps silent loss
    # visible (possible for unitigs glued across inverted-repeat
    # junctions within one transcript).
    input:
        rtab="results/vigs/unitig_caller_{idx}.rtab",
        probe="results/vigs/probe_{idx}.pyseer",
        seq="results/vigs/target_{idx}.seq.tsv"
    output:
        unitigs=temp("results/vigs/unitigs_{idx}.tsv"),
        blocks=temp("results/vigs/targets_{idx}.tsv"),
        fasta=temp("results/vigs/blocks_{idx}.fasta")
    log:
        "logs/vigs/unitigs_map_{idx}.log"
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
            -v k={params.k} \\
            -v imin={params.imin} -v imax={params.imax} \\
            -v mingc={params.mingc} -v maxgc={params.maxgc} -v maxhp={params.maxhp} \\
            -v seqfile={input.seq} -v probefile={input.probe} '
            function revcomp(s,   i, c, o) {{
                o = ""
                for (i = length(s); i >= 1; i--) {{
                    c = substr(s, i, 1)
                    o = o ((c == "A") ? "T" : (c == "C") ? "G" : (c == "G") ? "C" : (c == "T") ? "A" : c)
                }}
                return o
            }}
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
            function emit_block(s2,   L, W, j, cand, blk, bs, be) {{
                L = bend - bstart + 1
                blk = ""; bs = 0; be = 0
                if (L >= imin) {{
                    for (W = imax; W >= imin && blk == ""; W--) {{
                        if (W > L) continue
                        for (j = bstart; j <= bstart + L - W && blk == ""; j++) {{
                            cand = substr(s2, j - bstart + 1, W)
                            if (pass(cand)) {{ blk = cand; bs = j; be = j + W - 1 }}
                        }}
                    }}
                }}
                if (blk != "") {{
                    print tid, bs, be, length(blk), sprintf("%.3f", gcf(blk)) >> blocks_path
                    printf ">%s_p%d-%d_len%d_gc%.3f\\n%s\\n", tid, bs, be, length(blk), gcf(blk), blk >> fasta_path
                }}
            }}
            function emit_stretch(s2,   L) {{
                L = bend - bstart + 1
                print tid, bstart, bend, L, sprintf("%.3f", gcf(s2)) >> unitigs_path
                emit_block(s2)
            }}
            function locate(u,   s, p, cnt) {{
                # emit every occurrence of u in the target sequence
                s = tseq; p = index(s, u); cnt = 0
                while (p > 0) {{
                    bstart = p; bend = p + length(u) - 1
                    emit_stretch(u)
                    cnt++
                    # advance past the match start to find overlapping hits
                    p = index(substr(s, p + 1), u) + p
                }}
                return cnt
            }}
            BEGIN {{
                # single-record target sequence
                getline line < seqfile
                split(line, b, "\\t"); tid = b[1]; tseq = b[2]
                # probes present in the non-target transcriptome (a set;
                # --simple omits queries found nowhere, order is not kept)
                while ((getline line < probefile) > 0) {{
                    split(line, p, " "); pset[p[1]] = 1
                }}
                unitigs_path = "{output.unitigs}"
                blocks_path = "{output.blocks}"
                fasta_path = "{output.fasta}"
                tcol = ncol = 0
            }}
            {{
                if ($1 == "Unitig_sequence") {{
                    for (c = 2; c <= NF; c++) {{
                        if ($c == "nontarget") ncol = c
                        else if ($c != "Unitig_sequence") tcol = c
                    }}
                    if (!tcol || !ncol) {{
                        print "rtab header lacks expected colour columns" > "/dev/stderr"
                        exit 1
                    }}
                    next
                }}
                u = $1
                if (u !~ /^[ACGT]+$/) next
                if ($tcol != 1 || $ncol != 0) next
                L = length(u)
                s = (substr(u, 1, k) in pset) ? k + 1 : 1
                e = (substr(u, L - k + 1, k) in pset) ? L - k : L
                if (e - s + 1 < k) next
                u2 = substr(u, s, e - s + 1)
                hits = locate(u2)
                if (hits == 0) hits = locate(revcomp(u2))
                if (hits == 0) {{
                    print "unitig (len " length(u2) " bp) target-only but maps onto no whole transcript (likely glued across inverted-repeat junctions); dropped" > "/dev/stderr"
                }}
            }}
            END {{
                if (!tcol || !ncol) {{
                    print "rtab is empty: unitig-caller produced no unitigs" > "/dev/stderr"
                    exit 1
                }}
            }}' {input.rtab} 2> {log}
        """


rule vigs_concat:
    # Assembly of the per-target intermediates into the final output
    # files (same paths and schemas as before the unitig-caller swap).
    input:
        unitigs_list=expand("results/vigs/unitigs_{idx}.tsv",
                            idx=range(1, len(config["vigs_target_transcripts"]) + 1)),
        blocks_list=expand("results/vigs/targets_{idx}.tsv",
                           idx=range(1, len(config["vigs_target_transcripts"]) + 1)),
        fasta_list=expand("results/vigs/blocks_{idx}.fasta",
                          idx=range(1, len(config["vigs_target_transcripts"]) + 1)),
    output:
        unitigs="results/vigs/vigs_unitigs.tsv",
        blocks="results/vigs/vigs_targets.tsv",
        fasta="results/vigs/vigs_target_blocks.fasta"
    log:
        "logs/vigs/concat.log"
    conda: "../envs/vigs.yaml"
    params:
        # Same empty-config guard as vigs_split: an empty target list
        # must fail loudly here, since only this rule runs in that case.
        ids=lambda wc: " ".join(config["vigs_target_transcripts"])
    shell:
        """
        mkdir -p results/vigs logs/vigs
        n=$(printf '%s\\n' {params.ids} | grep -c .)
        if [ "$n" -eq 0 ]; then
            echo "vigs_target_transcripts is empty: set record IDs in config.yaml" >&2
            exit 1
        fi
        printf 'target\\tstretch_start\\tstretch_end\\tlength\\tgc\\n' > {output.unitigs}
        cat {input.unitigs_list} >> {output.unitigs}
        printf 'target\\tblock_start\\tblock_end\\tlength\\tgc\\n' > {output.blocks}
        cat {input.blocks_list} >> {output.blocks}
        cat {input.fasta_list} > {output.fasta}
        """