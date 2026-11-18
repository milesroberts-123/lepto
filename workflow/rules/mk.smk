# McDonald-Kreitman (MK) tests per gene: Dn/Ds from codon alignments of
# 1:1 orthologs (ingroup = Leptosiphon hap1, outgroup = Linanthus
# parryae), Pn/Ps from the per-pool grenedalf allele frequencies
# classified by translation (the ALT allele either changes the amino
# acid -> Pn or not -> Ps; 2-fold sites count, no degeneracy annotation).
# Standard MK only: alpha, neutrality index, DoS, Fisher exact test per
# gene per pool plus a pooled table; codeml pairwise dN/dS is reported
# alongside as a per-gene covariate (the MK counts do not use it).
#
# Products consumed from upstream branches:
#   results/orthofinder/proteomes/OrthoFinder   (orthofinder_all; the
#       results themselves sit in OrthoFinder/Results_<date>/)
#   results/featurecounts/hap1_expressed.{faa,gff3} (featurecounts_all)
#   results/braker/linanthus/braker_expressed.{faa,gff3} (braker_all)
#   results/vg/{variant}_vs_{backbone}/{ID}/grenedalf_results_frequency.csv
#       (vg_diversity_all -- run it first, same acceptance as freqk/ldx)
#
# grenedalf frequency column contract (see mk_polymorphism.py): col1
# CHROM, col2 POS, col3 REF, col4 ALT, col6 ALT count, col7 depth,
# col8 ALT frequency.


rule mk_ortholog_pairs:
    # 1:1 orthologs between the ingroup and outgroup proteomes; only .t1
    # transcripts survive (orthofinder proteomes are t1-filtered).
    input:
        directory("results/orthofinder/proteomes/OrthoFinder")
    output:
        "results/mk/ortholog_pairs.tsv"
    conda: "../envs/mk.yaml"
    shell:
        """
        mkdir -p results/mk
        python3 scripts/mk_ortholog_pairs.py \\
            --orthofinder results/orthofinder/proteomes/OrthoFinder \\
            --ingroup {config[mk_ingroup]} \\
            --outgroup {config[mk_outgroup]} \\
            --out {output}
        """


rule mk_cds_extract:
    # CDS per .t1 transcript plus a CDS-position -> genome map; the map
    # lets codon alignments be projected back onto reference coordinates
    # (polymorphism exclusion) and Pn/Ps be classified on the CDS.
    input:
        fasta=lambda wc: config["mk_species_annotations"][wc.species]["fasta"],
        gff3=lambda wc: config["mk_species_annotations"][wc.species]["gff3"],
    output:
        fna="results/mk/cds/{species}.fna",
        map="results/mk/cds/{species}_map.tsv",
    conda: "../envs/mk.yaml"
    params:
        # the ingroup's map must carry the same scaffold prefix grenedalf
        # uses (the vg backbone prefix), the outgroup's stays unprefixed
        chrom_prefix_arg="--chrom-prefix {prefix}".format(
            prefix=config["vg_ref_backbone"] + "_")
        if config["mk_ingroup"] == config["vg_ref_backbone"] else "",
    shell:
        """
        mkdir -p results/mk/cds
        python3 scripts/mk_extract_cds.py \\
            --gff {input.gff3} \\
            --fasta {input.fasta} \\
            --species {wildcards.species} \\
            --fna {output.fna} \\
            --map {output.map} \\
            {params.chrom_prefix_arg}
        """


rule mk_codon_align:
    # MAFFT protein alignment + pal2nal codon alignment per orthogroup,
    # parallel across orthogroups (xargs -P {threads}). Failures are
    # recorded as <orthogroup>.failed inside the output directory.
    input:
        pairs="results/mk/ortholog_pairs.tsv",
        ingroup_fna="results/mk/cds/{ingroup}.fna".format(ingroup=config["mk_ingroup"]),
        outgroup_fna="results/mk/cds/{outgroup}.fna".format(outgroup=config["mk_outgroup"]),
    output:
        directory("results/mk/alignments"),
    conda: "../envs/paml.yaml"
    params:
        ingroup=config["mk_ingroup"],
        outgroup=config["mk_outgroup"],
        ingroup_faa=lambda wc: config["mk_species_annotations"][config["mk_ingroup"]]["faa"],
        outgroup_faa=lambda wc: config["mk_species_annotations"][config["mk_outgroup"]]["faa"],
    shell:
        """
        mkdir -p results/mk/alignments
        rm -f results/mk/alignments/*.fna results/mk/alignments/*.failed
        awk 'NR > 1 {{print $1}}' {input.pairs} | sort -u | \\
            xargs -P {threads} -I {{}} python3 scripts/mk_align_pairs.py \\
                --pairs {input.pairs} \\
                --only {{}} \\
                --ingroup {params.ingroup} \\
                --outgroup {params.outgroup} \\
                --ingroup-fna {input.ingroup_fna} \\
                --outgroup-fna {input.outgroup_fna} \\
                --ingroup-faa {params.ingroup_faa} \\
                --outgroup-faa {params.outgroup_faa} \\
                --outdir results/mk/alignments
        test -n "$(ls results/mk/alignments/*.fna 2>/dev/null)"
        """


rule mk_frequency_filter:
    # Stream each pool's genome-wide grenedalf frequency table down to just
    # the polymorphic rows mk_polymorphism classifies (constant memory + C
    # speed in awk). Output: chrom,pos,ref,alt,alt_count,depth,alt_freq,
    # no header. The thresholds are the same config keys the counting step
    # uses, so the intermediate cannot drift; temp() keeps the extracts out
    # of the durable tree.
    input:
        freq=lambda wc: (
            "results/vg/{variant}_vs_{backbone}/{ID}/"
            "grenedalf_results_frequency.csv".format(
                variant=config["vg_ref_variant"],
                backbone=config["vg_ref_backbone"],
                ID=wc.ID))
    output:
        temp("results/mk/frequency/{ID}_polymorphic.csv")
    params:
        min_depth=config["mk_min_depth"],
        min_alt_count=config["mk_min_alt_count"],
    shell:
        """
        mkdir -p results/mk/frequency
        awk -F, -v d={params.min_depth} -v a={params.min_alt_count} \\
            'NR == 1 || $1 == "CHROM" {{ next }}
             ($6 + 0) >= a && ($7 + 0) >= d && ($8 + 0) > 0 && ($8 + 0) < 1 \\
                {{ print $1 "," $2 "," $3 "," $4 "," $6 "," $7 "," $8 }}' \\
            {input.freq} > {output}
        """


rule mk_polymorphism:
    # Pn/Ps per gene per pool: each polymorphic site is classified by
    # translating the reference codon and the codon carrying the ALT
    # allele (synonymous -> Ps, nonsynonymous/stop-gain -> Pn). No
    # degeneracy annotation is needed and 2-fold sites count.
    input:
        map="results/mk/cds/{ingroup}_map.tsv".format(ingroup=config["mk_ingroup"]),
        cds="results/mk/cds/{ingroup}.fna".format(ingroup=config["mk_ingroup"]),
        freq=expand("results/mk/frequency/{ID}_polymorphic.csv", ID=sample_ids),
    output:
        pooled="results/mk/polymorphism/pooled_polymorphism.tsv",
        sites="results/mk/polymorphism/polymorphic_sites.tsv",
        per_pool=expand("results/mk/polymorphism/{ID}_polymorphism.tsv",
                        ID=sample_ids),
    conda: "../envs/mk.yaml"
    params:
        freq_args=" ".join(
            "--freq {ID}=results/mk/frequency/{ID}_polymorphic.csv".format(ID=ID)
            for ID in sample_ids),
        min_depth=config["mk_min_depth"],
        min_alt_count=config["mk_min_alt_count"],
    shell:
        """
        mkdir -p results/mk/polymorphism
        python3 scripts/mk_polymorphism.py \\
            --map {input.map} \\
            --cds {input.cds} \\
            --min-depth {params.min_depth} \\
            --min-alt-count {params.min_alt_count} \\
            {params.freq_args} \\
            --outdir results/mk/polymorphism
        """


rule mk_divergence:
    # Fixed synonymous/nonsynonymous differences per gene from the codon
    # alignments; codons polymorphic in any pool are excluded.
    input:
        alignments="results/mk/alignments",
        map="results/mk/cds/{ingroup}_map.tsv".format(ingroup=config["mk_ingroup"]),
        sites="results/mk/polymorphism/polymorphic_sites.tsv",
    output:
        divergence="results/mk/divergence.tsv",
        failures="results/mk/divergence_failures.tsv",
    conda: "../envs/mk.yaml"
    params:
        ingroup=config["mk_ingroup"],
        outgroup=config["mk_outgroup"],
    shell:
        """
        python3 scripts/mk_divergence.py \\
            --alignments {input.alignments} \\
            --map {input.map} \\
            --polymorphic-sites {input.sites} \\
            --ingroup {params.ingroup} \\
            --outgroup {params.outgroup} \\
            --out {output.divergence} \\
            --failures {output.failures}
        """


rule mk_codeml:
    # Pairwise codeml (runmode -2) dN/dS per 1:1 ortholog; parallel across
    # orthogroups. Alignments with stops/frameshifts are reported instead
    # of run; per-gene rows are concatenated into one table.
    input:
        alignments="results/mk/alignments",
    output:
        "results/mk/dnds.tsv"
    conda: "../envs/paml.yaml"
    params:
        ingroup=config["mk_ingroup"],
        outgroup=config["mk_outgroup"],
    shell:
        """
        mkdir -p results/mk/codeml
        rm -f results/mk/codeml/*.tsv
        ls results/mk/alignments/*.fna | \\
            xargs -P {threads} -I {{}} python3 scripts/mk_codeml.py \\
                --alignment {{}} \\
                --ingroup {params.ingroup} \\
                --outgroup {params.outgroup} \\
                --outdir results/mk/codeml
        printf 'orthogroup\\tgene\\tdN\\tdS\\tomega\\tt\\tstatus\\n' > {output}
        cat $(ls results/mk/codeml/*.tsv | sort) >> {output}
        """


rule mk_stats:
    # Standard MK statistics (alpha, NI, DoS, Fisher exact + BH) per gene,
    # per pool and pooled.
    input:
        divergence="results/mk/divergence.tsv",
        pooled="results/mk/polymorphism/pooled_polymorphism.tsv",
        per_pool=expand("results/mk/polymorphism/{ID}_polymorphism.tsv",
                        ID=sample_ids),
    output:
        pooled="results/mk/mk_stats_pooled.tsv",
        per_pool=expand("results/mk/mk_stats_{ID}.tsv", ID=sample_ids),
    conda: "../envs/mk.yaml"
    params:
        poly_args=" ".join(
            "--polymorphism {ID}=results/mk/polymorphism/{ID}_polymorphism.tsv".format(ID=ID)
            for ID in sample_ids),
    shell:
        """
        python3 scripts/mk_stats.py \\
            --divergence {input.divergence} \\
            --pooled {input.pooled} \\
            {params.poly_args} \\
            --outdir results/mk
        """
