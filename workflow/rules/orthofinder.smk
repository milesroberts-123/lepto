# OrthoFinder pipeline: filter transcripts, prefix sequences by species,
# stage proteomes, and run OrthoFinder to produce orthogroups.

rule orthofinder_filter_t1:
    input:
        "results/featurecounts/{species}_expressed.faa"
    output:
        "results/orthofinder/filtered/{species}.t1.faa"
    conda: "../envs/orthofinder.yaml"
    shell:
        """
        mkdir -p results/orthofinder/filtered
        seqkit grep -n -r -p "\\.t1$" {input} > {output}
        """

rule orthofinder_unique_prefix:
    input:
        "results/orthofinder/filtered/{species}.t1.faa"
    output:
        "results/orthofinder/prefixed/{species}.faa"
    conda: "../envs/orthofinder.yaml"
    shell:
        """
        seqkit replace -p '^' -r '{wildcards.species}.' {input} -o {output}
        """

rule orthofinder_stage:
    # One staging rule for both proteome sources: in-repo species
    # (hap1/hap2 -> results/orthofinder/prefixed/<species>.faa) and
    # external proteomes from config. Previously two rules produced the
    # same output path, and the external one raised KeyError for hap1/hap2
    # whenever snakemake evaluated it (fatal for consumers under recent
    # snakemake versions).
    input:
        lambda wc: config["orthofinder_external_proteomes"].get(wc.species)
        or f"results/orthofinder/prefixed/{wc.species}.faa"
    output:
        "results/orthofinder/proteomes/{species}.faa"
    localrule: True
    shell:
        """
        mkdir -p results/orthofinder/proteomes
        cp {input} {output}
        """


rule orthofinder_run:
    input:
        expand("results/orthofinder/proteomes/{species}.faa",
               species=list(config["panther_input_fastas"].keys()) + list(config.get("orthofinder_external_proteomes", {}).keys()))
    output:
        directory("results/orthofinder/proteomes/OrthoFinder")
    conda: "../envs/orthofinder.yaml"
    shell:
        """
        orthofinder -f results/orthofinder/proteomes \\
            -t {threads} \\
            -a {threads} \\
            -S diamond \\
            -og 
        """
