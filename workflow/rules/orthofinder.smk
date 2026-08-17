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
    input:
        "results/orthofinder/prefixed/{species}.faa"
    output:
        "results/orthofinder/proteomes/{species}.faa"
    shell:
        """
        mkdir -p results/orthofinder/proteomes
        ln -sf $(realpath {input}) {output}
        """


rule orthofinder_stage_external:
    input:
        lambda wildcards: config["orthofinder_external_proteomes"][wildcards.species]
    output:
        "results/orthofinder/proteomes/{species}.faa"
    shell:
        """
        mkdir -p results/orthofinder/proteomes
        ln -sf $(realpath {input}) {output}
        """


rule orthofinder_run:
    input:
        expand("results/orthofinder/proteomes/{species}.faa",
               species=list(config["panther_input_fastas"].keys()) + list(config.get("orthofinder_external_proteomes", {}).keys()))
    output:
        "results/orthofinder/Orthogroups.tsv"
    conda: "../envs/orthofinder.yaml"
    shell:
        """
        orthofinder -f results/orthofinder/proteomes \\
            -t {threads} \\
            -a {threads} \\
            -S diamond \\
            -og 
        ln -sf $(ls -d results/orthofinder/proteomes/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv) {output}
        """
