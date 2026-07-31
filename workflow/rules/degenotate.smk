rule degenotate_run:
    input:
        gff=lambda wildcards: config["featurecounts_annotations"][wildcards.species],
        fasta=lambda wildcards: config["bwa_refs"][wildcards.species]
    output:
        "results/degenotate/{species}/degeneracy-all-sites.bed"
    conda: "../envs/degenotate.yaml"
    shell:
        """
        mkdir -p results/degenotate/{wildcards.species}
        degenotate.py -a {input.gff} -g {input.fasta} -o results/degenotate/{wildcards.species}
        """

rule degenotate_subset:
    input:
        "results/degenotate/{species}/degeneracy-all-sites.bed"
    output:
        zero="results/degenotate/{species}/degeneracy-zero-sites.bed",
        four="results/degenotate/{species}/degeneracy-four-sites.bed",
    shell:
        """
        awk '(($5 == 0))' {input} > {output.zero}
        awk '(($5 == 4))' {input} > {output.four}
        """
