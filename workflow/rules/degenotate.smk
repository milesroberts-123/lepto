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
