rule degenotate_run:
    input:
        gff="results/featurecounts/{species}_expressed.gff3",
        fasta=lambda wildcards: config["bwa_refs"][wildcards.species]
    output:
        "results/degenotate/{species}/degeneracy-all-sites.bed"
    conda: "../envs/degenotate.yaml"
    shell:
        """
        mkdir -p results/degenotate/{wildcards.species}
        degenotate.py -a {input.gff} -g {input.fasta} --overwrite -o results/degenotate/{wildcards.species}
        """

rule degenotate_subset:
    input:
        "results/degenotate/{species}/degeneracy-all-sites.bed"
    output:
        all="results/degenotate/{species}/degeneracy-cds-sites.bed",
        zero="results/degenotate/{species}/degeneracy-zero-sites.bed",
        four="results/degenotate/{species}/degeneracy-four-sites.bed",
    shell:
        """
        cut -f 1-3 {input} > {output.all}
        awk '(($5 == 0))' {input} | cut -f 1-3 > {output.zero}
        awk '(($5 == 4))' {input} | cut -f 1-3 > {output.four}
        """
