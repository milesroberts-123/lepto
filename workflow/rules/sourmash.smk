# sourmash k-mer sketches and repeat-robust mutation rate estimation.
# Chain: fasta -> sketch (per fasta) -> mutation_rate (s vs t pair).

rule sourmash_sketch:
    input:
        lambda wildcards: config["sourmash_fastas"][wildcards.fasta]
    output:
        "results/sourmash/{fasta}.sig"
    conda: "../envs/sourmash.yaml"
    params:
        k=config["sourmash_k"],
        scaled=config["sourmash_scaled"],
        sketch_mode=lambda wildcards: config["sourmash_sketch_modes"][wildcards.fasta]
    shell:
        """
        mkdir -p results/sourmash
        sourmash scripts sketch {input} \
            --sketch-mode {params.sketch_mode} \
            -o {output} \
            -k {params.k} \
            --scaled {params.scaled}
        """


rule sourmash_mutation_rate:
    input:
        s="results/sourmash/{s}.sig",
        t="results/sourmash/{t}.sig"
    output:
        "results/sourmash/mutation_rate_{s}_vs_{t}.txt"
    conda: "../envs/sourmash.yaml"
    params:
        estimator=config["sourmash_estimator"]
    shell:
        """
        sourmash scripts mutation_rate \
            --estimator {params.estimator} \
            --s-sig {input.s} \
            --t-sig {input.t} > {output}
        """


rule sourmash_all:
    input:
        "results/sourmash/mutation_rate_{s}_vs_{t}.txt".format(
            s=config["sourmash_s"], t=config["sourmash_t"]),
        expand("results/sourmash/{fasta}.sig", fasta=config["sourmash_fastas"].keys())
