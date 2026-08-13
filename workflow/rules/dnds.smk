rule mafft_align:
    input:
        config["dnds_aa_fasta"]
    output:
        "results/dnds/aa.aligned"
    conda: "../envs/paml.yaml"
    shell:
        """
        mkdir -p results/dnds
        linsi {input} > {output}
        """


rule pal2nal:
    input:
        aa="results/dnds/aa.aligned",
        cds=config["dnds_cds_fasta"]
    output:
        "results/dnds/cds.paml"
    conda: "../envs/paml.yaml"
    shell:
        """
        pal2nal.pl {input.aa} {input.cds} -nogap -output paml > {output}
        """


rule yn00:
    input:
        cds="results/dnds/cds.paml",
        ctl=config["dnds_yn00_ctl"]
    output:
        "results/dnds/yn_results.txt"
    conda: "../envs/paml.yaml"
    shell:
        """
        cp {input.ctl} results/dnds/yn00.ctl
        cd results/dnds && yn00 yn00.ctl
        """
