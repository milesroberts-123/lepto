# BRAKER3 genome annotation from RNA-seq.
# Chain: SRA download -> fastp trimming -> HISAT2 alignment -> BRAKER3
# (singularity container) gene prediction per reference genome.
# Requires --use-singularity; sra_download runs locally (login node has
# internet, compute nodes do not).


rule sra_download:
    output:
        r1="results/rna/{acc}_1.fastq",
        r2="results/rna/{acc}_2.fastq"
    conda: "../envs/sratools.yaml"
    shell:
        """
        mkdir -p results/rna
        prefetch -O results/rna {wildcards.acc}
        fasterq-dump -e {threads} -O results/rna results/rna/{wildcards.acc}/{wildcards.acc}.sra
        rm -rf results/rna/{wildcards.acc}
        """


rule fastp_rna:
    input:
        r1="results/rna/{acc}_1.fastq",
        r2="results/rna/{acc}_2.fastq"
    output:
        r1="results/rna/fastp/{acc}_1.fastq",
        r2="results/rna/fastp/{acc}_2.fastq",
        json="results/rna/fastp/{acc}.json"
    conda: "../envs/fastp.yaml"
    params:
        unqual_limit=config["fastp_unqual_limit"],
        min_len=config["fastp_min_len"],
        qual_thresh=config["fastp_qual_thresh"],
        window_length=config["fastp_window_length"],
        n_base_limit=config["fastp_n_base_limit"]
    shell:
        """
        mkdir -p results/rna/fastp
        fastp --thread {threads} \
            --n_base_limit {params.n_base_limit} \
            -u {params.unqual_limit} \
            -q {params.qual_thresh} \
            -l {params.min_len} \
            --cut_tail \
            --cut_tail_window_size {params.window_length} \
            --cut_tail_mean_quality {params.qual_thresh} \
            --json {output.json} \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2}
        """


rule hisat2_index:
    input:
        genome=lambda wildcards: braker_genomes[wildcards.ref]
    output:
        expand("results/hisat2/{{ref}}/{{ref}}.{n}.ht2", n=range(1, 9))
    conda: "../envs/hisat2.yaml"
    shell:
        """
        mkdir -p results/hisat2/{wildcards.ref}
        hisat2-build -p {threads} {input.genome} results/hisat2/{wildcards.ref}/{wildcards.ref}
        """


rule hisat2_align:
    input:
        r1="results/rna/fastp/{acc}_1.fastq",
        r2="results/rna/fastp/{acc}_2.fastq",
        idx=rules.hisat2_index.output
    output:
        "results/hisat2/{ref}/{acc}.bam"
    conda: "../envs/hisat2.yaml"
    shell:
        """
        hisat2 --dta -p {threads} -x results/hisat2/{wildcards.ref}/{wildcards.ref} \
            -1 {input.r1} -2 {input.r2} | \
            samtools sort -@ {threads} -o {output}
        """


def braker_bams(wildcards):
    return expand(
        "results/hisat2/{ref}/{acc}.bam",
        ref=wildcards.ref,
        acc=config["rna_sra_accessions"],
    )


rule braker_run:
    input:
        genome=lambda wildcards: braker_genomes[wildcards.ref],
        bams=braker_bams
    output:
        gtf="results/braker/{ref}/braker.gtf",
        gff3="results/braker/{ref}/braker.gff3",
        aa="results/braker/{ref}/braker.aa",
        codingseq="results/braker/{ref}/braker.codingseq"
    container: config["braker_container"]
    params:
        species=config["braker_species"],
        cfg=lambda wildcards: f"results/braker/{wildcards.ref}/augustus_config",
        workdir=lambda wildcards: f"results/braker/{wildcards.ref}/braker"
    shell:
        """
        mkdir -p results/braker/{wildcards.ref}
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        if [ ! -d {params.cfg} ]; then
            cp -r $AUGUSTUS_CONFIG_PATH {params.cfg}
        fi
        bams=$(echo {input.bams} | tr ' ' ',')
        braker.pl --species={params.species} \
            --genome={input.genome} \
            --bam=$bams \
            --threads={threads} \
            --gff3 \
            --softmasking \
            --AUGUSTUS_CONFIG_PATH=$(pwd)/{params.cfg} \
            --workingdir=$(pwd)/{params.workdir}
        cp {params.workdir}/braker.gtf {output.gtf}
        cp {params.workdir}/braker.gff3 {output.gff3}
        cp {params.workdir}/braker.aa {output.aa}
        cp {params.workdir}/braker.codingseq {output.codingseq}
        """
