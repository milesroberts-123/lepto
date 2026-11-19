# BRAKER3 genome annotation from RNA-seq and OrthoDB protein evidence.
# Chain: OrthoDB/SRA download -> fastp trimming -> HISAT2 alignment ->
# BRAKER3 (singularity container) gene prediction per reference genome.
# Requires --use-singularity; odb_download and sra_download need internet
# (login node has internet, compute nodes do not).


rule minimap2_braker_paf:
    input:
        backbone=config["minimap2_braker"]["backbone"],
        variant=config["minimap2_braker"]["variant"]
    output:
        "results/braker/{variant}_vs_{backbone}.paf"
    conda: "../envs/minimap2.yaml"
    shell:
        "minimap2 -t {threads} -x asm20 -f 0.02 {input.backbone} {input.variant} > {output}"


rule odb_download:
    # OrthoDB v12 Viridiplantae proteins used as BRAKER3 protein evidence
    # (switches GeneMark to ETP mode). Needs internet: run on a login node.
    output:
        fa="results/braker/odb12/Viridiplantae.fa"
    conda: "../envs/odb_download.yaml"
    params:
        url=config["braker_odb_url"],
        md5=config["braker_odb_md5"]
    shell:
        """
        mkdir -p $(dirname {output})
        wget -c -O {output}.gz {params.url}
        echo "{params.md5}  {output}.gz" | md5sum -c -
        gunzip {output}.gz
        """


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


def braker_protein(wildcards):
    # Protein evidence only when enabled; an empty list keeps the OrthoDB
    # download out of the DAG when braker_use_protein is false.
    if config["braker_use_protein"]:
        return ["results/braker/odb12/Viridiplantae.fa"]
    return []


def braker_protein_flag(wildcards):
    if config["braker_use_protein"]:
        return "--prot_seq=results/braker/odb12/Viridiplantae.fa"
    return ""


rule braker_run:
    input:
        genome=lambda wildcards: braker_genomes[wildcards.ref],
        bams=braker_bams,
        protein=braker_protein
    output:
        gtf="results/braker/{ref}/braker.gtf",
        gff3="results/braker/{ref}/braker.gff3",
        aa="results/braker/{ref}/braker.aa",
        codingseq="results/braker/{ref}/braker.codingseq"
    container: config["braker_container"]
    params:
        species=config["braker_species"],
        cfg=lambda wildcards: f"results/braker/{wildcards.ref}/augustus_config",
        workdir=lambda wildcards: f"results/braker/{wildcards.ref}/braker",
        protein_flag=braker_protein_flag
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
            {params.protein_flag} \
            --gff3 \
            --softmasking \
            --AUGUSTUS_CONFIG_PATH=$(pwd)/{params.cfg} \
            --workingdir=$(pwd)/{params.workdir}
        cp {params.workdir}/braker.gtf {output.gtf}
        cp {params.workdir}/braker.gff3 {output.gff3}
        cp {params.workdir}/braker.aa {output.aa}
        cp {params.workdir}/braker.codingseq {output.codingseq}
        """


# Expression filtering of BRAKER gene models: count reads per transcript
# over the braker.gff3 annotations with the same RNA-seq BAMs BRAKER used,
# then drop gene models with zero counts and their offspring from the GFF3
# and the protein FASTA. Adapts the featurecounts subworkflow above;
# transcript IDs are also emitted at gene level (stripped .tN suffix) so
# the top-level gene features survive the GFF3 filter.

rule braker_featurecounts_run:
    input:
        gff="results/braker/{ref}/braker.gff3",
        bams=lambda wildcards: expand(
            "results/hisat2/{ref}/{acc}.bam",
            ref=wildcards.ref,
            acc=config["rna_sra_accessions"],
        )
    output:
        "results/braker/{ref}/braker.counts"
    conda: "../envs/featurecounts.yaml"
    params:
        gff_feature=config["featurecounts_gff_feature"],
        paired="-p" if config["featurecounts_paired"] else ""
    shell:
        """
        mkdir -p $(dirname {output})
        featureCounts -a {input.gff} \
            -o {output} \
            -T {threads} \
            {params.paired} \
            -g {params.gff_feature} \
            {input.bams}
        """


rule braker_featurecounts_filter_expressed:
    input:
        "results/braker/{ref}/braker.counts"
    output:
        "results/braker/{ref}/braker_expressed_genes.txt"
    conda: "../envs/featurecounts.yaml"
    params:
        # last column is length; counts are the remaining columns
        n_bams=lambda wildcards: len(config["rna_sra_accessions"])
    shell:
        """
        awk -F'\t' 'NR > 2 && /^g/ {{
            sum = 0
            for (i = NF - {params.n_bams} + 1; i <= NF; i++) sum += $i
            if (sum > 0) {{
                print $1
                sub(/\\.t[0-9]+$/, "", $1)
                print $1
            }}
        }}' {input} | sort -u > {output}
        """


rule braker_featurecounts_filter_gff:
    input:
        glist="results/braker/{ref}/braker_expressed_genes.txt",
        gff="results/braker/{ref}/braker.gff3"
    output:
        "results/braker/{ref}/braker_expressed.gff3"
    conda: "../envs/featurecounts.yaml"
    shell:
        """
        awk -F'\\t' 'NR==FNR{{genes[$1]=1; next}}
            /^#/{{print; next}}
            {{
                match($9, /ID=([^;]+)/, id); match($9, /Parent=([^;]+)/, parent)
                if (id[1] in genes || parent[1] in genes) print
            }}' {input.glist} {input.gff} > {output}
        """


rule braker_filter_expressed_faa:
    input:
        genes="results/braker/{ref}/braker_expressed_genes.txt",
        faa="results/braker/{ref}/braker.aa"
    output:
        "results/braker/{ref}/braker_expressed.faa"
    conda: "../envs/featurecounts.yaml"
    shell:
        """
        # anchor patterns to full header match so g1 does not match g10.t1
        sed 's/.*/^&$/' {input.genes} > {input.genes}.patterns
        seqkit grep -n -r -f {input.genes}.patterns {input.faa} > {output}
        rm {input.genes}.patterns
        """
