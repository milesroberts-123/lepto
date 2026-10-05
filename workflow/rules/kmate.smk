# kMate: haplotype-mixture allele-frequency estimates for the P-group
# pools (see config kmate_groups). Alignment-free: pooled reads are
# counted as k-mers and an EM fits the hap1/hap2 founder mixture per
# unit (config kmate_unit, default chromosome), which is projected to
# per-variant allele frequencies through the panel's founder x variant
# matrix. The panel comes from freqk's norm.vcf.gz — the same filtered
# VCF freqk indexes but WITHOUT the '-v snps' cut, so indels stay in —
# with an all-reference hap1 backbone column appended and hap2
# haploidized. kmate-only targets therefore schedule the whole vg+freqk
# chain whenever the norm.vcf.gz temp is missing (same accepted coupling
# as freqk; run vg_diversity_all first). The 2-founder panel requires
# kmate_min_ac 1: the default --min-ac 2 empties the k-mer matrix (the
# production filter also drops ac==1 private k-mers, so with F=2 the
# keep window is exactly ac==1; verified e2e). kMate enforces one
# --chroms value per 'kmate run' invocation (--var-pa-prefix is
# per-chromosome), so runs are one job PER CHROM and a final concat
# rule assembles the per-pool TSV. Pool k-mers are counted once per
# sample into a temporary jellyfish DB (kmate_kmer_db), which is
# deleted after the concat; its --hash-size must be generous (pool 1 is
# 456 GB of fastq, config kmate_hash_size 32G). The wiki warns
# shared-filesystem DB reads are slow at scale — acceptable at 4 pools.

import os
import glob
from snakemake.exceptions import WorkflowError


def kmate_exclusion_input(wildcards):
    # Same W-pool exclusion as freqk_filter_vcf: mirrors
    # freqk_exclusion_input so the exclusion chain follows the
    # freqk_use_exclusion_filter toggle (no cost while it is off).
    if config["freqk_use_exclusion_filter"]:
        return "results/freqk/{variant}_vs_{backbone}/polymorphic_exclusion.bed"
    return []


def kmate_exclusion_stage(wildcards, input):
    # Keyed off the resolved input, not config; built by concatenation
    # so no format placeholders survive into the shell.
    if input.excl:
        return "| bcftools view -T ^{}".format(input.excl)
    return ""


rule kmate_panel_vcf:
    # Panel VCF: freqk's filtered VCF kept biallelic and QUAL-filtered
    # but WITHOUT the -v snps cut (kMate also gets indels), haploidized
    # and given the backbone as an all-reference founder column. awk
    # rewrites every record to two haploid GT columns: hap1=0
    # (all-reference backbone founder) and hap2=1 (variant founder,
    # carried at every record by construction). A minimal VCF header is
    # emitted because the reheadered source header is not preserved;
    # bgzip+tabix are required by the kmate builders.
    input:
        vcf="results/vg/{variant}_vs_{backbone}.vcf.gz",
        tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi",
        excl=kmate_exclusion_input
    output:
        panel="results/kmate/{variant}_vs_{backbone}/panel.vcf.gz",
        tbi="results/kmate/{variant}_vs_{backbone}/panel.vcf.gz.tbi"
    conda: "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/kmate/kmate_panel_vcf_{variant}_vs_{backbone}.bench"
    params:
        min_qual=config["freqk_min_qual"],
        excl_stage=kmate_exclusion_stage
    shell:
        """
        mkdir -p $(dirname {output.panel})
        bcftools norm -m -any {input.vcf} {params.excl_stage} \\
            | bcftools view -i 'QUAL>={params.min_qual}' -m2 -M2 \\
            | awk -F '\\t' -v OFS='\\t' \\
                'BEGIN{{print "##fileformat=VCFv4.2";
                        print "##FILTER=<ID=PASS,Description=\\"All filters passed\\">";
                        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"Genotype\\">";
                        print "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\thap1\\thap2"}} \\
                 !/^#/ {{print $1,$2,".",$4,$5,$6,".",".","GT",0,1}}' \\
            | bgzip -c > {output.panel}
        tabix -p vcf {output.panel}
        """


def kmate_chroms_from_checkpoint(checkpoint_output):
    # Checkpoint re-expansion: the chrom list is only known after
    # build_index has run, so downstream consumers re-evaluate against
    # whatever kmers.tsv.gz files the checkpoint produced.
    base = os.path.dirname(str(checkpoint_output))
    kmers = sorted(glob.glob(os.path.join(base, "index", "*_kmers.tsv.gz")))
    chroms = [os.path.basename(p)[:-len("_kmers.tsv.gz")] for p in kmers]
    if not chroms:
        raise WorkflowError(
            "kmate_build_index produced no *_kmers.tsv.gz under "
            + os.path.join(base, "index"))
    return chroms


checkpoint kmate_build_index:
    input:
        vcf="results/kmate/{variant}_vs_{backbone}/panel.vcf.gz",
        tbi="results/kmate/{variant}_vs_{backbone}/panel.vcf.gz.tbi",
        fasta="results/vg/{backbone}_prefixed.fasta",
        fai="results/vg/{backbone}_prefixed.fasta.fai"
    output:
        directory("results/kmate/{variant}_vs_{backbone}/index")
    conda: "../envs/kmate.yaml"
    benchmark:
        "benchmarks/kmate/kmate_build_index_{variant}_vs_{backbone}.bench"
    params:
        k=config["kmate_k"]
    shell:
        """
        kmate build-index \\
            --vcf {input.vcf} \\
            --ref {input.fasta} \\
            --out $(dirname {output})/index/panel \\
            -k {params.k} \\
            --haploid
        """


def kmate_kmer_pa_chroms(wildcards):
    chk = checkpoints.kmate_build_index.get(
        variant=wildcards.variant, backbone=wildcards.backbone)
    return kmate_chroms_from_checkpoint(chk.output[0])


rule kmate_build_kmer_pa:
    input:
        kmers="results/kmate/{variant}_vs_{backbone}/index/panel_{chrom}_kmers.tsv.gz",
        vcf="results/kmate/{variant}_vs_{backbone}/panel.vcf.gz",
        tbi="results/kmate/{variant}_vs_{backbone}/panel.vcf.gz.tbi",
        fasta="results/vg/{backbone}_prefixed.fasta",
        fai="results/vg/{backbone}_prefixed.fasta.fai"
    output:
        matrix="results/kmate/{variant}_vs_{backbone}/kmer_pa/kmer_pa_{chrom}.kmer_pa.npz",
        meta="results/kmate/{variant}_vs_{backbone}/kmer_pa/kmer_pa_{chrom}.meta.npz"
    conda: "../envs/kmate.yaml"
    benchmark:
        "benchmarks/kmate/kmate_build_kmer_pa_{variant}_vs_{backbone}_{chrom}.bench"
    params:
        min_ac=config["kmate_min_ac"],
        invariant_margin=config["kmate_invariant_margin"]
    shell:
        """
        kmate build-kmer-pa \\
            --kmers {input.kmers} \\
            --vcf {input.vcf} \\
            --ref {input.fasta} \\
            --chrom {wildcards.chrom} \\
            --out $(dirname {output.matrix})/kmer_pa_{wildcards.chrom} \\
            --treat-missing-as-n \\
            --filter-production \\
            --min-ac {params.min_ac} \\
            --invariant-margin {params.invariant_margin}
        """


rule kmate_build_var_pa:
    input:
        vcf="results/kmate/{variant}_vs_{backbone}/panel.vcf.gz"
    output:
        var_pa="results/kmate/{variant}_vs_{backbone}/var_pa/var_pa_{chrom}.var_pa.npz",
        called="results/kmate/{variant}_vs_{backbone}/var_pa/var_pa_{chrom}.var_called.npz",
        meta="results/kmate/{variant}_vs_{backbone}/var_pa/var_pa_{chrom}.meta.npz"
    conda: "../envs/kmate.yaml"
    benchmark:
        "benchmarks/kmate/kmate_build_var_pa_{variant}_vs_{backbone}_{chrom}.bench"
    shell:
        """
        kmate build-var-pa \\
            --vcf {input.vcf} \\
            --chrom {wildcards.chrom} \\
            --out $(dirname {output.var_pa})/var_pa_{wildcards.chrom}
        """


rule kmate_kmer_db:
    input:
        "results/fastp/{ID}.fastq"
    output:
        temp("results/kmate/kmerdb/{ID}.jf")
    conda: "../envs/kmate.yaml"
    benchmark:
        "benchmarks/kmate/kmate_kmer_db_{ID}.bench"
    params:
        hash_size=config["kmate_hash_size"]
    shell:
        """
        mkdir -p $(dirname {output})
        kmate build-kmer-db \\
            --reads {input} \\
            --out {output} \\
            --threads {threads} \\
            --hash-size {params.hash_size}
        """


rule kmate_run:
    # One invocation per chromosome (kMate requires exactly one --chroms
    # value with --var-pa-prefix). Reads come via the prebuilt jellyfish
    # DB (count-once path, identical AFs to re-counting, verified e2e).
    input:
        kmer_pa="results/kmate/{variant}_vs_{backbone}/kmer_pa/kmer_pa_{chrom}.kmer_pa.npz",
        kmer_meta="results/kmate/{variant}_vs_{backbone}/kmer_pa/kmer_pa_{chrom}.meta.npz",
        var_pa="results/kmate/{variant}_vs_{backbone}/var_pa/var_pa_{chrom}.var_pa.npz",
        var_called="results/kmate/{variant}_vs_{backbone}/var_pa/var_pa_{chrom}.var_called.npz",
        var_meta="results/kmate/{variant}_vs_{backbone}/var_pa/var_pa_{chrom}.meta.npz",
        reads="results/fastp/{ID}.fastq",
        jf="results/kmate/kmerdb/{ID}.jf"
    output:
        tsv=temp("results/kmate/{variant}_vs_{backbone}/{ID}_{chrom}_af.tsv"),
        hnpz="results/kmate/{variant}_vs_{backbone}/{ID}_{chrom}_af.h_per_chrom.npz"
    conda: "../envs/kmate.yaml"
    benchmark:
        "benchmarks/kmate/kmate_run_{variant}_vs_{backbone}_{ID}_{chrom}.bench"
    params:
        kmer_prefix="results/kmate/{variant}_vs_{backbone}/kmer_pa/kmer_pa",
        var_prefix="results/kmate/{variant}_vs_{backbone}/var_pa/var_pa",
        unit=config["kmate_unit"],
        out_prefix="results/kmate/{variant}_vs_{backbone}/{ID}_{chrom}_af"
    threads: 8
    shell:
        """
        kmate run \\
            --kmer-pa-prefix {params.kmer_prefix} \\
            --var-pa-prefix {params.var_prefix} \\
            --kmer-db {input.jf} \\
            --reads {input.reads} \\
            --sample {wildcards.ID} \\
            --out {params.out_prefix} \\
            --chroms {wildcards.chrom} \\
            --unit {params.unit}
        """


def kmate_concat_inputs(wildcards):
    # All per-chrom TSVs for this pool, re-expanded through the
    # checkpoint's chrom list; the jellyfish DB temp rides along so it
    # is deleted once every chrom run for the pool is done.
    chk = checkpoints.kmate_build_index.get(
        variant=wildcards.variant, backbone=wildcards.backbone)
    chroms = kmate_chroms_from_checkpoint(chk.output[0])
    return expand(
        "results/kmate/{variant}_vs_{backbone}/{ID}_{chrom}_af.tsv",
        variant=wildcards.variant,
        backbone=wildcards.backbone,
        ID=wildcards.ID,
        chrom=chroms)


rule kmate_concat:
    # Assemble the per-chrom frequency tables into the per-pool TSV
    # consumed by kmate_all; consumes the jellyfish DB temp so it is
    # deleted once all chrom runs for the pool are done.
    input:
        tsvs=kmate_concat_inputs,
        jf="results/kmate/kmerdb/{ID}.jf"
    output:
        "results/kmate/{variant}_vs_{backbone}/{ID}_af.tsv"
    benchmark:
        "benchmarks/kmate/kmate_concat_{variant}_vs_{backbone}_{ID}.bench"
    shell:
        """
        head -1 {input.tsvs[0]} > {output}
        for t in {input.tsvs}; do
            [ "$t" = "{input.tsvs[0]}" ] && continue
            tail -n +2 "$t" >> {output}
        done
        """


rule kmate_run_all:
    input:
        expand("results/kmate/{variant}_vs_{backbone}/{ID}_af.tsv",
               variant=[config["vg_ref_variant"]],
               backbone=[config["vg_ref_backbone"]],
               ID=p_sample_ids)
