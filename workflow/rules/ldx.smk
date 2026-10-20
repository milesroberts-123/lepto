# LDx: linkage disequilibrium (r^2) from pooled alignments (Feder,
# Petrov & Bergland 2012, PLoS ONE 7(11):e48588, doi
# 10.1371/journal.pone.0048588). LDx rebuilds short-range 2-locus
# haplotype tables from reads intersecting pairs of polymorphic sites
# (same read or read pair), reporting a naive "direct" r^2 from the
# intersecting reads plus an approximate-MLE r^2 anchored to allele
# frequencies from all reads, with CIs at +/-1.92 log-likelihood units
# on a 0.01 grid. LDx parses ONLY column 2 (POS) of its SNP list and
# ignores CHROM, so a run must see a single-scaffold SAM and a
# single-scaffold SNP list -- genome-wide inputs would silently pair
# positions across scaffolds. Hence the checkpoint scaffold list feeds
# a per-(pool, scaffold) slice chain (samtools view slice + SNP-list
# slice -> LDx -> per-pool merged table). All options must precede the
# two positional inputs, and -h sets the maximum site read depth
# (repeat/paralog guard), not help; depth is kept only for
# l < depth < h, and -i is strictly greater, so 11 means 12+ shared
# reads. Site-depth thresholds scale per pool from the mean per-site
# depth of its grenedalf frequency table: -l = max(floor, mean/
# divisor), -h = multiple x mean. The binary is provisioned outside
# the repo and referenced by config["ldx_binary"] (see AGENTS.md); the
# LDx step itself needs only system Perl (no conda env). Consuming the
# vg BAMs and grenedalf frequency CSVs couples ldx targets to the vg
# chain whenever those products are missing (run vg_diversity_all
# first), same acceptance as freqk/kmate.

import glob
from snakemake.exceptions import WorkflowError


def ldx_scaffolds(wildcards):
    # Checkpoint re-expansion: the scaffold list (from the backbone
    # prefixed dictionary) only exists once the checkpoint has run, so
    # downstream consumers re-evaluate the DAG at that point.
    chk = checkpoints.ldx_scaffold_list.get(
        variant=wildcards.variant, backbone=wildcards.backbone)
    scaffolds = []
    with open(str(chk.output[0])) as fh:
        for line in fh:
            line = line.strip()
            if line:
                scaffolds.append(line)
    if not scaffolds:
        raise WorkflowError(
            "ldx_scaffold_list produced an empty scaffold list from "
            + str(chk.output[0]))
    return scaffolds


checkpoint ldx_scaffold_list:
    # One LDx run per scaffold; the backbone dict's @SQ lines are the
    # authoritative scaffold list.
    input:
        dict="results/vg/{backbone}_prefixed.dict"
    output:
        "results/ldx/{variant}_vs_{backbone}/scaffolds.txt"
    shell:
        """
        mkdir -p $(dirname {output})
        grep '^@SQ' {input.dict} | cut -f2 | sed 's/^SN://' > {output}
        """


rule ldx_depth_mean:
    # Mean per-site depth for this pool from the grenedalf frequency
    # table (column 7 depth, header skipped, data-bearing sites only).
    # Scales the LDx min/max site-depth thresholds.
    input:
        freq="results/vg/{variant}_vs_{backbone}/{ID}/grenedalf_results_frequency.csv"
    output:
        "results/ldx/{variant}_vs_{backbone}/{ID}/depth_mean.txt"
    benchmark:
        "benchmarks/ldx/ldx_depth_mean_{variant}_vs_{backbone}_{ID}.bench"
    shell:
        """
        mkdir -p $(dirname {output})
        awk -F, 'NR > 1 && $7 > 0 {{ n++; s += $7 }} END {{ if (n > 0) printf "%.1f\\n", s/n }}' \\
            {input.freq} > {output}
        """


rule ldx_snplist:
    # Sites polymorphic within THIS pool (0 < alt freq < 1 in column 8,
    # with a depth floor on column 7): chrom + 1-based position, tab
    # separated, sorted; LDx reads the position column. LDx re-checks
    # depth and MAF against the reads, so the floor only drops dead
    # sites from the list.
    input:
        freq="results/vg/{variant}_vs_{backbone}/{ID}/grenedalf_results_frequency.csv"
    output:
        "results/ldx/{variant}_vs_{backbone}/{ID}/snplist.tsv"
    params:
        min_depth=config["ldx_snplist_min_depth"]
    benchmark:
        "benchmarks/ldx/ldx_snplist_{variant}_vs_{backbone}_{ID}.bench"
    shell:
        """
        mkdir -p $(dirname {output})
        awk -F, -v d={params.min_depth} 'NR > 1 && $7 >= d && $8 > 0 && $8 < 1 {{ print $1"\\t"$2 }}' \\
            {input.freq} | sort -k1,1 -k2,2n > {output}
        """


rule ldx_slice:
    # Per-(pool, scaffold) LDx inputs: primary paired alignments on
    # this scaffold (-F 0x904 drops unmapped/secondary/supplementary
    # reads because LDx ignores SAM flags) and the same scaffold's SNP
    # rows. The SNP list must stay two-column (chrom+pos): LDx parses
    # column 2 only, so a single-column position list is invisible to
    # it.
    input:
        bam="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam",
        bai="results/vg/{variant}_vs_{backbone}/{ID}/{ID}.bam.bai",
        snplist="results/ldx/{variant}_vs_{backbone}/{ID}/snplist.tsv"
    output:
        sam=temp("results/ldx/{variant}_vs_{backbone}/{ID}/{ID}.{scaffold}.sam"),
        snps=temp("results/ldx/{variant}_vs_{backbone}/{ID}/{ID}.{scaffold}.snps")
    conda: "../envs/bcftools.yaml"
    params:
        min_mapq=config["ldx_min_mapq"]
    benchmark:
        "benchmarks/ldx/ldx_slice_{variant}_vs_{backbone}_{ID}_{scaffold}.bench"
    shell:
        """
        samtools view -F 0x904 -q {params.min_mapq} {input.bam} {wildcards.scaffold} > {output.sam}
        awk -v s={wildcards.scaffold} '$1 == s' {input.snplist} > {output.snps}
        """


rule ldx_run:
    # Adaptive site-depth window from the pool's mean per-site depth:
    # -l = max(floor, mean/divisor), -h = multiple x mean (empty mean
    # file falls back to LDx defaults via depth=0 -> -h 0, which keeps
    # no sites). Options precede the positional inputs; a per-scaffold
    # tag is prepended because LDx output carries no scaffold column.
    input:
        sam="results/ldx/{variant}_vs_{backbone}/{ID}/{ID}.{scaffold}.sam",
        snps="results/ldx/{variant}_vs_{backbone}/{ID}/{ID}.{scaffold}.snps",
        depth="results/ldx/{variant}_vs_{backbone}/{ID}/depth_mean.txt"
    output:
        "results/ldx/{variant}_vs_{backbone}/{ID}/{ID}.{scaffold}.pairs.tsv"
    params:
        binary=config["ldx_binary"],
        min_depth_floor=config["ldx_min_depth_floor"],
        min_depth_divisor=config["ldx_min_depth_divisor"],
        max_depth_multiple=config["ldx_max_depth_multiple"],
        searchlen=config["ldx_searchlen"],
        baseq=config["ldx_baseq"],
        maf=config["ldx_maf"],
        min_intersect=config["ldx_min_intersect"]
    threads: 1
    benchmark:
        "benchmarks/ldx/ldx_run_{variant}_vs_{backbone}_{ID}_{scaffold}.bench"
    shell:
        """
        depth=$(cat {input.depth})
        depth=${{depth:-0}}
        low=$(awk -v d=$depth -v f={params.min_depth_floor} -v n={params.min_depth_divisor} \\
              'BEGIN {{ l = d/n; if (l < f) l = f; printf "%d", int(l) }}')
        high=$(awk -v d=$depth -v m={params.max_depth_multiple} \\
              'BEGIN {{ printf "%d", int(d*m) }}')
        perl {params.binary} \\
            -l $low -h $high -s {params.searchlen} -q {params.baseq} -a {params.maf} -i {params.min_intersect} \\
            {input.sam} {input.snps} \\
            | awk -v sc={wildcards.scaffold} '{{ print sc"\\t"$0 }}' > {output}
        """


def ldx_pair_inputs(wildcards):
    # All per-scaffold pair tables for this pool, re-expanded through
    # the checkpoint's scaffold list.
    scaffolds = ldx_scaffolds(wildcards)
    return expand(
        "results/ldx/{variant}_vs_{backbone}/{ID}/{ID}.{scaffold}.pairs.tsv",
        variant=wildcards.variant,
        backbone=wildcards.backbone,
        ID=wildcards.ID,
        scaffold=scaffolds)


rule ldx_merge:
    # Per-pool concatenated pair table: rows are already scaffold-tagged,
    # so a plain concat is sufficient (LDx emits no header).
    input:
        pairs=ldx_pair_inputs
    output:
        "results/ldx/{variant}_vs_{backbone}/{ID}.ldx_pairs.tsv"
    benchmark:
        "benchmarks/ldx/ldx_merge_{variant}_vs_{backbone}_{ID}.bench"
    shell:
        "cat {input} > {output}"