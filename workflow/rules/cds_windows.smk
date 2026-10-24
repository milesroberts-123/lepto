# Fixed, non-overlapping windows across an assembly with the number of
# CDS-annotated bases in each window. Chain: faidx (index redirected to
# results/, never written next to the fasta) -> GFF3 CDS extraction to
# 0-based BED with merger of isoform-overlapping CDS (each base counted
# once) -> bedtools makewindows | bedtools coverage -> chrom start end
# cds_bases (all windows kept, zeros included). The window width is a
# wildcard {wlen} embedded in the output filename
# cds_windows_w{wlen}.bed (grenedalf interval naming pattern), so
# results for different window sizes coexist; cds_windows_all expands
# over config["cds_windows_sizes"].


rule cds_windows_faidx:
    input:
        fasta=config["cds_windows_fasta"]
    output:
        "results/cds_windows/assembly.fai"
    conda: "../envs/cds_windows.yaml"
    benchmark:
        "benchmarks/cds_windows/cds_windows_faidx.bench"
    shell:
        "samtools faidx --fai-idx {output} {input.fasta}"


rule cds_windows_cds_bed:
    input:
        gff3=config["cds_windows_gff3"]
    output:
        "results/cds_windows/cds.merged.bed"
    conda: "../envs/cds_windows.yaml"
    benchmark:
        "benchmarks/cds_windows/cds_windows_cds_bed.bench"
    shell:
        """
        awk -F '\\t' -v OFS='\\t' '$1 !~ /^#/ && $3 == "CDS" {{print $1, $4 - 1, $5}}' {input} \
            | sort -k1,1 -k2,2n \
            | bedtools merge -i - > {output}
        """


rule cds_windows_windows:
    input:
        fai="results/cds_windows/assembly.fai",
        cds="results/cds_windows/cds.merged.bed"
    output:
        "results/cds_windows/cds_windows_w{wlen}.bed"
    conda: "../envs/cds_windows.yaml"
    benchmark:
        "benchmarks/cds_windows/cds_windows_windows_w{wlen}.bench"
    shell:
        """
        cut -f1,2 {input.fai} > {output}.genome
        bedtools makewindows -g {output}.genome -w {wildcards.wlen} \
            | bedtools coverage -a - -b {input.cds} \
            | awk -v OFS='\\t' '{{print $1, $2, $3, $5}}' > {output}
        rm {output}.genome
        """