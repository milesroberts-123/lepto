# --- sequence tube map from the minimap2 hap1-vs-hap2 alignment ---
# Static tube-map-style figures of the {variant}_vs_{backbone} whole-genome
# alignment: `vg construct -a -R chrom:start-end` builds a region graph from
# the full backbone fasta + paf-derived variant VCF (no slicing; VCF and
# fasta coordinates line up natively), whose reference path is the backbone
# and whose alt paths carry the variant assembly's alleles, which `vg viz`
# renders as a tube map.
# Notes: regions must be kb-scale — vg viz draws 1 bp per layout unit, so a
# 500 kb region is an unreadable single-line strip (e2e-verified); graph
# coordinates restart at the region start; `tubemap_region_pad` gives the
# graph breathing room around the region of interest.

def tubemap_parse_regions():
    # chrom_start_end wildcard value -> (chrom, padded_start, padded_end)
    if "_tubemap_regions" not in globals():
        global _tubemap_regions
        _tubemap_regions = {}
        pad = int(config.get("tubemap_region_pad", 0))
        for spec in config.get("tubemap_regions", []):
            chrom, coords = spec.rsplit(":", 1)
            start, end = (int(x) for x in coords.split("-"))
            if start < 1 or end <= start:
                raise ValueError("malformed tubemap region: " + spec)
            _tubemap_regions[f"{chrom}_{start}_{end}"] = (
                chrom, max(1, start - pad), end + pad)
    return _tubemap_regions


def tubemap_chunk_coords(wildcards):
    return tubemap_parse_regions()[str(wildcards.region)]


rule vg_tubemap_graph:
    input:
        fasta="results/vg/{backbone}_prefixed.fasta",
        fai="results/vg/{backbone}_prefixed.fasta.fai",
        vcfgz="results/vg/{variant}_vs_{backbone}.vcf.gz",
        tbi="results/vg/{variant}_vs_{backbone}.vcf.gz.tbi"
    output:
        xg="results/tubemap/{variant}_vs_{backbone}/{region}/graph.xg"
    params:
        chrom=lambda wc: tubemap_chunk_coords(wc)[0],
        start=lambda wc: tubemap_chunk_coords(wc)[1],
        end=lambda wc: tubemap_chunk_coords(wc)[2]
    threads: 1
    conda: "../envs/vg.yaml"
    shell:
        """
        mkdir -p $(dirname {output.xg})
        trap 'rm -f {output.xg}.vg' EXIT
        vg construct -r {input.fasta} -v {input.vcfgz} -a \
            -R {params.chrom}:{params.start}-{params.end} > {output.xg}.vg
        vg index -x {output.xg} {output.xg}.vg
        """


rule vg_tubemap_viz:
    input:
        xg="results/tubemap/{variant}_vs_{backbone}/{region}/graph.xg"
    output:
        svg="results/tubemap/{variant}_vs_{backbone}/{region}/tube_map.svg"
    threads: 1
    conda: "../envs/vg.yaml"
    shell:
        "vg viz -x {input.xg} --out {output.svg}"