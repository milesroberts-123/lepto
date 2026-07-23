rule panther_download:
    output:
        "resources/panther/PANTHER19.0_hmmscoring.tgz"
    params:
        url="https://data.pantherdb.org/ftp/panther_library/current_release/PANTHER19.0_hmmscoring.tgz"
    shell:
        """
        mkdir -p resources/panther
        curl -L -o {output} {params.url}
        """


rule panther_extract_hmm:
    input:
        "resources/panther/PANTHER19.0_hmmscoring.tgz"
    output:
        "results/panther/binHmm.h3m"
    shell:
        """
        mkdir -p results/panther
        tar -xzf {input} -C results/panther --strip-components=7 \
            target/famlib/rel/PANTHER19.0_altVersion/hmmscoring/PANTHER19.0/globals/binHmm.h3m
        """


rule panther_split_fasta:
    input:
        config["panther_input_fasta"]
    output:
        expand("results/panther/batches/batch_{batch}.faa",
               batch=range(1, config["panther_batches"] + 1))
    params:
        n_batches=config["panther_batches"]
    run:
        import os

        records = []
        with open(input[0]) as f:
            header = None
            seq_lines = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header is not None:
                        records.append((header, "".join(seq_lines)))
                    header = line
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if header is not None:
                records.append((header, "".join(seq_lines)))

        n = len(records)
        n_batches = params.n_batches
        batch_size = (n + n_batches - 1) // n_batches

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        for i in range(n_batches):
            start = i * batch_size
            end = min(start + batch_size, n)
            with open(output[i], "w") as out:
                for header, seq in records[start:end]:
                    out.write(f"{header}\n{seq}\n")


rule panther_hmmsearch:
    input:
        hmm="results/panther/binHmm.h3m",
        fasta="results/panther/batches/batch_{batch}.faa"
    output:
        tblout="results/panther/hmmsearch/batch_{batch}.tblout",
        domtblout="results/panther/hmmsearch/batch_{batch}.domtblout",
        out="results/panther/hmmsearch/batch_{batch}.out"
    conda: "../envs/hmmer.yaml"
    params:
        cpu=config["panther_hmmsearch_cpu"],
        evalue=config["panther_hmmsearch_evalue"],
        dome=config["panther_hmmsearch_dome"],
        ince=config["panther_hmmsearch_ince"],
        incdome=config["panther_hmmsearch_incdome"]
    shell:
        """
        mkdir -p results/panther/hmmsearch
        if [ -s {input.fasta} ]; then
            hmmsearch --cpu {params.cpu} \
                -E {params.evalue} \
                --domE {params.dome} \
                --incE {params.ince} \
                --incdomE {params.incdome} \
                --tblout {output.tblout} \
                --domtblout {output.domtblout} \
                -o {output.out} \
                {input.hmm} {input.fasta}
        else
            touch {output.tblout} {output.domtblout} {output.out}
        fi
        """


rule panther_all:
    input:
        expand("results/panther/hmmsearch/batch_{batch}.tblout",
               batch=range(1, config["panther_batches"] + 1))
