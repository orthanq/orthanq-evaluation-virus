#quantify lineages by kallisto

#Step 0: concatenate all fastas into one
rule concat_fastas:
    input:
        expand("results/orthanq/candidates/sequences/{lineage}.fasta", lineage=orthanq_sequences)
    output:
        "results/orthanq/candidates/sequences/sarscov2_all_nextstrain.fasta"
    log:
        "logs/concat_fastas/concatenation.log"
    shell:
        "cat {input} > {output} 2> {log}"

rule kallisto_index:
    input:
        fasta="results/orthanq/candidates/sequences/sarscov2_all_nextstrain.fasta",
    output:
        index="results/kallisto_index/sarscov2_all_nexstrain.idx",
    params:
        extra="",  # optional parameters
    log:
        "logs/kallisto_index/kallisto_index.log",
    threads: config["kallisto_index_threads"]
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input}"

rule kallisto_quant:
    input:
        fastq=get_fastq_input,
        index="results/kallisto_index/sarscov2_all_nexstrain.idx",
    output:
        dir=directory("results/kallisto/quant_results_{sample}"),
        tsv="results/kallisto/quant_results_{sample}/abundance.tsv"
    params:
        extra="",
    log:
        "logs/kallisto_quant/{sample}.log",
    threads: config["kallisto_quant_threads"]
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto quant -i {input.index} {input.fastq} -o {output.dir}"