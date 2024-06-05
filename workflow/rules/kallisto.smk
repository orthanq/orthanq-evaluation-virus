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
    wrapper:
        "v3.11.0/bio/kallisto/index"

rule kallisto_quant:
    input:
        fastq=get_fastq_input,
        index="results/kallisto_index/sarscov2_all_nexstrain.idx",
    output:
        directory("results/kallisto/quant_results_{sample}"),
    params:
        extra="",
    log:
        "logs/kallisto_quant/{sample}.log",
    threads: 1
    wrapper:
        "v3.11.0/bio/kallisto/quant"
