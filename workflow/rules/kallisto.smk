#quantify lineages by kallisto

rule kallisto_index:
    input:
        lineages=config["viral_lineages_fasta"]
    output:
        index="results/kallisto_index/viral_lineages.idx",
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
        fastq=get_trimmed_fastq_input,
        index="results/kallisto_index/viral_lineages.idx",
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
    benchmark:    
        "benchmarks/kallisto_quant/{sample}.tsv" 
    shell:
        "kallisto quant -i {input.index} {input.fastq} -o {output.dir}"