rule kallisto_index:
    input:
        lineages=config["viral_lineages_fasta"]
    output:
        index=get_viral_lineages_path(),
    params:
        extra="",
    log:
        "logs/kallisto_index/kallisto_index.log",
    threads: config["kallisto_index_threads"]
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input}"

rule kallisto_quant:
    input:
        fastq=get_processed_fastq_input if not config["simulate_given"] and not config["simulate_pandemics"] else get_raw_fastq_input,
        index=get_viral_lineages_path(),
    output:
        dir=directory("results/kallisto/quant_results_{sample}"),
        tsv="results/kallisto/quant_results_{sample}/abundance.tsv",
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
        "kallisto quant -i {input.index} {input.fastq} -t {threads} -o {output.dir}"