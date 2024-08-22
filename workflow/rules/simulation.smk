rule download_country_data:
    output:
        sequence_counts_usa_gz="results/download_sequence_counts/SeqCountsUSA.tsv.gz",
        sequence_counts_usa="results/download_sequence_counts/SeqCountsUSA.tsv"
    log:
        "logs/download_country_data/download_country_data.log"
    shell:
        "wget -c https://data.nextstrain.org/files/workflows/forecasts-ncov/gisaid/pango_lineages/usa.tsv.gz -O {output.sequence_counts_usa_gz} && "
        " gzip -d {output.sequence_counts_usa_gz} -c > {output.sequence_counts_usa}"

rule create_simulation_input:
    input:
        sequence_counts_usa="results/download_sequence_counts/SeqCountsUSA.tsv"
    output:
        # sequence_counts_usa_raw="results/simulation_input/SeqCountsUSAraw.csv",
        sequence_counts_usa_final="results/simulation_input/SeqCountsUSAFinal.csv",
    log:
        "logs/create_simulation_input/create_simulation_input.log"
    conda:
        "../envs/data_wrangle.yaml"
    benchmark:
        "benchmarks/create_simulation_input/create_simulation_input.tsv" 
    script:
        "../scripts/create_simulation_input.py"

checkpoint create_sample_compositions:
    input:
        sequence_counts_usa_final="results/simulation_input/SeqCountsUSAFinal.csv",
    output:
        synthetic_samples = expand("results/simulation_input/SimulatedSample{num}.csv", num=num_list)
    log:
        "logs/create_sample_compositions/create_sample_compositions.log"
    benchmark:
        "benchmarks/create_sample_compositions/create_sample_compositions.tsv" 
    params: n_of_samples=config["number_of_samples"]
    script:
        "../scripts/create_sample_compositions.py"

rule mason:
    input:
        ref="resources/lineages/{lineage}.fasta"
    output:
        read1="results/mason/{lineage}_1.fq",
        read2="results/mason/{lineage}_2.fq",
        alignment="results/mason/{lineage}.sam"
    log:
        "logs/mason/{lineage}.log"
    conda:
        "../envs/mason.yaml"
    benchmark:    
        "benchmarks/mason/{lineage}.tsv" 
    shell:
        "mason_simulator -ir {input.ref} -n 50000 --illumina-read-length 150 -o {output.read1} -or {output.read2} --out-alignment {output.alignment} --read-name-prefix {wildcards.lineage} 2> {log}"

rule get_fractions:
    input:
        fq1="results/mason/{lineage}_1.fq",
        fq2="results/mason/{lineage}_2.fq" 
    output:
        out_fq1="results/fractions/{sample}-{lineage}-{num}_1.fq",
        out_fq2="results/fractions/{sample}-{lineage}-{num}_2.fq"
    log:
        "logs/seqtk/{sample}-{lineage}-{num}.log",
    params:
        "{num}"
    conda:
        "../envs/seqtk.yaml"
    benchmark:    
        "benchmarks/get_fractions/{sample}-{lineage}-{num}.tsv" 
    shell:
        "seqtk sample -s100 {input.fq1} {params} > {output.out_fq1}; "
        "seqtk sample -s100 {input.fq2} {params} > {output.out_fq2} 2> {log}"

rule concat_fractions: 
    input:
        # fq1=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_1.fq",
        #     zip,
        #     lineage=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['lineage'],
        #     num=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['num_reads']
        #     ),
        # fq2=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_2.fq",
        #     zip,
        #     lineage=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['lineage'],
        #     num=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['num_reads']
        #     ),
        fq1=lambda wildcards: get_concat_fractions_input(wildcards)[0],
        fq2=lambda wildcards: get_concat_fractions_input(wildcards)[1]

        # fq1="results/fractions/{{sample}}-{lineage}-{num}_1.fq",
        # fq2="results/fractions/{{sample}}-{lineage}-{num}_2.fq"
    output:
        out_fq1="results/mixed/{sample}_1.fastq",
        out_fq2="results/mixed/{sample}_2.fastq"
    log:
        "logs/mixed/{sample}.log",
    benchmark:    
        "benchmarks/concat_fractions/{sample}.tsv" 
    shell:
        "cat {input.fq1} > {output.out_fq1}; "
        "cat {input.fq2} > {output.out_fq2}"

