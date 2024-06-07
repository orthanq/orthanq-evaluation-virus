rule get_dataset:
    output: 
        directory("results/nextstrain/sarscov2_dataset")
    params: 
        data="nextstrain/sars-cov-2/wuhan-hu-1/orfs" #this is the name of the dataest provided by nextclade
    conda: 
        "../envs/nextclade.yaml"
    log: 
        "logs/nextclade/get_dataset/get_dataset.log",
    shell:
        "nextclade dataset get --name {params} --output-dir {output} 2> {log}"

rule run_nextclade:
    input:
        sample_fasta="results/consensus/{sample}.fasta",
        dataset="results/nextstrain/sarscov2_dataset"
    output: 
        directory("results/nextstrain/results/{sample}")
    conda: 
        "../envs/nextclade.yaml"
    log: 
        "logs/nextclade/run_nextclade/{sample}.log",
    shell:
        "nextclade run --input-dataset {input.dataset} --output-all={output} {input.sample_fasta}"