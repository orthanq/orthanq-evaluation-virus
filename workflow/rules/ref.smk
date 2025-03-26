rule get_ref_sequence:
    output:
        "results/ref/reference_sequence.fasta",
    log:
        "logs/get_ref_sequence.log",
    conda:
        "../envs/entrez.yaml"
    params: accession_id=config["reference_accession"]
    benchmark:    
        "benchmarks/ncbi_datasets/get_ref_sequence.tsv" 
    shell:
        "esearch -db nucleotide -query {params} | efetch -format fasta > {output} 2> {log}"
