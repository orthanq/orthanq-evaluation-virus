rule get_ref_sequence:
    output:
        get_ref_seq_path(),
    log:
        "logs/get_ref_sequence.log",
    conda:
        "../envs/entrez.yaml"
    params: accession_id=config["reference_accession"]
    benchmark:    
        "benchmarks/ncbi_datasets/get_ref_sequence.tsv" 
    shell:
        "esearch -db nucleotide -query {params} | efetch -format fasta > {output} 2> {log}"

rule genome_faidx:
    input:
        ref,
    output:
        ref_fai,
    log:
        "logs/genome-faidx.log",
    wrapper:
        "v7.2.0/bio/samtools/faidx"