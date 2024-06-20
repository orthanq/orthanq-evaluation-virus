rule get_sra:
    output:
        "results/sra/{accession}_1.fastq.gz",
        "results/sra/{accession}_2.fastq.gz",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "v2.3.2/bio/sra-tools/fasterq-dump"