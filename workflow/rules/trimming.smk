rule cutadapt:
    input:
        get_fastq_input,
    output:
        fastq1="results/trimmed/{sample}.1.fastq",
        fastq2="results/trimmed/{sample}.2.fastq",
        qc="results/trimmed/{sample}.qc.txt",
    params:
        adapters="",
        extra="-q 30,30",
    log:
        "logs/cutadapt/{sample}.log",
    threads: 4  # set desired number of threads here
    wrapper:
        "v3.13.8/bio/cutadapt/pe"
