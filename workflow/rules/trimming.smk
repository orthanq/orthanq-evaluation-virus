rule cutadapt:
    input:
        get_fastq_input,
    output:
        fastq1="results/trimmed/{sample}.1.fastq",
        fastq2="results/trimmed/{sample}.2.fastq",
        qc="results/trimmed/{sample}.qc.txt",
    params:
        #adapters obtained from: https://github.com/baymlab/wastewater_analysis/blob/9c82d9a63ca5718d1dcc7e73cc0e8931e161338a/auxiliary_data/adapters.fa#L3
        adapters="-g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -G AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT" #the reverse complement adapter (ADAPTER1_rc) is found in 5' ends in both pairs and the regular one (ADAPTER1) is in the 3`
        # extra="-q 30,30",
    log:
        "logs/cutadapt/{sample}.log",
    threads: 4  # set desired number of threads here
    wrapper:
        "v3.13.8/bio/cutadapt/pe"
