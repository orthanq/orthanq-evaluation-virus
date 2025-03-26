rule fastqc:
    input:
        lambda wildcards: get_raw_fastq_input(wildcards)[0],  # Extracts the first (zeroth index) input file
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    resources:
        mem_mb=1024,
    wrapper:
        "v5.5.2/bio/fastqc"

rule fastp_pe_wo_trimming:
    input:
        sample=get_raw_fastq_input
    output:
        html="results/report/pe_wo_trimming/{sample}.html",
        json="results/report/pe_wo_trimming/{sample}.json"
    log:
        "logs/fastp/pe_wo_trimming/{sample}.log"
    params:
        extra=""
    threads: 10
    wrapper:
        "v5.9.0/bio/fastp"

rule fastp_pe:
    input:
        sample=get_raw_fastq_input
    output:
        trimmed=["results/fastp-trimmed/pe/{sample}.1.fastq", "results/fastp-trimmed/pe/{sample}.2.fastq"],
        json="results/fastp-trimmed/pe/{sample}.json",
        html="results/fastp-trimmed/pe/{sample}.html",
    log:
        "logs/fastp/pe/{sample}.log"
    params:
        extra="--detect_adapter_for_pe"
    threads: 10
    wrapper:
        "v5.9.0/bio/fastp"

rule bwa_index:
    input:
        "results/ref/reference_sequence.fasta"
    output:
        idx=multiext("results/ref/reference_sequence.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    wrapper:
        "v5.5.2/bio/bwa/index" 

rule align_bwa:
    input:
        ref="results/ref/reference_sequence.fasta",
        idx=rules.bwa_index.output,
        f1=lambda wildcards: get_raw_fastq_input(wildcards)[0],
        f2=lambda wildcards: get_raw_fastq_input(wildcards)[1]
    output:
        "results/bwa_aligned/{sample}.bam"
    log:
        "logs/prep_bwa/{sample}.log",
    threads: 20
    params: idx=lambda w, input: os.path.splitext(input.idx[0])[0]
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem -t {threads} {params.idx} {input.f1} {input.f2} "
        " | samtools view -bh | samtools sort > {output} 2> {log}"

rule ampliconclip:
    input:
        bam="results/bwa_aligned/{sample}.bam",
        primers=get_primers
    output:
        bam="results/clipped/{sample}.trimmed.bam"
    log:
        "logs/ampliconclip/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools ampliconclip -b {input.primers} {input.bam} -o {output.bam} 2> {log}"
 
rule extract_fastqs:
    input:
        "results/clipped/{sample}.trimmed.bam" 
    output:
        f1="results/extracted_reads/{sample}.forward.fastq",
        f2="results/extracted_reads/{sample}.reverse.fastq",
        singles="results/extracted_reads/{sample}.singles.fastq"
    log:
        "logs/extract_fastqs/{sample}.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        "samtools fastq -1 {output.f1} -2 {output.f2} -s {output.singles} <(samtools sort -n {input}) 2> {log}"