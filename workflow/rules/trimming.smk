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
        extra=fastp_params()
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
    threads: 30
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


##remove LTRs - required for the labmix sample
rule align_minimap2:
    input:
        ref="results/ref/hiv_reference_sequence.fasta",
        f1="results/sra/SRR961514_1.fastq.gz",
        f2="results/sra/SRR961514_2.fastq.gz"
    output:
        "results/minimap2_aligned/SRR961514_sorted.bam"
    log:
        "logs/prep_minimap2/SRR961514.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax sr {input.ref} {input.f1} {input.f2} | samtools view -bh | samtools sort > {output} 2> {log}"

rule index_bam:
    input:
        "results/minimap2_aligned/SRR961514_sorted.bam"
    output:
        "results/minimap2_aligned/SRR961514_sorted.bam.bai"
    log:
        "logs/labmix/index.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        "samtools index {input}"

rule exclude_repeat_reads:
    input:
        bam="results/minimap2_aligned/SRR961514_sorted.bam",
        bai="results/minimap2_aligned/SRR961514_sorted.bam.bai",
        bed="resources/labmix/repeat_regions.bed"
    output:
        "results/filtered/SRR961514_filtered.bam"
    log:
        "logs/exclude_repeat_reads.log"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -v -abam {input.bam} -b {input.bed} > {output} 2> {log}
        """

rule bam_to_fastq:
    input:
        "results/filtered/SRR961514_filtered.bam"
    output:
        fq1="results/fastq/SRR961514_filtered.fq1",
        fq2="results/fastq/SRR961514_filtered.fq2",
    log:
        "logs/bam_to_fastq.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        "samtools fastq -1 {output.fq1} -2 {output.fq2} -s /dev/null -n {input}  2> {log}"