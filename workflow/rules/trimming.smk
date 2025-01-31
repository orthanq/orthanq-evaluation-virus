# rule cutadapt:
#     input:
#         get_fastq_input,
#     output:
#         fastq1="results/trimmed/{sample}.1.fastq",
#         fastq2="results/trimmed/{sample}.2.fastq",
#         qc="results/trimmed/{sample}.qc.txt",
#     params:
#         #adapters obtained from: https://github.com/baymlab/wastewater_analysis/blob/9c82d9a63ca5718d1dcc7e73cc0e8931e161338a/auxiliary_data/adapters.fa#L3
#         adapters="-g AGATGTGTATAAGAGACAG -a CTGTCTCTTATACACATCT -G AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT" #the reverse complement adapter (ADAPTER1_rc) is found in 5' ends in both pairs and the regular one (ADAPTER1) is in the 3`
#         # extra="-q 30,30",
#     log:
#         "logs/cutadapt/{sample}.log",
#     threads: 4  # set desired number of threads here
#     wrapper:
#         "v3.13.8/bio/cutadapt/pe"

rule trimmomatic:
    input:
        fastqs=get_fastq_input,
        adapters="resources/adapters.fa"
    output:
        f1="results/trimmomatic/{sample}_1.fastq",
        f1_s="results/trimmomatic/{sample}_s1.fastq",
        f2="results/trimmomatic/{sample}_2.fastq",
        f2_s="results/trimmomatic/{sample}_2_s2.fastq",
    log:
        "logs/trimmomatic/{sample}.log",
    threads: 20
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE -threads {threads} {input.fastqs} {output.f1} {output.f1_s} {output.f2} {output.f2_s} ILLUMINACLIP:{input.adapters}:2:30:10:2:keepBothReads " 
        "SLIDINGWINDOW:4:15 MINLEN:36 LEADING:3 TRAILING:3 > {log} 2>&1"

rule bwa_index:
    input:
        "results/ref/reference_sequence.fasta"
    output:
        idx=multiext("results/ref/reference_sequence", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index.log"
    params:
        prefix="results/ref/reference_sequence",
    wrapper:
        "v5.5.2/bio/bwa/index" 

rule align_bwa:
    input:
        ref="results/ref/reference_sequence.fasta",
        ref_idx=multiext("results/ref/reference_sequence", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        f1="results/trimmomatic/{sample}_1.fastq",
        f2="results/trimmomatic/{sample}_2.fastq",
    output:
        "results/trimmomatic/{sample}.bam"
    log:
        "logs/prep_bwa/{sample}.log",
    threads: 20
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem -L 100 -t {threads} {input.ref} {input.f1} {input.f2} "
        " | samtools view -bh | samtools sort > {output} 2> {log}"
    
rule ivar:
    input:
        bam="results/trimmomatic/{sample}.bam",
        primers="resources/primer.bedpe"
    output:
        dir=directory("results/ivar/{sample}"),
        file="results/ivar/{sample}.bam"
    log:
        "logs/ivar/{sample}.log"
    conda:
        "../envs/ivar.yaml"
    shell:
        "ivar trim -i {input.bam} -b {input.primers} -p {output.dir} -e > {log} 2>&1"

rule jvarkit:
    input:
        "results/ivar/{sample}.bam"
    output:
        "results/jvarkit/{sample}.bam"
    log:
        "logs/jvarkit/{sample}.log"
    conda:
        "../envs/jvarkit.yaml"
    shell:
        "jvarkit biostar84452 --samoutputformat BAM <(samtools sort {input} > {output} 2> {log}"

rule extract_fastqs:
    input:
        "results/jvarkit/{sample}.bam" 
    output:
        f1="results/jvarkit/{sample}.forward.fastq",
        f2="results/jvarkit/{sample}.reverse.fastq",
        singles="results/jvarkit/{sample}.singles.fastq"
    log:
        "logs/extract_fastqs/{sample}.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        "samtools fastq -1 {output.f1} -2 {output.f2} -s {output.singles} <(samtools sort -n {input}) 2> {log}"