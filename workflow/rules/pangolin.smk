#align reads with bwa, sort and call variants and create consensus fasta sequences.
#Then use pangolin for finding the abundant lineage

# align reads by bwa
rule bwa_index:
    input:
        "results/orthanq/candidates/reference.fasta"
    output:
        idx=multiext("results/bwa-index/MN908947", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index/MN908947.log"
    params:
        prefix=lambda w, output: os.path.splitext(output[0])[0],
        algorithm="bwtsw",
    wrapper:
        "v2.0.0/bio/bwa/index" 

rule bwa_mem:
    input:
        reads = get_fastq_input,
        idx = multiext("results/bwa-index/MN908947", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "results/bwa_alignment/{sample}_mapped.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:    
        "benchmarks/bwa_mem/{sample}.tsv"
    params:
        index="results/bwa-index/MN908947",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",             
        sort_order="coordinate", 
    threads: config["bwa_threads"]
    wrapper:
        "v2.0.0/bio/bwa/mem"

rule samtools_index_bwa:
    input:
        "results/bwa_alignment/{sample}_mapped.bam"
    output:
        "results/bwa_alignment/{sample}_mapped.bai"
    log:
        "logs/samtools_index_bwa/{sample}.log"
    threads: 10
    wrapper:
        "v2.0.0/bio/samtools/index"

#samtools mpileup does not support bcf/vcf output anymore. it's deprecated. use bcftools mpileup instead.
rule bcftools_mpileup:
    input:
        alignments=["results/bwa_alignment/{sample}_mapped.bam"],
        ref="results/orthanq/candidates/reference.fasta",  # this can be left out if --no-reference is in options
        index="results/orthanq/candidates/reference.fasta.fai",
    output:
        pileup="results/bcftools_mpileup/{sample}.pileup.bcf",
    params:
        uncompressed_bcf=False,
        # extra="--max-depth 100 --min-BQ 15",
    log:
        "logs/bcftools_mpileup/{sample}.log",
    wrapper:
        "v3.11.0/bio/bcftools/mpileup"

rule bcftools_call:
    input:
        pileup="results/bcftools_mpileup/{sample}.pileup.bcf",
    output:
        calls="results/bcftools_call/{sample}.vcf.gz",
    params:
        uncompressed_bcf=False,
        caller="-m",  # valid options include -c/--consensus-caller or -m/--multiallelic-caller
        extra="-v -Oz --write-index --ploidy 1 --prior 0.001"
    log:
        "logs/bcftools_call/{sample}.log",
    wrapper:
        "v3.11.0/bio/bcftools/call"

rule generate_consensus:
    input:
        vcf="results/bcftools_call/{sample}.vcf.gz",
        ref="results/orthanq/candidates/reference.fasta"
    output:
        "results/consensus/{sample}.fasta"
    log:
        "logs/bcftools_consensus/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    benchmark:
        "benchmarks/evaluation/{sample}.tsv" 
    shell:
        """
        bcftools consensus -f {input.ref} {input.vcf} > {output} 2> {log}
        """

rule pangolin:
    input:
        "results/consensus/{sample}.fasta",
    output:
        "results/pangolin/{sample}.csv"
    log:
        "logs/pangolin/{sample}.log"
    conda:
        "../envs/pangolin.yaml"
    benchmark:
        "benchmarks/evaluation/{sample}.tsv" 
    shell:
        "pangolin {input} --outfile {output} 2> {log}"