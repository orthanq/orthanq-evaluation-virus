rule orthanq_candidates_generic:
    input:
        genome="results/ref/reference_sequence.fasta",
        lineages=config["viral_lineages_fasta"]
    output:
        candidates="results/orthanq/candidates/candidates.vcf",
        candidates_folder=directory("results/orthanq/candidates/"),
    log:
        "logs/orthanq_candidates/candidates_virus.log",
    conda:
        "../envs/orthanq_dev.yaml"
    priority: 50
    benchmark:    
        "benchmarks/orthanq_candidates/orthanq_candidates.tsv" 
    shell:
        "LD_LIBRARY_PATH=$CONDA_PREFIX/lib /projects/koesterlab/orthanq/orthanq/target/release/orthanq candidates virus --genome {input.genome} --lineages {input.lineages} --output {output.candidates_folder} 2> {log}"

rule orthanq_preprocess:
    input:
        candidates="results/orthanq/candidates/candidates.vcf", 
        reads=get_processed_fastq_input if not config["simulate_given"] and not config["simulate_pandemics"] else get_raw_fastq_input,
        genome="results/ref/reference_sequence.fasta"
    output: 
        bcf="results/orthanq/preprocess/{sample}.bcf",
        bam="results/orthanq/preprocess/{sample}_sorted.bam"
    log:
        "logs/orthanq_preprocess/{sample}.log",
    conda:
        "../envs/orthanq_dev.yaml"
    benchmark:    
        "benchmarks/orthanq_preprocess/{sample}.tsv" 
    shell:
        "LD_LIBRARY_PATH=$CONDA_PREFIX/lib /projects/koesterlab/orthanq/orthanq/target/release/orthanq preprocess virus --genome {input.genome} --candidates {input.candidates} --output {output.bcf} --reads {input.reads[0]} {input.reads[1]} --output-bam 2> {log}"

#wrappers should be used once they are ready
rule orthanq_quantify:
    input:
        haplotype_variants="results/orthanq/candidates/candidates.vcf", 
        haplotype_calls="results/orthanq/preprocess/{sample}.bcf"
    output:
        tsv="results/orthanq/calls/{sample}/{sample}.csv",
        solutions="results/orthanq/calls/{sample}/viral_solutions.json",
        final_solution="results/orthanq/calls/{sample}/final_solution.json",
        lp_solution="results/orthanq/calls/{sample}/lp_solution.json",
    log:
        "logs/orthanq_call/{sample}.log"
    conda:
        "../envs/orthanq_dev.yaml"
    params:
        prior="uniform"
    resources: 
        mem_mb=5000
    threads: 20
    benchmark:    
        "benchmarks/orthanq_quantify/{sample}.tsv"
    shell:
        "LD_LIBRARY_PATH=$CONDA_PREFIX/lib /projects/koesterlab/orthanq/orthanq/target/release/orthanq call virus --haplotype-variants {input.haplotype_variants} --haplotype-calls {input.haplotype_calls} --prior {params.prior} --output {output.tsv} 2> {log}"
