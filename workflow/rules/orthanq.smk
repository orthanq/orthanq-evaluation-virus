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
        "../envs/orthanq.yaml"
    priority: 50
    benchmark:    
        "benchmarks/orthanq_candidates/orthanq_candidates.tsv" 
    shell:
        "orthanq candidates virus generic --genome {input.genome} --lineages {input.lineages} --output {output.candidates_folder} 2> {log}"

rule orthanq_preprocess:
    input:
        candidates="results/orthanq/candidates/candidates.vcf", #just to make sure the file is generated
        candidates="results/orthanq/candidates",
        reads=get_trimmed_fastq_input,
        genome="results/ref/reference_sequence.fasta"
    output: "results/orthanq/preprocess/{sample}.bcf",
    log:
        "logs/orthanq_preprocess/{sample}.log",
    conda:
        "../envs/orthanq.yaml"
    benchmark:    
        "benchmarks/orthanq_preprocess/{sample}.tsv" 
    shell:
        "orthanq preprocess virus --genome {input.genome} --candidates {input.candidates} --output {output} --reads {input.reads[0]} {input.reads[1]} 2> {log}"

#wrappers should be used once they are ready
rule orthanq_quantify:
    input:
        candidates="results/orthanq/candidates/candidates.vcf", #just to make sure the file is generated
        candidates_folder="results/orthanq/candidates",
        haplotype_calls="results/orthanq/preprocess/{sample}.bcf"
    output:
        tsv="results/orthanq/calls/{sample}/{sample}.tsv",
        solutions="results/orthanq/calls/{sample}/viral_solutions.json",
        final_solution="results/orthanq/calls/{sample}/final_solution.json",
        lp_solution="results/orthanq/calls/{sample}/lp_solution.json",
    log:
        "logs/orthanq_call/{sample}.log"
    conda:
        "../envs/orthanq.yaml"
    params:
        prior="uniform"
    resources: 
        mem_mb=5000
    benchmark:    
        "benchmarks/orthanq_quantify/{sample}.tsv"
    shell:
        "orthanq call virus --candidates-folder {input.candidates_folder} --haplotype-calls {input.haplotype_calls} --enable-equivalence-class-constraint --prior {params.prior} --output {output.tsv} 2> {log}"
