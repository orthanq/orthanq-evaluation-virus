#wrappers should be used once they are ready
rule orthanq_candidates:
    output:
        candidates_folder=directory("results/orthanq/candidates"),
        candidates="results/orthanq/candidates/candidates.vcf",
        reference_fasta="results/orthanq/candidates/reference.fasta",
        reference_fasta_idx="results/orthanq/candidates/reference.fasta.fai",
        sequences=expand("results/orthanq/candidates/sequences/{lineage}.fasta", lineage=orthanq_sequences), #all sequences should readily be there before simulation
    log:
        "logs/orthanq_candidates/candidates_virus.log",
    conda:
        "../envs/orthanq.yaml"
    priority: 50
    benchmark:    
        "benchmarks/orthanq_candidates/orthanq_candidates.tsv" 
    shell:
        "orthanq candidates virus --output {output.candidates_folder} 2> {log}"

rule orthanq_preprocess:
    input:
        candidates_folder="results/orthanq/candidates",
        reads=get_fastq_input
    output: "results/orthanq/preprocess/{sample}.bcf",
    log:
        "logs/orthanq_preprocess/{sample}.log",
    conda:
        "../envs/orthanq.yaml"
    benchmark:    
        "benchmarks/orthanq_preprocess/{sample}.tsv" 
    shell:
        "orthanq preprocess virus --candidates-folder {input.candidates_folder} --output {output} --reads {input.reads[0]} {input.reads[1]} 2> {log}"

#wrappers should be used once they are ready
rule orthanq_quantify:
    input:
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
        "orthanq call virus --candidates-folder {input.candidates_folder} --haplotype-calls {input.haplotype_calls} --prior {params.prior} --enable-equivalence-class-constraint --output {output.tsv} 2> {log}"

#the following rules uses orthanq from local build, comment out when necessary (testing purposes)
# rule orthanq_preprocess:
#     input:
#         candidates_folder="results/orthanq/candidates",
#         reads=get_fastq_input
#     output: "results/orthanq/preprocess/{sample}.bcf",
#     log:
#         "logs/orthanq_preprocess/{sample}.log",
#     conda:
#         "../envs/orthanq_cargo.yaml"
#     benchmark:    
#         "benchmarks/orthanq_preprocess/{sample}.tsv" 
#     shell:
#         "/home/hamdiyeuzuner/Documents/orthanq/target/release/orthanq preprocess virus --candidates-folder {input.candidates_folder} --output {output} --reads {input.reads[0]} {input.reads[1]} 2> {log}"

# #wrappers should be used once they are ready
# rule orthanq_quantify:
#     input:
#         candidates_folder="results/orthanq/candidates",
#         haplotype_calls="results/orthanq/preprocess/{sample}.bcf"
#     output:
#         tsv="results/orthanq/calls/{sample}/{sample}.tsv",
#         solutions="results/orthanq/calls/{sample}/viral_solutions.json",
#         final_solution="results/orthanq/calls/{sample}/final_solution.json",
#         lp_solution="results/orthanq/calls/{sample}/lp_solution.json",
#     log:
#         "logs/orthanq_call/{sample}.log"
#     conda:
#         "../envs/orthanq_cargo.yaml"
#     params:
#         prior="uniform"
#     resources: 
#         mem_mb=5000
#     benchmark:    
#         "benchmarks/orthanq_quantify/{sample}.tsv"
#     shell:
#         "/home/hamdiyeuzuner/Documents/orthanq/target/release/orthanq call virus --candidates-folder {input.candidates_folder} --haplotype-calls {input.haplotype_calls} --prior {params.prior} --enable-equivalence-class-constraint --output {output.tsv} 2> {log}"
