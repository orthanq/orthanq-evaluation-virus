#wrappers should be used once they are ready
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
        "/home/hamdiyeuzuner/Documents/orthanq/target/release/orthanq preprocess virus --candidates-folder {input.candidates_folder} --output {output} --reads {input.reads[0]} {input.reads[1]} 2> {log}"

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
    benchmark:    
        "benchmarks/orthanq_quantify/{sample}.tsv"
    shell:
        "/home/hamdiyeuzuner/Documents/orthanq/target/release/orthanq call virus --candidates-folder {input.candidates_folder} --haplotype-calls {input.haplotype_calls} --prior {params.prior} --output {output.tsv} 2> {log}"
