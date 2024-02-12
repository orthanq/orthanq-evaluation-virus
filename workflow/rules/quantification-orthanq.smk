#wrappers should be used once they are ready
rule orthanq_candidates:
    output:
        directory("results/orthanq/candidates"),
        candidates="results/orthanq/candidates/candidates.vcf",
    log:
        "logs/orthanq_candidates/candidates_virus.log",
    conda:
        "../envs/orthanq.yaml"
    shell:
        "/home/uzuner/Documents/orthanq/target/debug/orthanq candidates virus --output {output} 2> {log}"

rule orthanq_preprocess:
    input:
        candidates="results/orthanq/candidates",
        reads=get_fastq_input
    output: "results/orthanq/preprocess/{sample}.bcf",
    log:
        "logs/orthanq_preprocess/{sample}.log",
    conda:
        "../envs/orthanq.yaml"
    shell:
        "/home/uzuner/Documents/orthanq/target/debug/orthanq preprocess virus --candidates-folder {input.candidates} --output {output} --reads {input.reads[0]} {input.reads[1]} 2> {log}"

rule orthanq_quantify:
    input:
        candidates="results/orthanq/candidates",
        haplotype_calls="results/orthanq/preprocess/{sample}.bcf"
    output:
        directory("results/orthanq/calls/{sample}")
    log:
        "logs/orthanq_call/{sample}.log"
    conda:
        "../envs/orthanq.yaml"
    params:
        prior="uniform"
    shell:
        "/home/uzuner/Documents/orthanq/target/debug/orthanq call virus --candidates-folder {input.candidates} --haplotype-calls {input.haplotype_calls} --prior {params.prior} --output {output} 2> {log}"
