rule orthanq_candidates_hiv:
    input:
        genome=get_ref_seq_path(),
        genome_fai=get_ref_seq_path() + ".fai",
        lineages="results/hiv_genomes/hiv_final.fasta"
    output:
        candidates="results/orthanq/candidates/hiv/candidates.vcf",
        candidates_folder=directory("results/orthanq/candidates/hiv"),
    log:
        "logs/orthanq_candidates/hiv/candidates.log",
    conda:
        "../envs/orthanq_virus_eval.yaml"
    priority: 50
    benchmark:    
        "benchmarks/orthanq_candidates/hiv/orthanq_candidates.tsv" 
    shell:
        "orthanq candidates virus --genome {input.genome} --lineages {input.lineages} --output {output.candidates_folder} 2> {log}"

rule vembrane_subsample:
    input:
        "results/orthanq/candidates/hiv/candidates.vcf"
    output:
        "results/orthanq/candidates/hiv/subsampled_candidates.vcf"
    log:
        "logs/subsample_varaints.log"
    conda:
        "../envs/vembrane.yaml"
    shell:
        '''
        vembrane filter --context "import random; random.seed(42)" "random.random() < 0.1" {input} > {output} 2> {log}
        '''