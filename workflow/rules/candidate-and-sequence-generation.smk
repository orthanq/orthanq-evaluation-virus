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
        "/home/hamdiyeuzuner/Documents/orthanq/target/release/orthanq candidates virus --output {output.candidates_folder} 2> {log}"
