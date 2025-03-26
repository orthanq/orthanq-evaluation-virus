module uncovar_pipeline:
    snakefile:
        github("IKIM-Essen/uncovar", path="workflow/Snakefile", commit="b16b3e3067e2d188227f80f1a637f76b5a8fdfab")
    config: 
        config

use rule * from uncovar_pipeline as uncovar_*

rule pangolin:
    input:
        f"results/{DATE}/polishing/bcftools-illumina/{{sample}}.fasta"
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
