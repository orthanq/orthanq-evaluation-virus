# rule art_simulation:
#     input:
#         "resources/lineages/{lineage}.fasta"
#     output:
#         "results/art/{lineage}_1.fq",
#         "results/art/{lineage}_2.fq"
#     log:
#         "logs/art/{lineage}.log",
#     #threads: config["threads"]["art"]
#     conda:
#         "../envs/art.yaml"
#     params: config["f_coverage"]
#     shell:
#         "art_illumina -ss HS25 -i {input} -p -l 150 -s 100 -m 200 -f {params} --noALN --rndSeed 31303889 -o"
#         " results/art/{wildcards.lineage}_ 2> {log}"

# rule wgsim:
#     input:
#         ref="resources/lineages/{lineage}.fasta"
#     output:
#         read1="results/art/{lineage}_1.fq",
#         read2="results/art/{lineage}_2.fq"
#     log:
#         "logs/wgsim/{lineage}.log"
#     params:
#         "-X 0 -R 0 -r 0 -h"
#     wrapper:
#         "v3.3.6/bio/wgsim"


rule mason:
    input:
        ref="resources/lineages/{lineage}.fasta"
    output:
        read1="results/art/{lineage}_1.fq",
        read2="results/art/{lineage}_2.fq"
    log:
        "logs/mason/{lineage}.log"
    conda:
        "../envs/mason.yaml"
    benchmark:    
        "benchmarks/mason/{lineage}.tsv" 
    shell:
        "mason_simulator -ir {input.ref} -n 10000 -o {output.read1} -or {output.read2} --read-name-prefix {wildcards.lineage}"

rule get_fractions:
    input:
        fq1="results/art/{lineage}_1.fq",
        fq2="results/art/{lineage}_2.fq" 
    output:
        out_fq1="results/fractions/{sample}-{lineage}-{num}_1.fq",
        out_fq2="results/fractions/{sample}-{lineage}-{num}_2.fq"
    log:
        "logs/seqtk/{sample}-{lineage}-{num}.log",
    params:
        "{num}"
    conda:
        "../envs/seqtk.yaml"
    benchmark:    
        "benchmarks/get_fractions/{sample}-{lineage}-{num}.tsv" 
    shell:
        "seqtk sample -s100 {input.fq1} {params} > {output.out_fq1}; "
        "seqtk sample -s100 {input.fq2} {params} > {output.out_fq2} 2> {log}"
 
rule concat_fractions: 
    input:
        fq1=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_1.fq",
            zip,
            lineage=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['lineage'],
            num=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['num_reads']
            ),
        fq2=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_2.fq",
            zip,
            lineage=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['lineage'],
            num=simulated_sample.loc[simulated_sample['sample_name'] == wc.sample]['num_reads']
            ),
    output:
        out_fq1="results/mixed/{sample}_1.fq",
        out_fq2="results/mixed/{sample}_2.fq"
    log:
        "logs/mixed/{sample}.log",
    benchmark:    
        "benchmarks/concat_fractions/{sample}.tsv" 
    shell:
        "cat {input.fq1} > {output.out_fq1}; "
        "cat {input.fq2} > {output.out_fq2}"

