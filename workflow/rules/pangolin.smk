# #align reads with bwa, sort and call variants and create consensus fasta sequences.
# #Then use pangolin for finding the abundant lineage

# # align reads by bwa
# rule bwa_index:
#     input:
#         "results/orthanq/candidates/reference.fasta"
#     output:
#         idx=multiext("results/bwa-index/MN908947", ".amb", ".ann", ".bwt", ".pac", ".sa")
#     log:
#         "logs/bwa_index/MN908947.log"
#     params:
#         prefix=lambda w, output: os.path.splitext(output[0])[0],
#         algorithm="bwtsw",
#     wrapper:
#         "v2.0.0/bio/bwa/index" 

# rule bwa_mem:
#     input:
#         reads = get_fastq_input,
#         idx = multiext("results/bwa-index/MN908947", ".amb", ".ann", ".bwt", ".pac", ".sa")
#     output:
#         "results/bwa_alignment/{sample}_mapped.bam"
#     log:
#         "logs/bwa_mem/{sample}.log"
#     benchmark:    
#         "benchmarks/bwa_mem/{sample}.tsv"
#     params:
#         index="results/bwa-index/MN908947",
#         extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
#         sorting="samtools",             
#         sort_order="coordinate", 
#     threads: config["bwa_threads"]
#     wrapper:
#         "v2.0.0/bio/bwa/mem"

# rule samtools_index_bwa:
#     input:
#         "results/bwa_alignment/{sample}_mapped.bam"
#     output:
#         "results/bwa_alignment/{sample}_mapped.bai"
#     log:
#         "logs/samtools_index_bwa/{sample}.log"
#     threads: 10
#     wrapper:
#         "v2.0.0/bio/samtools/index"

# #samtools mpileup does not support bcf/vcf output anymore. it's deprecated. use bcftools mpileup instead.
# rule bcftools_mpileup:
#     input:
#         alignments=["results/bwa_alignment/{sample}_mapped.bam"],
#         ref="results/orthanq/candidates/reference.fasta",  # this can be left out if --no-reference is in options
#         index="results/orthanq/candidates/reference.fasta.fai",
#     output:
#         pileup="results/bcftools_mpileup/{sample}.pileup.bcf",
#     params:
#         uncompressed_bcf=False,
#         # extra="--max-depth 100 --min-BQ 15",
#     log:
#         "logs/bcftools_mpileup/{sample}.log",
#     wrapper:
#         "v3.11.0/bio/bcftools/mpileup"

# rule bcftools_call:
#     input:
#         pileup="results/bcftools_mpileup/{sample}.pileup.bcf",
#     output:
#         calls="results/bcftools_call/{sample}.vcf.gz",
#     params:
#         uncompressed_bcf=False,
#         caller="-m",  # valid options include -c/--consensus-caller or -m/--multiallelic-caller
#         extra="-v -Oz --write-index --ploidy 1 --prior 0.001"
#     log:
#         "logs/bcftools_call/{sample}.log",
#     wrapper:
#         "v3.11.0/bio/bcftools/call"

# rule vembrane_filter:
#     input:
#         vcf="results/orthanq/preprocess/{sample}.bcf",
#     output:
#         vcf="results/orthanq/filtered_variants/{sample}.bcf"
#     params:
#         expression="FORMAT['AF'][SAMPLES[0]] > 0.5",
#         extra=""
#     log:
#         "logs/{sample}.log"
#     wrapper:
#         "v3.12.0/bio/vembrane/filter"

# rule bgzip:
#     input:
#         "results/orthanq/filtered_variants/{sample}.bcf",
#     output:
#         "results/orthanq/filtered_variants/{sample}.bcf.gz",
#     params:
#         extra="", # optional
#     threads: 1
#     log:
#         "logs/bgzip/{sample}.log",
#     wrapper:
#         "v3.12.0/bio/bgzip"


# rule tabix:
#     input:
#         "results/orthanq/filtered_variants/{sample}.bcf.gz",
#     output:
#         "results/orthanq/filtered_variants/{sample}.bcf.gz.tbi",
#     log:
#         "logs/tabix/{sample}.log",
#     params:
#         # pass arguments to tabix (e.g. index a vcf)
#         "-p bcf",
#     wrapper:
#         "v3.12.0/bio/tabix/index"


# rule generate_consensus:
#     input:
#         vcf="results/orthanq/filtered_variants/{sample}.bcf.gz",
#         vcf_tbi="results/orthanq/filtered_variants/{sample}.bcf.gz.tbi",
#         ref="results/orthanq/candidates/reference.fasta"
#     output:
#         "results/consensus/{sample}.fasta"
#     log:
#         "logs/bcftools_consensus/{sample}.log",
#     conda:
#         "../envs/bcftools.yaml"
#     benchmark:
#         "benchmarks/evaluation/{sample}.tsv" 
#     shell:
#         """
#         bcftools consensus -f {input.ref} {input.vcf} > {output} 2> {log}
#         """

# rule ivar_consensus:
#     input:
#         bam="results/bwa_alignment/{sample}_mapped.bam",
#         ref="results/orthanq/candidates/reference.fasta"
#     output:
#         "results/ivar_consensus/{sample}.fa"
#     log:
#         "logs/ivar_consensus/{sample}.log"
#     conda:
#         "../envs/ivar.yaml"
#     benchmark:
#         "benchmarks/ivar/{sample}.tsv"
#     shell:
#         "samtools mpileup -A -Q 0 {input.bam} | ivar consensus -p {output} -t 0.7 2> {log}"

# rule clone_uncovar_pipeline:
#     output:
#         touch("results/clone_uncovar_done.txt"),
#         uncovar_dir=directory("uncovar"),
#     shell:
#         "git clone https://github.com/IKIM-Essen/uncovar.git {output.uncovar_dir}"

rule create_sample_sheet_unicovar:
    input:
        fq1=unicovar_inputs[0],
        fq2=unicovar_inputs[1],
        template=UNCOVAR_SAMPLE_SHEET,
        # uncovar_check="results/clone_uncovar_done.txt"
    output:
        sample_sheet="uncovar/config/pep/samples2.csv"
    log:
        "logs/create_sample_sheet_unicovar.log"
    script:
        "../scripts/create_unicovar_sheet.py"

rule update_configs_uncovar:
    input:
        config=UNCOVAR_CONFIG,
        pepfile=UNCOVAR_PEP_CONFIG,
        sample_sheet="uncovar/config/pep/samples2.csv",
        # uncovar_check="results/clone_uncovar_done.txt"
    output:
        main_config="uncovar/config/config2.yaml",
        pep_config="uncovar/config/pep/config2.yaml"
    log:
        "logs/uncovar/change_sample_sheet_path.log"
    params:
        pepfile_path=lambda w, output: output.pep_config,
        sample_sheet_path=lambda w, input: input.sample_sheet
    conda:
        "../envs/yaml.yaml"
    script:
        "../scripts/change_sample_sheet_path.py"

ruleorder: touch_uncovar_results > execute_uncovar > transfer_results
rule touch_uncovar_results:  
    output:
        output_files="uncovar/results/{date}/polishing/bcftools-illumina/{sample}.fasta"
    shell:
        "touch {output.output_files}"

rule execute_uncovar:
    input:
        sample_sheet="uncovar/config/pep/samples2.csv",
        main_config="uncovar/config/config2.yaml",
        pep_config="uncovar/config/pep/config2.yaml",
        # uncovar_results_touch="results/touch_uncovar_results_done.txt"
    output:
        touch("results/pangolin/unconvar_execution_done.txt")
    params: cores=2,
        output_files=expand("results/{date}/polishing/bcftools-illumina/SimulatedSample{num}.fasta", date=DATE, num=num_list)
    log:
        "logs/uncovar/execute_workflow.log"
    shell:
        "cd uncovar && "
        " snakemake -p --configfile config/config2.yaml --sdm conda --cores {params.cores} {params.output_files} --rerun-incomplete"

rule transfer_results:
    input:
        uncovar_results="uncovar/results/{date}/polishing/bcftools-illumina/{sample}.fasta",
        uncovar_execution_check="results/pangolin/unconvar_execution_done.txt" #required for the rule execute_uncovar to be executed before transfer_results
    output:
        "results/{date}/polishing/bcftools-illumina/{sample}.fasta"
    shell:
        "cp {input.uncovar_results} {output}"

rule pangolin:
    input:
        f"results/{DATE}/polishing/bcftools-illumina/{{sample}}.fasta"
    output:
        "results/pangolin/{sample}_{date}.csv"
    log:
        "logs/pangolin/{sample}_{date}.log"
    conda:
        "../envs/pangolin.yaml"
    benchmark:
        "benchmarks/evaluation/{sample}_{date}.tsv" 
    shell:
        "pangolin {input} --outfile {output} 2> {log}"