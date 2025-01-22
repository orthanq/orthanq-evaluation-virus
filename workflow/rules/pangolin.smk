
## following two rules are necessary for dynamic import of the config file ( and sample sheet that's being filled up during workflow execution)
## currently, we fill up the sample sheet of uncovar manually.

# rule create_sample_sheet_unicovar:
#     input:
#         fq1=unicovar_inputs[0],
#         fq2=unicovar_inputs[1],
#         template="config/pep/samples.csv",
#     output:
#         sample_sheet="resources/uncovar_config/pep/samples_filled.csv"
#     log:
#         "logs/create_sample_sheet_unicovar.log"
#     script:
#         "../scripts/create_unicovar_sheet.py"

# rule update_configs_uncovar:
#     input:
#         main_config="resources/uncovar_config/config.yaml",
#         pep_config="resources/uncovar_config/pep/config.yaml",
#         sample_sheet="resources/uncovar_config/pep/samples_filled.csv",
#     output:
#         new_main_config="resources/uncovar_config/config_filled.yaml",
#         new_pep_config="resources/uncovar_config/pep/pep_filled.yaml"
#     log:
#         "logs/uncovar/change_sample_sheet_path.log"
#     params:
#         new_pep_config_path=lambda w, output: output.new_pep_config,
#         sample_sheet_path=lambda w, input: input.sample_sheet
#     conda:
#         "../envs/yaml.yaml"
#     script:
#         "../scripts/change_sample_sheet_path.py"


module uncovar_pipeline:
    snakefile:
        github("IKIM-Essen/uncovar", path="workflow/Snakefile", commit="b16b3e3067e2d188227f80f1a637f76b5a8fdfab")
        # "/projects/koesterlab/orthanq/orthanq-evaluation-virus/uncovar/workflow/Snakefile"
    config: 
        config
    # prefix: "uncovar/"

use rule * from uncovar_pipeline as uncovar_*

# if config["simulate_pandemics"] or config["simulate_given"]:
#     rule ensure_uncovar_output:
#         input:
#             f"results/{DATE}/polishing/bcftools-illumina/{{sample}}-{{coverage}}.fasta"
# else:
#     rule ensure_uncovar_output:
#         input:
#             f"results/{DATE}/polishing/bcftools-illumina/{{sample}}.fasta"

##for simulation cases
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
