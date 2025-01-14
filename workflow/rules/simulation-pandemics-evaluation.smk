coverage=["100x", "1000x"]
rule scatter_plot:
    input:
        orthanq_prediction=expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/SimulatedSample{num}-{coverage}.csv", num=num_list, coverage=coverage),
        simulation=expand("results/simulation_input/SimulatedSample{num}.csv", num=num_list),
    output:
        plot=expand("results/evaluation/scatter_plot_{coverage}.svg",coverage=coverage),
    log:
        "logs/scatter_plot/scatter_plot.log"
    conda:
        "../envs/altair.yaml"
    benchmark:
        "benchmarks/scatter_plot/scatter_plot.tsv" 
    script:
        "../scripts/scatter_plot.py"

#the downloaded file is updated daily by UCSC
rule download_clade_to_pangolin:
    output:
        "results/clade_to_lineage/clade_to_lineages.tsv.gz"
    log:
        "logs/download_and_prepare_clade_to_pangolin/download.log"
    shell:
        "wget -c https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.metadata.tsv.gz -O {output} " 

rule prepare_clade_to_pangolin:
    input:
        "results/clade_to_lineage/clade_to_lineages.tsv.gz"
    output: 
        "results/clade_to_lineage/clade_to_lineages.tsv"
    log:
        "logs/download_and_prepare_clade_to_pangolin/prepare.log"
    shell:
        "gunzip -c {input} | cut -f 10,11 | sort -u > {output} 2> {log}" 

rule abundant_lineage_validation:
    input:
        simulation=expand("results/simulation_input/SimulatedSample{num}.csv", num=num_list),
        clade_to_lineage="results/clade_to_lineage/clade_to_lineages.tsv",
        orthanq=expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/SimulatedSample{num}-{coverage}.csv", num=num_list, coverage=coverage),
        pangolin=expand("results/pangolin/SimulatedSample{num}-{coverage}.csv", num=num_list, coverage=coverage),
        nextclade=expand("results/nextstrain/results/SimulatedSample{num}-{coverage}/nextclade.tsv", num=num_list, coverage=coverage),
        kallisto=expand("results/kallisto/quant_results_SimulatedSample{num}-{coverage}/abundance.tsv", num=num_list, coverage=coverage),
    output:
        validation=expand("results/evaluation/validation_{coverage}.tsv",coverage=coverage)
    log:
        "logs/evaluation/validation.log"
    conda:
        "../envs/altair.yaml"
    benchmark:
        "benchmarks/evaluation/validation.tsv" 
    script:
        "../scripts/abundant_lineage_validation.py"