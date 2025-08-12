coverage=["100x", "1000x"]
rule scatter_plot:
    input:
        orthanq_prediction=expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/SimulatedSample{num}-{coverage}.csv", num=num_list, coverage=coverage),
        kallisto_prediction=expand("results/kallisto/quant_results_SimulatedSample{num}-{coverage}/abundance.tsv", num=num_list, coverage=coverage),
        simulation=expand("results/simulation_input/SimulatedSample{num}.csv", num=num_list),
    output:
        # orthanq_svg=expand("results/evaluation-pandemics/plots/orthanq/scatter_plot_{coverage}.svg",coverage=coverage),
        # orthanq_html=expand("results/evaluation-pandemics/plots/orthanq/scatter_plot_{coverage}.html",coverage=coverage),
        # kallisto_svg=expand("results/evaluation-pandemics/plots/kallisto/scatter_plot_{coverage}.svg",coverage=coverage),
        # kallisto_html=expand("results/evaluation-pandemics/plots/kallisto/scatter_plot_{coverage}.html",coverage=coverage),
        orthanq_svg_100x="results/evaluation-pandemics/plots/orthanq/scatter_plot_100x.svg",
        orthanq_html_100x=report("results/evaluation-pandemics/plots/orthanq/scatter_plot_100x.html",
        caption="../report/scatterplot_pandemic.rst",
        category="Simulation scatter plots", subcategory="100x", labels={
            "name": "orthanq",
            "type": "html"
         }),
        kallisto_svg_100x="results/evaluation-pandemics/plots/kallisto/scatter_plot_100x.svg",
        kallisto_html_100x=report("results/evaluation-pandemics/plots/kallisto/scatter_plot_100x.html",
        caption="../report/scatterplot_pandemic.rst",
        category="Simulation scatter plots", subcategory="100x", labels={
            "name": "kallisto",
            "type": "html"
        }),
        orthanq_svg_1000x="results/evaluation-pandemics/plots/orthanq/scatter_plot_1000x.svg",
        orthanq_html_1000x=report("results/evaluation-pandemics/plots/orthanq/scatter_plot_1000x.html",
        caption="../report/scatterplot_pandemic.rst",
        category="Simulation scatter plots", subcategory="1000x", labels={
            "name": "orthanq",
            "type": "html"
         }),
        kallisto_svg_1000x="results/evaluation-pandemics/plots/kallisto/scatter_plot_1000x.svg",
        kallisto_html_1000x=report("results/evaluation-pandemics/plots/kallisto/scatter_plot_1000x.html",
        caption="../report/scatterplot_pandemic.rst",
        category="Simulation scatter plots", subcategory="1000x", labels={
            "name": "kallisto",
            "type": "html"
        })
    log:
        "logs/evaluation-pandemics/scatter_plot.log"
    conda:
        "../envs/altair.yaml"
    benchmark:
        "benchmarks/scatter_plot/scatter_plot.tsv" 
    script:
        "../scripts/scatter_plot_pandemics_simulation.py"

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
        validation=expand("results/evaluation-pandemics/tables/abundant_lineage_validation_{coverage}.tsv",coverage=coverage)
    log:
        "logs/evaluation-pandemics/abundant_lineage.log"
    conda:
        "../envs/altair.yaml"
    benchmark:
        "benchmarks/evaluation/validation.tsv" 
    script:
        "../scripts/abundant_lineage_validation.py"

rule find_similarities:
    input:
        simulation=expand("results/simulation_input/SimulatedSample{num}.csv", num=num_list),
        candidates="results/orthanq/candidates/candidates.vcf",
    output:
        table="results/evaluation-pandemics/tables/lineage_similarities.csv"
    log:
        "logs/evaluation-pandemics/similarities.log"
    conda:
        "../envs/altair.yaml"
    script:
        "../scripts/find_similarities.py"