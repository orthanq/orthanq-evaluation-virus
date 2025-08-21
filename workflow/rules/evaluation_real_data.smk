rule plot_barchart:
    input:
        orthanq=expand("results/orthanq/calls/{sample}/{sample}.csv", sample=samples["sra"]),
        kallisto=expand("results/kallisto/quant_results_{sample}/abundance.tsv", sample=samples["sra"]),
        pangolin=expand("results/pangolin/{sample}.csv", sample=samples["sra"]),
        nextclade = expand("results/nextstrain/results/{sample}/nextclade.csv", sample=samples["sra"]),
        truth="resources/truth_both.csv"
    output:
        plot_svg="results/evaluation-real-data/plots/stacked_barchart.svg",
        plot_html=report("results/evaluation-real-data/plots/stacked_barchart.html", 
            htmlindex="index.html", category="Co-infection evaluation", 
            subcategory="plot",labels={
            "name": "stacked bar chart",
            "type": "html"
        },
        caption="../report/stackedbarchart.rst"),
        table="results/evaluation-real-data/tables/all_tools_predictions.csv",
    log:
        "logs/evaluation-real-data/plot_stacked_barchart.log"
    conda:
        "../envs/altair.yaml"
    script:
        "../scripts/plot_stacked_barchart.py" 

rule datavzrd_tool_predictions:
    input:
        config="resources/datavzrd/tool_predictions.yaml",
        tool_predictions="results/evaluation-real-data/tables/all_tools_predictions.csv",
    output:
        report(
            directory("results/evaluation-real-data/datavzrd-report/all_tool_predictions"),
            htmlindex="index.html", category="Co-infection evaluation", subcategory="table", labels={
            "name": "all tool predictions",
            "type": "table"
        }),
    log:
        "logs/datavzrd/tool_predictions.log",
    wrapper:
        "v7.2.0/utils/datavzrd"