rule plot_barchart:
    input:
        orthanq=expand("results/orthanq/calls/{sample}/{sample}.csv", sample=samples["sra"]),
        kallisto=expand("results/kallisto/quant_results_{sample}/abundance.tsv", sample=samples["sra"]),
        pangolin=expand("results/pangolin/{sample}.csv", sample=samples["sra"]),
        nextclade = expand("results/nextstrain/results/{sample}/nextclade.csv", sample=samples["sra"]),
        truth="resources/truth_both.csv"
    output:
        plot_svg="results/plots/stacked_barchart.svg",
        plot_html="results/plots/stacked_barchart.html",
        table="results/tables/all_tools_predictions.csv",
    log:
        "logs/plots/plot_stacked_barchart.log"
    conda:
        "../envs/altair.yaml"
    script:
        "../scripts/plot_stacked_barchart.py" 

