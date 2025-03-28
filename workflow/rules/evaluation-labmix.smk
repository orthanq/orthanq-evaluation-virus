rule scatter_plot:
    input:
        orthanq_prediction="results/orthanq/calls/SRR961514/SRR961514.csv",
        kallisto_prediction="results/kallisto/quant_results_SRR961514/abundance.tsv"
    output:
        orthanq_svg="results/evaluation-hiv/plots/orthanq/scatter_plot.svg",
        orthanq_html="results/evaluation-hiv/plots/orthanq/scatter_plot.html",
        kallisto_svg="results/evaluation-hiv/plots/kallisto/scatter_plot.svg",
        kallisto_html="results/evaluation-hiv/plots/kallisto/scatter_plot.html",
    log:
        "logs/evaluation-hiv/scatter_plot.log"
    conda:
        "../envs/altair.yaml"
    script:
        "../scripts/scatter_plot_hiv.py"