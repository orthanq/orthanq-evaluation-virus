rule scatter_plot:
    input:
        orthanq_prediction="results/orthanq/calls/SRR961514/SRR961514.csv",
    output:
        svg="results/evaluation-hiv/plots/scatter_plot.svg",
        html="results/evaluation-hiv/plots/scatter_plot.html"
    log:
        "logs/evaluation-hiv/scatter_plot.log"
    conda:
        "../envs/altair.yaml"
    script:
        "../scripts/scatter_plot_hiv.py"