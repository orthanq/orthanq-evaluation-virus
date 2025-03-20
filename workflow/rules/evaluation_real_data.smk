rule scatter_plot:
    input:
        orthanq_prediction=expand("results/orthanq/calls/{sample}/{sample}.csv",sample=samples["sra"]),
        truth="resources/PRJNA809680_truth.csv",
    output:
        plot="results/evaluation/scatter_plot_PRJNA809680.svg"
    log:
        "logs/scatter_plot/PRJNA809680/scatter_plot.log"
    conda:
        "../envs/altair.yaml"
    script:
        "../scripts/scatter_plot_real_data.py"
