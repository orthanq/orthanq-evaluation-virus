rule scatter_plot:
    input:
        orthanq_prediction=expand("results/orthanq/calls/SimulatedSample{num}/SimulatedSample{num}.tsv", num=num_list),
        simulation=expand("results/simulation_input/SimulatedSample{num}.csv", num=num_list),
    output:
        plot="results/evaluation/scatter_plot.svg",
        # foo="foo.csv"
    log:
        "logs/scatter_plot/scatter_plot.log"
    conda:
        "../envs/altair.yaml"
    benchmark:
        "benchmarks/scatter_plot/scatter_plot.tsv" 
    script:
        "../scripts/scatter_plot.py"