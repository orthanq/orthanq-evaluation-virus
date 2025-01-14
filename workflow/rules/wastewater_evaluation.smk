rule evaluate_wastewater:
    input:
        orthanq = expand("results/orthanq/calls/{sample}/{sample}.csv", sample=samples["sra"]),
        metadata="results/project_metadata/sra_collection_dates.tsv"
    output:
        validation_B117="results/evaluation/wastewater_validation_B117.svg",
        validation_B1526="results/evaluation/wastewater_validation_B1526.svg"
    log:
        "logs/evaluation/wastewater_validation.log"
    conda:
        "../envs/altair.yaml"
    benchmark:
        "benchmarks/evaluation/wastewater_validation.tsv" 
    script:
        "../scripts/wastewater_validation.py"