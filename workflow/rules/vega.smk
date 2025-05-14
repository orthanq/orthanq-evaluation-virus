rule vg2svg_orthanq:
    input:
        solutions="results/orthanq/calls/{sample}/viral_solutions.json",
        lp_solution="results/orthanq/calls/{sample}/lp_solution.json",
        final_solution="results/orthanq/calls/{sample}/final_solution.json"
    output:
        solutions=report("results/orthanq/calls/{sample}/viral_solutions.html", 
        caption="../report/orthanq_plots.rst",
        category="Orthanq detailed solutions", subcategory="{sample}",labels={
            "figure": "3-field solutions"
        }),
        lp_solution="results/orthanq/calls/{sample}/lp_solution.html",
        final_solution=report("results/orthanq/calls/{sample}/final_solution.html",
        caption="../report/orthanq_plots.rst",
        category="Orthanq detailed solutions", subcategory="{sample}", labels={
            "figure": "final solution"
        })        
    log:
        "logs/vg2svg/orthanq/{sample}.log",
    conda:
        "../envs/vega.yaml"
    shell:
        "vl-convert vl2html --input {input.solutions} --output {output.solutions} 2> {log} && "
        "vl-convert vl2html --input {input.lp_solution} --output {output.lp_solution} 2>> {log} && "
        "vl-convert vl2html --input {input.final_solution} --output {output.final_solution} 2>> {log}"