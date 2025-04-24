import polars as pl
import pandas as pd
import altair as alt
import numpy as np
import pyarrow
import os
import json
import sys
import collections
import scipy.stats

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    def create_scatter_plot(data):
        scatterplot = alt.Chart(data).mark_point().transform_calculate(jitter="random()").encode(
            x='Actual:Q',
            y='Predicted:Q',
            color=alt.value("black"),
            tooltip=alt.Tooltip(['Actual', 'Predicted'])
        )

        # add a diagonal line
        line = pd.DataFrame({
        'Actual': [0.0, 1.0],
        'Predicted':  [0.0, 1.0],
        })

        line_plot = alt.Chart(line).mark_line(color='red').encode(
        x= 'Actual',
        y= 'Predicted',
        )

        plot = (scatterplot + line_plot).properties(
            width=200,
            height=200
        )
        return plot

    ## ORTHANQ ##
    print("----- ORTHANQ -----")

    #orthanq results
    orthanq_input = snakemake.input.orthanq_prediction

    #fractions in the truth
    truth = {"896": 0.226, "HXB2": 0.10, "JRCSF": 0.296, "NL43": 0.269, "YU2": 0.109}

    #read table
    results = pd.read_csv(orthanq_input)

    ##find the records that have the same density
    best_odds = [1, 1.0, 1.00]
    best_results = results[results['odds'].isin(best_odds)]
    print("---best_results---")
    print(best_results)

    #remove density and odds columns
    best_results = best_results.drop(columns=['density', 'odds'])
    
    predicted_fractions = best_results.iloc[0].tolist()

    #checked by eye that orthanq column headers with haplotypes and truth align each other
    print("truth: ",list(truth.values()))
    print("prediction: ",predicted_fractions)

    orthanq_data = pd.DataFrame({
        'Actual': list(truth.values()),
        'Predicted': predicted_fractions,
    })
  
    orthanq_plot = create_scatter_plot(orthanq_data)

    #export to svg and html
    orthanq_plot.save(snakemake.output.orthanq_svg)
    orthanq_plot.save(snakemake.output.orthanq_html)

    #calculate spearman correlation
    spearman=scipy.stats.spearmanr(list(truth.values()), predicted_fractions)   # Spearman's rho
    print("spearman",spearman)

    ## KALLISTO ##
    print("----- KALLISTO -----")

    #kallisto results
    kallisto_input = snakemake.input.kallisto_prediction

    #read table
    results = pd.read_csv(kallisto_input, sep="\t")
    print("results: ",results)

    #read table
    predicted_fractions_list = dict(zip(results["target_id"], results["tpm"] / results["tpm"].sum()))

    # Sort fractions_dict and truths dict so that the values align each other
    sorted_fractions_dict = dict(sorted(predicted_fractions_list.items()))
    sorted_truth = dict(sorted(truth.items()))
    print("sorted_fractions_dict: ", sorted_fractions_dict)
    print("sorted_truth: ", sorted_truth)

    #create data
    kalllisto_data = pd.DataFrame({
        'Actual': list(sorted_truth.values()),
        'Predicted': list(sorted_fractions_dict.values())
    })

    kallisto_plot = create_scatter_plot(kalllisto_data)

    #export to svg and html
    kallisto_plot.save(snakemake.output.kallisto_svg)
    kallisto_plot.save(snakemake.output.kallisto_html)

    #calculate spearman correlation
    spearman=scipy.stats.spearmanr(list(sorted_truth.values()), list(sorted_fractions_dict.values()))   # Spearman's rho
    print("spearman",spearman)