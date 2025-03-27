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

    #orthanq results
    orthanq_input = snakemake.input.orthanq_prediction

    #fractions in the truth
    truth = {"896": 0.226, "HXB2": 0.10, "JR-CSF": 0.296, "NL4-3": 0.269, "YU2": 0.109}

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

    print("truth: ",list(truth.values()))
    print("prediction: ",predicted_fractions)
    data = pd.DataFrame({
        'Actual': list(truth.values()),
        'Predicted': predicted_fractions,
    })

    # scatter plot
    scatterplot = alt.Chart(data).mark_point().transform_calculate(jitter="random()").encode(
        x='Actual:Q',
        y='Predicted:Q',
        # xOffset='jitter:Q',
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

    plot = scatterplot + line_plot    

    #export to svg and html
    plot.save(snakemake.output.svg)
    plot.save(snakemake.output.html)

    #calculate spearman correlation
    spearman=scipy.stats.spearmanr(list(truth.values()), predicted_fractions)   # Spearman's rho
    print("spearman",spearman)