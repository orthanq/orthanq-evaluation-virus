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

    #the following are the predicted lineages by orthanq. Here we want to create a scatterplot for orthanq. 
    #For Kallisto, we would like to select these five lineages in the scatterplot and also add the remaining fractions to it.
    #the rename map ise used to map the accession ids of those lineages to lineage names used in the original paper.
    rename_map = {
    "AF324493.2": "NL43",
    "K03455.1": "HXB2",
    "M38429.1": "JRCSF",
    "M93258.1": "YU2",
    "U39362.2": "896"}

    #fractions in the truth
    truth = {"896": 0.226, "HXB2": 0.10, "JRCSF": 0.296, "NL43": 0.269, "YU2": 0.109}


    def create_scatter_plot(data):
        scatterplot = alt.Chart(data).mark_point().transform_calculate(jitter="random()").encode(
            x='Actual:Q',
            y='Predicted:Q',
            color=alt.value("black"),
            tooltip=alt.Tooltip(['Actual', 'Predicted'])
        )

        # add a diagonal line
        line = pd.DataFrame({
        'Actual': [0.0, 0.4],
        'Predicted':  [0.0, 0.4],
        })

        line_plot = alt.Chart(line).mark_line(color='red').encode(
        x= alt.X('Actual',scale=alt.Scale(domain=[0.0, 0.4])),
        y= alt.Y('Predicted', scale=alt.Scale(domain=[0.0, 0.4])),
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

    # read table
    results = pd.read_csv(orthanq_input)

    # find the records that have the same density
    best_odds = [1, 1.0, 1.00]
    best_results = results[results['odds'].isin(best_odds)]
    print("---best_results---")
    print(best_results)

    # remove density and odds columns
    best_results = best_results.drop(columns=['density', 'odds'])
    
    # select only the relevant target_id columns
    target_cols = list(rename_map.keys())

    # get the first row's values for those columns
    row_values = best_results.loc[0, target_cols]

    # rename the keys and convert to dictionary
    renamed_fractions = {
        rename_map[col]: row_values[col] for col in target_cols
    }

    # sort the dictionary; no need to sort the truth as it is already sorted.
    sorted_renamed_fractions = dict(sorted(renamed_fractions.items()))

    # print the dicts for easy debugging
    print("renamed_fractions", sorted_renamed_fractions)
    print("truth: ",truth)

    orthanq_data = pd.DataFrame({
        'Actual': list(truth.values()),
        'Predicted': list(sorted_renamed_fractions.values()),
    })
  
    orthanq_plot = create_scatter_plot(orthanq_data)

    # export to svg and html
    orthanq_plot.save(snakemake.output.orthanq_svg)
    orthanq_plot.save(snakemake.output.orthanq_html)

    # calculate spearman correlation
    spearman=scipy.stats.spearmanr(list(truth.values()), list(sorted_renamed_fractions.values()))   # Spearman's rho
    print("spearman",spearman)

    ## KALLISTO ##
    print("----- KALLISTO -----")

    # kallisto results
    kallisto_input = snakemake.input.kallisto_prediction
 
    # define the list of accessions to be used according to: AF324493.2: NL4-3, K03455.1: HXB2, M38429.1: JRCSF, M93258.1: YU-2, U39362.2: 896
    target_filter_list = ["AF324493.2", "K03455.1", "M38429.1", "M93258.1", "U39362.2"]  
    
    # read table
    results = pd.read_csv(kallisto_input, sep="\t")
    print("results: ",results)

    #calculate fraction estimates of each accession
    predicted_fractions_dict = dict(zip(
        results["target_id"],
        results["tpm"] / results["tpm"].sum()
    ))

    print("predicted_fractions_dict", predicted_fractions_dict)

    # entries that are in the list
    matched_results = {
        k: v for k, v in predicted_fractions_dict.items() if k in target_filter_list
    }
    print("matched_results", matched_results)

    # rename keys using the mapping
    renamed_matched_predictions = {
            rename_map.get(k, k): v for k, v in matched_results.items()}
    print("renamed_matched_predictions", renamed_matched_predictions)

    # Sort fractions_dict and truths dict so that the values align each other
    sorted_fractions_dict = dict(sorted(renamed_matched_predictions.items()))
    sorted_truth = dict(sorted(truth.items()))

    print("sorted_fractions_dict: ", sorted_fractions_dict)
    print("sorted_truth: ", sorted_truth)

    # add zero points to the truth just to add them to the plot for transparency

    # first find entries that are **not** in the initial list (unmatched)
    unmatched_results = {
        k: v for k, v in predicted_fractions_dict.items() if k not in target_filter_list
    }
    
    # then, find nonzero truths and add as many 0.0 fractions to yield a final list of truth fractions (to show in the plot)
    nonzero_truth = list(sorted_truth.values())
    zero_fractions = [0.0] * len(unmatched_results)
    all_truth = nonzero_truth + zero_fractions
    print("all_truth", all_truth)

    #add remaining fractions to yield final predicted fractions
    matched_fractions = list(sorted_fractions_dict.values())
    all_predicted_fractions = matched_fractions + list(unmatched_results.values())

    #create data
    kalllisto_data = pd.DataFrame({
        'Actual': all_truth,
        'Predicted': all_predicted_fractions
    })

    kallisto_plot = create_scatter_plot(kalllisto_data)

    #export to svg and html
    kallisto_plot.save(snakemake.output.kallisto_svg)
    kallisto_plot.save(snakemake.output.kallisto_html)

    #calculate spearman correlation
    spearman=scipy.stats.spearmanr(all_truth, all_predicted_fractions)   # Spearman's rho
    print("spearman",spearman)