#we throw all data points together with the fractions and each sample is read alphabetically (e.g. SimulatedSample1, SimulatedSample2 etc.) for both predictions and simulations
#here it is assumed that two scatterplots for two coverages are generated, 100x and 1000x.

# import pandas as pd
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

    #simulation input
    simulation = snakemake.input.simulation

    #1-) collect all simulated fractions
    sample_composition = {}

    for file in simulation:
        #retrieve sample name
        sample_name = os.path.basename(file).split(".")

        #read table
        results = pd.read_csv(file)

        #select fraction and lineage column
        fractions = results['fraction'].tolist()
        lineages = results['lineage'].tolist()

        #record sample composition
        fractions_per_lineage = {lineage: fraction for fraction, lineage in zip(fractions, lineages)}
        sample_composition[sample_name[0]] = fractions_per_lineage
    print("sample_composition", sample_composition)

    # #2-) collect all predicted fractions of orthanq, separately by coverages 
    predicted_fractions_100x=[]
    predicted_fractions_1000x=[] 

    for file in orthanq_input:
        #get sample and coverage
        splitted = os.path.basename(file).split(".")[0].split("-") #first split is to get rid of the extension
        sample = splitted[0]
        cov = splitted[1]
        print("sample:", sample)
        print("cov", cov)
        
        #set correct dictionary in case of differing coverages
        if cov == "100x": 
            predicted_fractions_final = predicted_fractions_100x
        elif cov == "1000x":
            predicted_fractions_final = predicted_fractions_1000x
        
        #read table
        results = pd.read_csv(file)

        ##find the records that have the same density
        best_odds = [1, 1.0, 1.00]
        best_results = results[results['odds'].isin(best_odds)]
        print("---best_results---")
        print(best_results)

        #remove density and odds columns
        best_results = best_results.drop(columns=['density', 'odds'])

        #names of predicted haplotypes in the results
        predicted_haplotypes = best_results.columns.tolist()
        print("predicted_haplotypes", predicted_haplotypes)

        #find true haplotype fractions for this sample
        simulated_fractions = sample_composition[sample]
        print("simulated_fractions: ",simulated_fractions)

        rmse_values = {}
        for (i,row) in best_results.iterrows():

            fractions_real = []
            fractions_predicted = []

            #loop over each lineage,fraction in the simulation and check if the prediction contains the lineage, if yes, take the RMSE
            #if not, record it as 0.0 in the prediction.
            for (l,f) in simulated_fractions.items():
                fractions_real.append(simulated_fractions[l])
                if l in predicted_haplotypes:
                    fractions_predicted.append(row[l])
                else:
                    fractions_predicted.append(0.0)
            
            print("fractions_real: ", fractions_real)
            print("fractions_predicted: ", fractions_predicted)
            
            #calculate rmse of this solution
            rmse = np.sqrt(np.mean((np.array(fractions_real) - np.array(fractions_predicted)) ** 2))
            rmse_values[rmse] = [fractions_real, fractions_predicted]

        print("rmse values:", rmse_values)
        #find the smallest key
        smallest_key = min(rmse_values.keys())        
        print("smallest_key", smallest_key)

        #extend the final prediction list by the values belonging to this rmse value, first index belongs to prediction
        predicted_fractions_final.extend(rmse_values[smallest_key][1])

    print("predicted_fractions_100x",predicted_fractions_100x)
    print("predicted_fractions_1000x",predicted_fractions_1000x)

    #svg paths
    sorted_output_svg = sorted(snakemake.output.svg, key=lambda x: int(x.split('_')[-1].replace('x', '').replace('.svg', '')))
    plot_100x_svg = sorted_output_svg[0]
    plot_1000x_svg = sorted_output_svg[1]

    #svg paths
    sorted_output_html = sorted(snakemake.output.html, key=lambda x: int(x.split('_')[-1].replace('x', '').replace('.html', '')))
    plot_100x_html = sorted_output_html[0]
    plot_1000x_html = sorted_output_html[1]

    for cov in ["100x", "1000x"]:
        print("cov: ", cov)
        if cov == "100x":
            predicted_fractions = predicted_fractions_100x
        elif cov == "1000x":
            predicted_fractions = predicted_fractions_1000x

        #find x and y values for the plots
        x_values = [value for sample in sample_composition.values() for value in sample.values()]
        y_values = predicted_fractions
        print("x_values: ", x_values)
        print("y_values", y_values)
        
        data = pd.DataFrame({
            'Actual': x_values,
            'Predicted': y_values,
        })
        print(data)

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
        if cov == "100x":
            plot.save(plot_100x_svg)
            plot.save(plot_100x_html)
        elif cov == "1000x":
            plot.save(plot_1000x_svg)
            plot.save(plot_1000x_html)

        #calculate spearman correlation
        spearman=scipy.stats.spearmanr(x_values, y_values)   # Spearman's rho
        print("spearman",spearman)
        
        # # df.write_csv(snakemake.output.foo)
