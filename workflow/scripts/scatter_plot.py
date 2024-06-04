#scatterplot for all points
#in this plot we throw all data points together with the fractions 
#and each sample is read alphabetically (e.g. SimulatedSample1, SimulatedSample2 etc.) for both predictions and simulations
#in the end, each datapoint for the corresponding lineage is paired to each other.

# import pandas as pd
import polars as pl
import pandas as pd
import altair as alt
import pyarrow
import os
import json
import sys
import collections
import scipy.stats

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    def flatten(xss):
        return [x for xs in xss for x in xs]

    #orthanq results
    orthanq_input = snakemake.input.orthanq_prediction

    #simulation input
    simulation = snakemake.input.simulation


    #1-) collect all simulated fractions

    #initialize a list to collect simulated fractions: lineage_name: fraction
    simulated_fractions={}

    #keep record of number of fractions for each sample: sample: lineage_names
    lineages_simulations={}

    for file in simulation:
        #retrieve sample name
        sample_name = os.path.basename(file).split(".")
        print(sample_name[0])

        #read table
        results = pl.read_csv(file)

        #select fraction and lineage column
        fractions = results.select('fraction').rows()
        lineages = results.select('lineage').rows()

        #loop over fractions and lineages
        for f,l in zip(fractions, lineages):
            print(f[0], l[0])
            lineage = l[0]
            fraction = f[0]
            if lineage in simulated_fractions:
                fracs_for_key = simulated_fractions[lineage]
                fracs_for_key.append(fraction)
                simulated_fractions[lineage] = fracs_for_key
            else:
                simulated_fractions[lineage] = [fraction]
        
        #flatten the list of lists
        flattened = flatten(lineages)

        #lastly add number of fractions 
        lineages_simulations[sample_name[0]] = flattened

    print(lineages_simulations)

    #2-) collect all predicted fractions of orthanq

    #initialize a list to collect predicted fractions
    predicted_fractions={}

    for file in orthanq_input:

        #retrieve sample name
        sample_name = os.path.basename(file).split(".")
        print(sample_name[0])

        #read table
        results = pl.read_csv(file)

        #first row
        best_results = results[0]
        list_of_first_row = best_results.rows(named=True)[0]

        #remove density and odds columns
        list_of_first_row.pop('density', None)
        list_of_first_row.pop('odds', None)

        for k,v in list_of_first_row.items():
            print(k,v)
            if k in predicted_fractions:
                fracs_for_key = predicted_fractions[k]
                fracs_for_key.append(v)
                predicted_fractions[k] = fracs_for_key
            else:
                predicted_fractions[k] = [v]

        #*some lineages are not present in orthanq predictions
        #append 0.0 fractions as long as the number of fractions in simulated sample
        sample_lineages = lineages_simulations[sample_name[0]]
        
        for lin in sample_lineages:
            lineages_in_results = list_of_first_row.keys()
            print("lin:", lin)
            print("lineages_in_results:",lineages_in_results)
            if not lin in lineages_in_results:
                if lin in predicted_fractions:
                    fracs_for_lin = predicted_fractions[lin]
                    fracs_for_lin.append(0.0)
                    predicted_fractions[lin] = fracs_for_lin
                else:
                    predicted_fractions[lin] = [0.0]

    print("simulated_fractions", simulated_fractions)  
    print(len(predicted_fractions))
    print("predicted_fractions",predicted_fractions)

    #now sort both dictionaries
    ordered_simulated_fractions = dict(sorted(simulated_fractions.items()))
    ordered_predicted_fractions = dict(sorted(predicted_fractions.items()))

    print("ordered_simulated_fractions", ordered_simulated_fractions)  
    print("ordered_predicted_fractions",ordered_predicted_fractions)

    #*some lineages are present in simulated fractions
    #remove them (ask johannes)
    ordered_predicted_fractions_copy = ordered_predicted_fractions.copy() #an explicit copy
    for k_o, v_o in ordered_predicted_fractions.items():
        if k_o not in ordered_simulated_fractions:
            ordered_predicted_fractions_copy.pop(k_o, None)
    print("ordered_predicted_fractions_copy:", ordered_predicted_fractions_copy)

    #convert to dataframe
    data = {"predicted_fractions": flatten(ordered_predicted_fractions_copy.values()), "simulated_fractions": flatten(ordered_simulated_fractions.values())}
    print(data)
    df = pl.DataFrame(data)
    print(df)

    #throw all fractions
    
    # scatter plot
    scatterplot = alt.Chart(df).mark_point().transform_calculate(jitter="random()").encode(
        x='simulated_fractions:Q',
        y='predicted_fractions:Q',
        # xOffset='jitter:Q',
        color=alt.value("black"),
    )

    #add a diagonal line
    line = pd.DataFrame({
    'simulated_fractions': [0.0, 1.0],
    'predicted_fractions':  [0.0, 1.0],
    })

    line_plot = alt.Chart(line).mark_line(color='red').encode(
    x= 'simulated_fractions',
    y= 'predicted_fractions',
    )

    plot = scatterplot + line_plot

    #export 
    plot.save(snakemake.output.plot)

    #calculate spearman correlation
    spearman=scipy.stats.spearmanr(flatten(ordered_predicted_fractions_copy.values()), flatten(ordered_simulated_fractions.values()))   # Spearman's rho
    print(spearman)
    
    # df.write_csv(snakemake.output.foo)
