#validation with the most abundant lineages? Analysis in terms of call rate and accuracy.

import polars as pl
import pandas as pd
import os
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #orthanq results
    orthanq_input = snakemake.input.orthanq_prediction

    #simulation input
    simulation = snakemake.input.simulation

    #initialize orthanq_predictions
    orthanq_predictions=[]

    #collect the abundandant lineage in predictions
    for file in orthanq_input:
        #read table
        results = pl.read_csv(file)

        #get best results
        best_results = results[0]
        list_of_first_row = best_results.rows(named=True)[0]
        print(list_of_first_row)

        #remove density and odds columns
        list_of_first_row.pop('density', None)
        list_of_first_row.pop('odds', None)

        #find the largest abundant lineage
        max_value=max(list_of_first_row, key=list_of_first_row.get)
        print("max_value:", max_value)
        orthanq_predictions.append(max_value)

    print(orthanq_predictions)

    #collect the abundant lineage in simulations
    simulation_largest_lineages=[]

    for file in simulation:
        #read table
        results = pl.read_csv(file)

        #find the largest fraction
        largest_fraction=0.0
        largest_fraction_lineage=""
        for row in results.rows(named=True):
            if row['fraction'] > largest_fraction:
                largest_fraction=row['fraction']
                largest_fraction_lineage=row['lineage']
        
        simulation_largest_lineages.append(largest_fraction_lineage)
    print(simulation_largest_lineages)

    #we see that orthanq always has a call, because we never saw a case without in viral lineage quantification
    #but if that's not the case, raise an exception
    if len(orthanq_predictions)!=len(simulation_largest_lineages):
        raise Exception("call rate cannot be 1 because of different lengths")
    call_rate=1.0
    n_correct=0
    #calculate call rate and accuracy
    for pred,sim in zip(orthanq_predictions,simulation_largest_lineages):
        if pred==sim:
            n_correct+=1

    accuracy=n_correct/len(orthanq_predictions)*100 #it can also be the lenght of simulated samples, but since they're the same, it doesn't matter which of them is used here
    print(accuracy)
    #convert to dataframe and write
    data = {"n_samples": len(orthanq_predictions),"call_rate": call_rate, "accuracy": accuracy}
    df = pl.DataFrame(data)
    print(df)
    df.write_csv(snakemake.output.validation)