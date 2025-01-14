#validation with the most abundant lineages? Analysis in terms of call rate and accuracy.

import polars as pl
import pandas as pd
import os
import sys
import functools as ft
import re

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #natural sorting (the sample names contains both string and number hence sorted() cannot sort them properly)
    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]

    def extract_ground_truth(simulation):
        #collect the abundant lineage in simulations
        abundant_lineage_per_sample={}

        for file in simulation:
            sample_name = os.path.basename(file).split(".")
            print(sample_name)
            #read table
            results = pd.read_csv(file)

            #find the largest fraction
            largest_fraction=0.0
            largest_fraction_lineage=""
            for (i,row) in results.iterrows():
                print("row", row)
                if row['fraction'] > largest_fraction:
                    largest_fraction=row['fraction']
                    largest_fraction_lineage=row['lineage']
            abundant_lineage_per_sample[sample_name[0]] = largest_fraction_lineage
            print("sample_name[0]",sample_name[0])
        return abundant_lineage_per_sample

    def extract_orthanq(list_of_input):
        print("ORTHANQ")
        #collect the abundandant lineage in predictions
        #initialize orthanq predictions
        predictions_100x={}
        predictions_1000x={}

        for file in list_of_input:
            splitted = os.path.basename(file).split(".")[0].split("-") #first split is to get rid of the extension
            sample = splitted[0]
            cov = splitted[1]
            print("sample:", sample)
            print("cov", cov)

            #set correct dictionary in case of differing coverages
            if cov == "100x": 
                predictions = predictions_100x
            elif cov == "1000x":
                predictions = predictions_1000x
        
            #read table
            results = pd.read_csv(file)

            ##find the records that have the same density
            best_odds = [1, 1.0, 1.00]
            best_results = results[results['odds'].isin(best_odds)]
            print("---best_results---")
        
            #remove density and odds columns
            best_results_filtered = best_results.drop(columns=['density', 'odds'])

            #names of haplotypes in the results
            lineages = best_results_filtered.columns.tolist()
            print("lineages", lineages)

            abundant_lineages = []

            for (i,row) in best_results_filtered.iterrows():
                largest_fraction=0.0
                largest_fraction_lineage=""
                for lin in lineages:
                    if row[lin] > largest_fraction:
                        largest_fraction=row[lin]
                        largest_fraction_lineage=lin
                abundant_lineages.append(largest_fraction_lineage)
            predictions[sample]=abundant_lineages
        return (predictions_100x, predictions_1000x)
    

    def extract_pangolin(list_of_input):
        print("PANGOLIN")
        predictions_100x = {}
        predictions_1000x = {}
        for file in list_of_input:
            splitted = os.path.basename(file).split(".")[0].split("-") #first split is to get rid of the extension
            sample = splitted[0]
            cov = splitted[1]
            print("sample:", sample)
            print("cov", cov)
        
            #set correct dictionary in case of differing coverages
            if cov == "100x": 
                predictions = predictions_100x
            elif cov == "1000x":
                predictions = predictions_1000x
        
            #extract the prediction in 'lineage' column; pangolin only outputs one row
            print("ch1")
            df=pd.read_csv(file)
            predictions[sample] = df.loc[0,"lineage"]
        return (predictions_100x, predictions_1000x)
    
    def extract_kallisto(list_of_input):
        print("KALLISTO")
        predictions_100x = {}
        predictions_1000x = {}
        for file in list_of_input:
            splitted = os.path.basename(os.path.dirname(file)).removeprefix("quant_results_").split("-") #first split is to get rid of the extension
            sample = splitted[0]
            cov = splitted[1]
            print("sample:", sample)
            print("cov", cov)
        
            #set correct dictionary in case of differing coverages
            if cov == "100x": 
                predictions = predictions_100x
            elif cov == "1000x":
                predictions = predictions_1000x

            df=pd.read_csv(file, sep="\t")

            #find the row for column 'target_id' (lineage) that has the largest value for column 'tpm'
            sorted_df = df.sort_values(by='tpm', ascending=False).reset_index(drop=True)
            print(sorted_df)
            lineage = sorted_df.loc[0,'target_id']
            print(lineage)
            predictions[sample] = lineage
        return (predictions_100x, predictions_1000x)
    
    def extract_nextclade(list_of_input):
        print("NEXTCLADE")
        predictions_100x = {}
        predictions_1000x = {}        
        for file in list_of_input:
            splitted = os.path.basename(os.path.dirname(file)).split("-") #first split is to get rid of the extension
            sample = splitted[0]
            cov = splitted[1]
            print("sample:", sample)
            print("cov", cov)
        
            #set correct dictionary in case of differing coverages
            if cov == "100x": 
                predictions = predictions_100x
            elif cov == "1000x":
                predictions = predictions_1000x

            #read the file into a df
            df=pd.read_csv(file, sep="\t")
            #extract the prediction in 'clade' column 
            predictions[sample] = df.loc[0,"clade"]
        return (predictions_100x, predictions_1000x)
    
    def validate_nextclade(clade_to_lineage, predictions, ground_truth):
        if predictions.keys() != ground_truth.keys():
            raise Exception("call rate cannot be 1 because of different lengths")
        clade_to_lineage=pd.read_csv(clade_to_lineage, header=None, sep="\t")

        #first rename the columns
        clade_to_lineage.columns = ['clade', 'pango']

        sorted_keys = sorted(predictions.keys(), key=lambda x: int(x.replace("SimulatedSample", "")))
        call_rate=1.0
        n_correct=0
        #calculate accuracy
        for key in sorted_keys:
            truth = ground_truth[key]
            pred = predictions[key]

            #convert pred from clade to pango
            # pango_lineages = clade_to_lineage.loc[clade_to_lineage['clade'] == pred, 'pango'].tolist()
            pango_lineages = clade_to_lineage.loc[clade_to_lineage['clade'].str.contains(pred, case=False, na=False), 'pango'].tolist()
            print("nextclade prediction:", pred)
            print("pango_lineages:",pango_lineages)

            if truth in pango_lineages:
                n_correct+=1
        accuracy=n_correct/len(predictions)*100 #it can also be the length of simulated samples, but since they're the same, it doesn't matter which of them is used here
        print(accuracy)

        #convert to dataframe and write
        validation = pd.DataFrame({"n_samples": [len(predictions)], "pangolin_call_rate": [call_rate], "pangolin_accuracy": [accuracy]})

        return validation
    
    def validate_tool(tool_name, predictions, ground_truth):
        if predictions.keys() != ground_truth.keys():
            raise Exception("call rate cannot be 1 because of different lengths")
        sorted_keys = sorted(predictions.keys(), key=lambda x: int(x.replace("SimulatedSample", "")))
        call_rate=1.0
        n_correct=0
        #calculate accuracy
        for key in sorted_keys:
            truth = ground_truth[key]
            preds = predictions[key]
            if truth in preds:
                n_correct+=1
        accuracy=n_correct/len(predictions)*100
        print("accuracy", accuracy)
        #create column names with tools
        call_rate_column=tool_name+"_"+"call_rate"
        accuracy_column=tool_name+"_"+"accuracy"

        #convert to dataframe and write
        tool_validation = pd.DataFrame({"n_samples": [len(predictions)], call_rate_column: [call_rate], accuracy_column: [accuracy]})

        return tool_validation

    print("---extracting ground truth---")
    simulation = list(snakemake.input.simulation)
    ground_truth = extract_ground_truth(simulation)
    print(ground_truth)

    #validate orthanq
    orthanq_input = list(snakemake.input.orthanq)
    (orthanq_predictions_100x,orthanq_predictions_1000x)=extract_orthanq(orthanq_input)
    orthanq_validation_100x = validate_tool("orthanq", orthanq_predictions_100x, ground_truth)
    orthanq_validation_1000x = validate_tool("orthanq", orthanq_predictions_1000x, ground_truth)

    #validate pangolin
    pangolin_input = list(snakemake.input.pangolin)
    (pangolin_predictions_100x,pangolin_predictions_1000x)=extract_pangolin(pangolin_input)
    pangolin_validation_100x=validate_tool("pangolin", pangolin_predictions_100x, ground_truth)        
    pangolin_validation_1000x=validate_tool("pangolin", pangolin_predictions_1000x, ground_truth)        
    
    #validate kallisto
    kallisto_input = list(snakemake.input.kallisto)
    (kallisto_predictions_100x, kallisto_predictions_1000x)=extract_kallisto(kallisto_input)
    kallisto_validation_100x=validate_tool("kallisto", kallisto_predictions_100x, ground_truth)
    kallisto_validation_1000x=validate_tool("kallisto", kallisto_predictions_1000x, ground_truth)

    #nextclade
    nextclade_input = list(snakemake.input.nextclade)
    (nextclade_predictions_100x, nextclade_predictions_1000x)=extract_nextclade(nextclade_input)
    nextclade_validation_100x=validate_nextclade(snakemake.input.clade_to_lineage, nextclade_predictions_100x, ground_truth)
    nextclade_validation_1000x=validate_nextclade(snakemake.input.clade_to_lineage, nextclade_predictions_1000x, ground_truth)

    sorted_output = sorted(snakemake.output.validation, key=lambda x: int(x.split('_')[-1].replace('x', '').replace('.tsv', '')))
    [path_100x, path_1000x] = sorted_output
    print("path_100x", path_100x)
    print("path_1000x", path_1000x)

    for cov in ["100x", "1000x"]:
        if cov == "100x":
            dfs = [orthanq_validation_100x, pangolin_validation_100x, kallisto_validation_100x]
            dfs_final = ft.reduce(lambda left, right: pd.merge(left, right, on='n_samples', how='left'), dfs)
            dfs_final.to_csv(path_100x, index=False)
        elif cov == "1000x":
            dfs = [orthanq_validation_1000x, pangolin_validation_1000x, kallisto_validation_1000x]
            dfs_final = ft.reduce(lambda left, right: pd.merge(left, right, on='n_samples', how='left'), dfs)
            dfs_final.to_csv(path_1000x, index=False)
