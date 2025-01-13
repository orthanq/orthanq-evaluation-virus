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
    
    def validate_tool(tool_name, predictions, ground_truth):

        if predictions.keys() != ground_truth.keys():
            raise Exception("call rate cannot be 1 because of different lengths")
        call_rate=1.0
        sorted_keys = sorted(predictions.keys(), key=lambda x: int(x.replace("SimulatedSample", "")))
        n_correct=0
        #calculate call rate and accuracy
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

    orthanq_input = list(snakemake.input.orthanq)
    (orthanq_predictions_100x,orthanq_predictions_1000x)=extract_orthanq(orthanq_input)
    print("orthanq_predictions_100x: ",orthanq_predictions_100x)
    print("orthanq_predictions_1000x: ",orthanq_predictions_1000x)

    #validate orthanq
    orthanq_validation_100x = validate_tool("orthanq", orthanq_predictions_100x, ground_truth)
    orthanq_validation_1000x = validate_tool("orthanq", orthanq_predictions_1000x, ground_truth)

    print("orthanq_validation_100x",orthanq_validation_100x)
    print("orthanq_validation_1000x",orthanq_validation_1000x)

    #############################################
    
    #############################################
    #pangolin####################################
    #############################################

    # #assign list of pangolin results to a variable
    # pangolin_input = list(snakemake.input.pangolin)

    # def extract_pangolin(list_of_input):
    #     #initialize a list to collect predictions
    #     predictions=[]

    #     #loop over the files

    #     #sort the list of files
    #     list_of_input.sort(key=natural_keys)
    #     print("sorted list: ", list_of_input)

    #     for file in list_of_input:

    #         #predictions
    #         #read the file into a df
    #         df=pl.read_csv(file)

    #         #pangolin only output one row
    #         df=df[0]

    #         #extract the prediction in 'lineage' column 
    #         predictions.append(df[0,"lineage"])

    #     return predictions
    # #extract pangolin predictions
    # pangolin_predictions=extract_pangolin(pangolin_input)
    # print("pangolin predictions: ", pangolin_predictions)

    # #validate pangolin (convert pango lineages to nextstrain, 
    # #give a score to pango prediction if the predicted lineage has the correspondence that matches the ground truth)
    # clade_to_lineage=snakemake.input.clade_to_lineage

    # def validate_pangolin(pangolin_predictions, clade_to_lineage, sorted_ground_truth, coverage_number):
    #     #duplicate ground truth because there there are two coverage (100x and 1000x) samples
    #     duplicated_ground_truth = [x for x in sorted_ground_truth for _ in range(coverage_number)]

    #     #convert pango lineages to nextstrain clades
    #     clade_to_lineage=pl.read_csv(clade_to_lineage, has_header=False, separator="\t")
    #     print(clade_to_lineage)

    #     #first rename the columns
    #     clade_to_lineage=clade_to_lineage.rename({"column_1": "clade", "column_2": "pango"})

    #     #loop over pangolin predictions to calculate accuracy
    #     n_correct=0
    #     for idx,pango_lineage in enumerate(pangolin_predictions):
    #         #find the rows with pango_lineage in clade_to_lineage
    #         filtered=clade_to_lineage.filter(pl.col("pango") == pango_lineage)
    #         print(filtered)
            
    #         #create a list for all corresponding clades (can be more than one)
    #         corresp_clades=list(filtered['clade'])
    #         print(corresp_clades)

    #         #calculate call rate and accuracy
    #         truth_at_index=ground_truth[idx]
    #         print("truth_at_index: ", truth_at_index)
    #         print("corresp_clades: ", corresp_clades)
    #         if truth_at_index in corresp_clades:
    #             print("test")
    #             n_correct+=1
    #     call_rate=1.0
    #     accuracy=n_correct/len(pangolin_predictions)*100 #it can also be the length of simulated samples, but since they're the same, it doesn't matter which of them is used here
    #     print(accuracy)

    #     #convert to dataframe and write
    #     pangolin_validation = pl.DataFrame({"n_samples": len(pangolin_predictions), "pangolin_call_rate": call_rate, "pangolin_accuracy": accuracy})

    #     return pangolin_validation
    
    # pangolin_validation=validate_pangolin(pangolin_predictions, clade_to_lineage, ground_truth, snakemake.params.coverage_number)        
    # # pangolin_validation=validate_tool("pangolin", pangolin_predictions, ground_truth)
    # print("pangolin validation: ", pangolin_validation)

    # #############################################
    
    # #############################################
    # #nextclade###################################
    # #############################################    

    # #assign list of pangolin results to a variable
    # nextclade_input = list(snakemake.input.nextclade)
    # print(nextclade_input)
    # def extract_nextclade(list_of_input):
    #     #initialize a list to collect predictions
    #     predictions=[]

    #     #loop over the files
        
    #     #sort the list of files
    #     list_of_input.sort(key=natural_keys)
    #     print("sorted list: ", list_of_input)

    #     for file in list_of_input:

    #         #predictions
    #         #read the file into a df
    #         df=pl.read_csv(file, separator="\t")

    #         #nextclade only output one row
    #         df=df[0]

    #         #extract the prediction in 'clade' column 
    #         predictions.append(df[0,"clade"])

    #     return predictions
    # #extract nextclade predictions
    # nextclade_predictions=extract_nextclade(nextclade_input)
    # print("nextclade predictions: ", nextclade_predictions)

    # #validate pangolin
    # nextclade_validation=validate_tool("nextclade", nextclade_predictions, ground_truth, snakemake.params.coverage_number)
    # print("nextclade validation: ", nextclade_validation)

    # #############################################
    
    # #############################################
    # #kallisto####################################
    # ############################################# 

    # #assign list of pangolin results to a variable
    # kallisto_input = list(snakemake.input.kallisto)
    # print(kallisto_input)

    # def extract_kallisto(list_of_input):
    #     print("KALLISTO")
    #     #initialize a list to collect predictions
    #     predictions=[]

    #     #loop over the files

    #     #sort the list of files
    #     list_of_input.sort(key=natural_keys)
    #     print("sorted list: ", list_of_input)

    #     for file in list_of_input:

    #         #predictions
    #         #read the file into a df
    #         df=pl.read_csv(file, separator="\t")

    #         #find the row for column 'target_id' (lineage) that has the largest value for column 'tpm'
    #         sorted_df=df.sort('tpm')
    #         print(sorted_df)
    #         lineage=sorted_df.select(pl.last('target_id'))
    #         lineage=lineage[0,'target_id']
    #         print(lineage)
    #         predictions.append(lineage)
    #     return predictions
    # #extract kallisto predictions
    # kallisto_predictions=extract_kallisto(kallisto_input)
    # print("kallisto predictions: ", kallisto_predictions)

    # #validate kallisto
    # kallisto_validation=validate_tool("kallisto", kallisto_predictions, ground_truth, snakemake.params.coverage_number)
    # print("kallisto validation: ", kallisto_validation)

    sorted_output = sorted(snakemake.output.validation, key=lambda x: int(x.split('_')[-1].replace('x', '').replace('.tsv', '')))
    [path_100x, path_1000x] = sorted_output
    print("path_100x", path_100x)
    print("path_1000x", path_1000x)

    for cov in ["100x", "1000x"]:
        #export 
        if cov == "100x":
            dfs = [orthanq_validation_100x]
            dfs_final = ft.reduce(lambda left, right: left.join(right, on='n_samples'), dfs)
            dfs_final.to_csv(path_100x)
        elif cov == "1000x":
            dfs = [orthanq_validation_1000x]
            dfs_final = ft.reduce(lambda left, right: left.join(right, on='n_samples'), dfs)
            dfs_final.to_csv(path_1000x)

    # all_tools.write_csv(snakemake.output.validation)
