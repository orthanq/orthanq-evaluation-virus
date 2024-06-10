#validation with the most abundant lineages? Analysis in terms of call rate and accuracy.

import polars as pl
import pandas as pd
import os
import sys
import functools as ft

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    ###########################################
    #ground truth##############################
    ###########################################

    #simulation input
    simulation = snakemake.input.simulation

    def extract_ground_truth(simulation):
        #collect the abundant lineage in simulations
        simulation_largest_lineages=[]

        for file in sorted(simulation):
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
        return simulation_largest_lineages
    
    ground_truth = extract_ground_truth(simulation)
    print(ground_truth)

    def validate_tool(tool_name, sorted_tool_predictions, sorted_ground_truth):
        #we see that all tools always have a call, because we never saw a case without in viral lineage quantification
        #but if that's not the case, raise an exception
        if len(sorted_tool_predictions)!=len(sorted_ground_truth):
            raise Exception("call rate cannot be 1 because of different lengths")
        call_rate=1.0
        n_correct=0
        #calculate call rate and accuracy
        for pred,sim in zip(sorted_tool_predictions,sorted_ground_truth):
            if pred==sim:
                n_correct+=1

        accuracy=n_correct/len(sorted_tool_predictions)*100 #it can also be the length of simulated samples, but since they're the same, it doesn't matter which of them is used here
        print(accuracy)

        #create column names with tools
        call_rate_column=tool_name+"_"+"call_rate"
        accuracy_column=tool_name+"_"+"accuracy"

        #convert to dataframe and write
        tool_validation = pl.DataFrame({"n_samples": len(sorted_tool_predictions), call_rate_column: call_rate, accuracy_column: accuracy})

        return tool_validation
    
    ############################################
    
    ############################################
    #orthanq####################################
    ############################################
    
    #orthanq results
    orthanq_input = snakemake.input.orthanq

    def extract_orthanq(list_of_input):
        #initialize orthanq predictions
        predictions=[]

        #collect the abundandant lineage in predictions
        for file in sorted(list_of_input):
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
            predictions.append(max_value)
        return predictions
    
    #extract orthanq results
    orthanq_predictions=extract_orthanq(orthanq_input)
    print("orthanq_predictions: ",orthanq_predictions)

    #validate orthanq
    orthanq_validation = validate_tool("orthanq", orthanq_predictions, ground_truth)
    print("orthanq_validation",orthanq_validation)

    #############################################
    
    #############################################
    #pangolin####################################
    #############################################

    #assign list of pangolin results to a variable
    pangolin_input = snakemake.input.pangolin

    def extract_pangolin(list_of_input):
        #initialize a list to collect predictions
        predictions=[]

        #loop over the files
        for file in sorted(list_of_input):

            #predictions
            #read the file into a df
            df=pl.read_csv(file)

            #pangolin only output one row
            df=df[0]

            #extract the prediction in 'lineage' column 
            predictions.append(df[0,"lineage"])

        return predictions
    #extract pangolin predictions
    pangolin_predictions=extract_pangolin(pangolin_input)
    print("pangolin predictions: ", pangolin_predictions)

    #validate pangolin (convert pango lineages to nextstrain, 
    #give a score to pango prediction if the predicted lineage has the correspondence that matches the ground truth)
    clade_to_lineage=snakemake.input.clade_to_lineage

    def validate_pangolin(pangolin_predictions, clade_to_lineage, sorted_ground_truth):

        #convert pango lineages to nextstrain clades
        clade_to_lineage=pl.read_csv(clade_to_lineage, has_header=False, separator="\t")
        print(clade_to_lineage)

        #first rename the columns
        clade_to_lineage=clade_to_lineage.rename({"column_1": "clade", "column_2": "pango"})

        #loop over pangolin predictions to calculate accuracy
        n_correct=0
        for idx,pango_lineage in enumerate(pangolin_predictions):
            #find the rows with pango_lineage in clade_to_lineage
            filtered=clade_to_lineage.filter(pl.col("pango") == pango_lineage)
            print(filtered)
            
            #create a list for all corresponding clades (can be more than one)
            corresp_clades=list(filtered['clade'])
            print(corresp_clades)

            #calculate call rate and accuracy
            truth_at_index=ground_truth[idx]
            print("truth_at_index: ", truth_at_index)
            print("corresp_clades: ", corresp_clades)
            if truth_at_index in corresp_clades:
                print("test")
                n_correct+=1
        call_rate=1.0
        accuracy=n_correct/len(pangolin_predictions)*100 #it can also be the length of simulated samples, but since they're the same, it doesn't matter which of them is used here
        print(accuracy)

        #convert to dataframe and write
        pangolin_validation = pl.DataFrame({"n_samples": len(pangolin_predictions), "pangolin_call_rate": call_rate, "pangolin_accuracy": accuracy})

        return pangolin_validation
    
    pangolin_validation=validate_pangolin(pangolin_predictions, clade_to_lineage, ground_truth)        
    # pangolin_validation=validate_tool("pangolin", pangolin_predictions, ground_truth)
    print("pangolin validation: ", pangolin_validation)

    #############################################
    
    #############################################
    #nextclade###################################
    #############################################    

    #assign list of pangolin results to a variable
    nextclade_input = snakemake.input.nextclade
    print(nextclade_input)
    def extract_nextclade(list_of_input):
        #initialize a list to collect predictions
        predictions=[]

        #loop over the files
        for file in sorted(list_of_input):

            #predictions
            #read the file into a df
            df=pl.read_csv(file, separator="\t")

            #nextclade only output one row
            df=df[0]

            #extract the prediction in 'clade' column 
            predictions.append(df[0,"clade"])

        return predictions
    #extract nextclade predictions
    nextclade_predictions=extract_nextclade(nextclade_input)
    print("nextclade predictions: ", nextclade_predictions)

    #validate pangolin
    nextclade_validation=validate_tool("nextclade", nextclade_predictions, ground_truth)
    print("nextclade validation: ", nextclade_validation)

    #############################################
    
    #############################################
    #kallisto####################################
    ############################################# 

    #assign list of pangolin results to a variable
    kallisto_input = snakemake.input.kallisto
    print(kallisto_input)

    def extract_kallisto(list_of_input):
        print("KALLISTO")
        #initialize a list to collect predictions
        predictions=[]

        #loop over the files
        for file in sorted(list_of_input):

            #predictions
            #read the file into a df
            df=pl.read_csv(file, separator="\t")

            #find the row for column 'target_id' (lineage) that has the largest value for column 'tpm'
            sorted_df=df.sort('tpm')
            print(sorted_df)
            lineage=sorted_df.select(pl.last('target_id'))
            lineage=lineage[0,'target_id']
            print(lineage)
            predictions.append(lineage)
        return predictions
    #extract kallisto predictions
    kallisto_predictions=extract_kallisto(kallisto_input)
    print("kallisto predictions: ", kallisto_predictions)

    #validate kallisto
    kallisto_validation=validate_tool("kallisto", kallisto_predictions, ground_truth)
    print("kallisto validation: ", kallisto_validation)

    #merge all validation dataframes
    dfs = [orthanq_validation, pangolin_validation, nextclade_validation, kallisto_validation]
    df_final = ft.reduce(lambda left, right: left.join(right, on='n_samples'), dfs)
    print(df_final)

    # all_tools.write_csv(snakemake.output.validation)
    df_final.write_csv(snakemake.output.validation)
