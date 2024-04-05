# This script take Per Country Case Counts data as input. Then it does the following as listed: 
#
# 1-) It extracts only the "canonical" (e.g. Alpha, Beta etc. and not 'others', etc.) Nextstrain clades.
# 2-) It calculates frequency of each clade for all records, and it does that by counting the total 
# number of cases that are included in the previous step and dividing number of cases for each clade to this number.
# 3-) It removes records with case numbers less than 100, removes clades with 0.0 freq from each record.
# 4-) It only includes records not more than 6 clades and the remaining records are sampled to create a dataset with N (=15).

import pandas as pd
from pathlib import Path
import json
import os
from random import sample 

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #read the json and load its values 
    jsonFile = open(snakemake.input.per_country_data_json, 'r')
    values = json.load(jsonFile)
    jsonFile.close()

    #take country distributions and convert to df
    distributions = values['regions'][0]['distributions']
    filtered = [x for x in distributions if (x['country'] == 'USA')]
    usa_dist = filtered[0]["distribution"]
    usa_dist_df = pd.DataFrame(usa_dist)

    #write distributions to tabular file
    usa_dist_df.to_csv(snakemake.output.per_country_data_csv_raw, index=False)

    #iterate over the dataframe to calculate frequencies for each clade
    new_usa_dist = []
    for index, row in usa_dist_df.iterrows():
        clades = row['stand_estimated_cases']
        n_total_cases = row['stand_total_cases']

        #zero total cases are filtered out
        if n_total_cases == 0:
            continue
        else:
            new_stand_clade_cases = 0
            new_clades_with_case_numbers = {}
            for clade in clades:
                #insert only records that contain alpha, beta, gamma, delta, omicron, kappa, epsilon, eta, iota, lambda, mu, omicron
                #sum case numbers as the new total number of cases
                clades_to_look = ["Alpha", "Beta", "Gamma", "Delta", "Omicron", "Kappa", "Epsilon", "Eta", "Iota", "Lambda", "Mu", "Omicron"]

                if any(query_clade in clade for query_clade in clades_to_look):
                    new_stand_clade_cases += clades[clade]
                    new_clades_with_case_numbers[clade] = clades[clade]
            
            if new_stand_clade_cases != 0:
                #calculate frequencies for each clade
                new_clades_with_freqs = {}
                for clade in new_clades_with_case_numbers:
                    freq = new_clades_with_case_numbers[clade]/new_stand_clade_cases
                    new_clades_with_freqs[clade] = freq

                #finally append to the main dictionary
                new_usa_dist.append({'percent_total_cases': row['percent_total_cases'], 
                                    'total_sequences': row['total_sequences'],
                                    'week': row['week'],
                                    'new_stand_estimated_cases': new_clades_with_freqs, 
                                    'stand_total_cases': new_stand_clade_cases})
    new_usa_dist_df = pd.DataFrame(new_usa_dist)
    print(new_usa_dist_df)

    #then, remove stand_total_cases less than 100
    for index, row in new_usa_dist_df.iterrows():
        if row['stand_total_cases'] < 100:
            new_usa_dist_df.drop(index=[index], inplace=True)

    print(new_usa_dist_df)

    #then, remove the clade from clades dict if it has 0.0
    for index, row in new_usa_dist_df.iterrows():
        clades = row['new_stand_estimated_cases']
        for clade in clades.copy():
            if clades[clade] == 0.0:
                clades.pop(clade, None)
    print(new_usa_dist_df)

    #then, only select rows with clades size not more than 6
    for index, row in new_usa_dist_df.iterrows():
        clades = row['new_stand_estimated_cases']
        if len(clades) > 6:
            new_usa_dist_df.drop(index=[index], inplace=True)

    print(new_usa_dist_df)

    #sample 20 sets
    list_of_freqs = []
    for index, row in new_usa_dist_df.iterrows():
        clades = row['new_stand_estimated_cases']
        list_of_freqs.append(clades.values())

    n_of_samples=15 #can be configured

    print(sample(list_of_freqs,n_of_samples)) 

    #write distributions to tabular file after updating the clade values with frequencies
    new_usa_dist_df.to_csv(snakemake.output.per_country_data_csv_final, index=False)
