#this script samples fractions with the input from the previous rule 'create_simulation_input' and creates sample compositions
#in separate tables and creates fasta files for each of them (?).

import pandas as pd
from pathlib import Path
from random import sample 
import numpy as np
import random

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f 

    n_of_samples=snakemake.params.n_of_samples

    #read input table
    df = pd.read_csv(snakemake.input.sequence_counts_usa_final)

    #sample N sets
    list_of_distributions = []
    list_of_distributions = [eval(row['clades_frequencies']) for _, row in df.iterrows()]
    random.seed(1313)

    #first sample some timepoints 
    sampled_values = sample([list(d.values()) for d in list_of_distributions],n_of_samples)

    #normalization function to normalize fractions in the final composition
    def normalize_dict_values(data):
        total = sum(data.values())
        return {key: value / total for key, value in data.items()}
    
    #create separate dataframes for each composition and write to csv
    #to do that, simulate combinations of 3 lineages sampled randomly propotional to case count abundances.
    
    sample_name_counter = 1
    for i,d in enumerate(list_of_distributions):
        for s in sampled_values:
            #when finding the index for the sample, create a csv for each synthetic sample
            if sorted(d.values()) == sorted(s):
                sample_dict_list = []
                print(list(d.keys()))
                print(s)
                #set random seed for numpy 
                np.random.seed(42)

                sampled_lineages = np.random.choice(list(d.keys()), p=s, size=3, replace=False)
                print("sampled_lineages:",sampled_lineages)
                
                #create sub dictionary with sampled lineages and normalize fractions
                subset_dict = {key: d[key] for key in sampled_lineages if key in d}
                normalized_subset_dict=normalize_dict_values(subset_dict)
                print("normalized_subset_dict:",normalized_subset_dict)

                for l in normalized_subset_dict:
                    sample_dict = {}
                    sample_dict['lineage'] = str(l.split(" ")[0])
                    sample_dict['fraction'] = normalized_subset_dict[l]
                    sample_dict['sample'] = "SimulatedSample" + str(sample_name_counter)
                    sample_dict_list.append(sample_dict)
                        
                sample_df = pd.DataFrame(sample_dict_list)
                sample_df.to_csv(snakemake.output.synthetic_samples[sample_name_counter-1], index=False, doublequote=False)

                sample_name_counter += 1