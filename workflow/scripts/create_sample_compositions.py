#this script samples fractions with the input from the previous rule 'create_simulation_input' and creates sample compositions
#in separate tables and creates fasta files for each of them (?).

import pandas as pd
from pathlib import Path
from random import sample 
import random

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f 

    n_of_samples=snakemake.params.n_of_samples

    #read input table
    df = pd.read_csv(snakemake.input.per_country_data_csv_final)

    #sample N sets
    list_of_distributions = []
    list_of_distributions = [eval(row['new_stand_estimated_cases']) for _, row in df.iterrows()]

    random.seed(1313)

    sampled_values = sample([d.values() for d in list_of_distributions],n_of_samples)

    #create separate dataframes for each composition and write to csv
    sample_name_counter = 1
    for i,d in enumerate(list_of_distributions):
        for s in sampled_values:
            #when finding the index for the sample, create a csv for each synthetic sample
            if sorted(d.values()) == sorted(s):
                sample_dict_list = []
                for key in d:
                    sample_dict = {}
                    # print(str(key))
                    # print(str(key).split(" ")[0])
                    sample_dict['lineage'] = str(key.split(" ")[0])
                    sample_dict['fraction'] = d[key]
                    sample_dict['sample'] = "SimulatedSample" + str(sample_name_counter)
                    sample_dict_list.append(sample_dict)
                sample_df = pd.DataFrame(sample_dict_list)
                sample_df.to_csv(snakemake.output.synthetic_samples[sample_name_counter-1], index=False, doublequote=False)

                sample_name_counter += 1