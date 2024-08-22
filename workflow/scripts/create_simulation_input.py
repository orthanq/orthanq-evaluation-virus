# This script take Per Country Case Counts data as input. Then it does the following as listed: 
#
# 1-) It calculates frequency of each clade for all records, and it does that by counting the total 
# number of cases that are included in the previous step and dividing number of cases for each clade to this number.
# 2-) It removes records with case numbers less than 50
# 3-) Then, it sorts frequencies and only takes lineages with first n numbers.
# 4-) TODO: Finally, we should probably add some other fraction to sum up to 1.

import pandas as pd
import polars as pl

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #read csv
    counts = pl.read_csv(snakemake.input.sequence_counts_usa, separator="\t")
    # counts = pl.read_csv("usa.tsv", separator="\t")

    #firstly sum number of sequences for each date and clade (because there are sequences from different states in usa)
    summed_df = counts.group_by(['date', 'clade']).agg([
        pl.col('sequences').sum()
    ])

    #then, group by each date and aggregate clade-sequence pairs to a list for each clade

    #function to convert keys and values to dictionary
    def lists_to_dict(keys, values):
        return dict(zip(keys, values))

    #now group by 'date' and aggregate clade and sequences
    agg_df = (
        summed_df.group_by('date')
        .agg([
            pl.col('clade'),
            pl.col('sequences')
        ])
        .with_columns(
            pl.struct(['clade', 'sequences']).map_elements(
                lambda row: lists_to_dict(row['clade'], row['sequences']),
                return_dtype=pl.Object
            ).alias('dict')
        )
    )

    #create a column for total number of sequences by collecting sequences
    agg_df_total = agg_df.with_columns(
        pl.col('sequences').map_elements(lambda x: sum(x),return_dtype=pl.Object).alias('total_sequences')
    )

    #write to csv
    # agg_df_total.write_csv(snakemake.output.sequence_counts_usa_raw)

    #now, iterate over the dataframe to calculate frequencies for each clade
    usa_freq_dist = []
    for row in agg_df_total.iter_rows(named=True):
        clades = row['dict']
        n_total_cases = row['total_sequences']
        #zero total cases are filtered out
        if n_total_cases == 0:
            continue
        else:
            new_total_cases = 0
            clades_counts_dict = {}
            for clade in clades:
                new_total_cases += clades[clade]
                clades_counts_dict[clade] = clades[clade]
            
            if new_total_cases != 0:
                #calculate frequencies for each clade
                clades_freqs = {}
                for clade in clades_counts_dict:
                    freq = clades_counts_dict[clade]/new_total_cases
                    clades_freqs[clade] = freq

                #sort the dictionary and take first n elements
                n=3
                sorted_dict = dict(sorted(clades_freqs.items(), key=lambda item: item[1], reverse=True)[:n])

                #finally append to the main dictionary
                usa_freq_dist.append({'total_sequences': row['total_sequences'],
                                    'week': row['date'],
                                    'clades_frequencies': sorted_dict})
                
    #converting directly to polars somehow creates issues with clades_freqs, it can't use the dict type for some reason.
    new_usa_dist_df = pd.DataFrame(usa_freq_dist)


    ##NOW DO SOME FILTERING

    #then, remove total_case_numbers less than 50
    for index, row in new_usa_dist_df.iterrows():
        if row['total_sequences'] < 50:
            new_usa_dist_df.drop(index=[index], inplace=True)

    #write to csv
    # new_usa_dist_df.to_csv("final_filter.csv", index=False)
    # print(new_usa_dist_df)
    new_usa_dist_df.to_csv(snakemake.output.sequence_counts_usa_final, index=False)

