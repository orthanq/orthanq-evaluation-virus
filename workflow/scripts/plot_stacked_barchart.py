import pandas as pd
import altair as alt
import numpy as np
import matplotlib.pyplot as plt
import pyarrow
import os

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    def get_orthanq_results(file_paths):
        #create a dataframe with sample names, fractions and predictions in categories below.
        orthanq_predictions = pd.DataFrame(columns=('sample', 'prediction', 'WHO_name', 'fraction', 'tool'))

        for file in file_paths:
            #find sample name
            sample = os.path.basename(file).split(".")[0]
            print("sample",sample)
                    
            #read table
            results = pd.read_csv(file)

            #find the row with best odds
            best_odds = [1, 1.0, 1.00]
            best_results = results[results['odds'].isin(best_odds)]
        
            #remove density and odds columns
            best_results_filtered = best_results.drop(columns=['density', 'odds'])

            #both studies (PRJNA804575 and PRJNA809680) contain only one row with the best density.
            #loop over the lineages and fractions and group them into "Delta", "Omicron", "Other" and "Recombinant".
            if len(best_results_filtered) == 1:
                best_results_dict = best_results_filtered.iloc[0].to_dict()  # Convert to a dictionary
                print(best_results_dict)
                #Delta includes B.1.617.2 and all lineages starting with AY (CITE: Bolze et al 2022)
                #Omicron lineages refer to BA starting variants, see https://covariants.org/ for other omicron lineages.
                #B.1.1.161 & B.1.189 (early pandemic lineages), B (very early ancestral lineage)
                #X starting lineages refer to recombinant.
                #All orthanq results are checked by eye for an unassigned lineage.
                for lineage in best_results_dict:
                    if best_results_dict[lineage] != 0.0:
                        #DELTA 
                        if lineage.startswith("AY") or lineage == "B.1.617.2":
                            new_row = pd.DataFrame([[sample, lineage, "delta", best_results_dict[lineage], 'orthanq']], columns=['sample', 'prediction', 'WHO_name', 'fraction', 'tool'])
                            orthanq_predictions = pd.concat([orthanq_predictions, new_row], ignore_index=True)
                        #OMICRON 
                        elif lineage.startswith("BA"):
                            new_row = pd.DataFrame([[sample, lineage, "omicron", best_results_dict[lineage], 'orthanq']], columns=['sample', 'prediction', 'WHO_name', 'fraction', 'tool'])
                            orthanq_predictions = pd.concat([orthanq_predictions, new_row], ignore_index=True)
                        #RECOMBINANT 
                        elif lineage.startswith("X"):
                            new_row = pd.DataFrame([[sample, lineage, "recombinant", best_results_dict[lineage], 'orthanq']], columns=['sample', 'prediction', 'WHO_name', 'fraction', 'tool'])
                            orthanq_predictions = pd.concat([orthanq_predictions, new_row], ignore_index=True)
                        #OTHER
                        else:
                            new_row = pd.DataFrame([[sample, lineage, "other", best_results_dict[lineage], 'orthanq']], columns=['sample', 'prediction', 'WHO_name', 'fraction', 'tool'])
                            orthanq_predictions = pd.concat([orthanq_predictions, new_row], ignore_index=True)
            else:
                print("DataFrame has multiple rows or is empty.")
        return orthanq_predictions  

    def get_pangolin_results(file_paths):
        data = []
        
        for file in file_paths:
            # Extract sample name from filename
            sample = os.path.basename(file).split('.')[0]
            
            # Read the CSV file
            df = pd.read_csv(file)
            
            # Extract relevant columns
            for _, row in df.iterrows():
                prediction = row['lineage']
                WHO_name = row['scorpio_call']
                
                # Adjust WHO_name 
                if WHO_name == "Probable Omicron (Unassigned)" or WHO_name == "Omicron (Unassigned)" or WHO_name == "Omicron (BA.1-like)":
                    WHO_name = "omicron"
                elif WHO_name == "Probable Delta (Unassigned)" or WHO_name == "Delta (B.1.617.2-like)":
                    WHO_name = "delta"
                
                data.append([sample, prediction, WHO_name, 1.0, 'pangolin'])
        
        # Create DataFrame
        result_df = pd.DataFrame(data, columns=['sample', 'prediction', 'WHO_name', 'fraction', 'tool'])
        
        return result_df

    def get_nextclade_results(file_paths):
        data = []
        
        for file in file_paths:
            # Extract the sample name from the folder name (parent directory)
            sample = os.path.basename(os.path.dirname(file))
            
            # Read the CSV file
            df = pd.read_csv(file, delimiter=';')  # Nextclade uses semicolon as delimiter
            
            # Extract relevant columns
            for _, row in df.iterrows():
                prediction = row['Nextclade_pango']
                WHO_name = row['clade_who'].lower()
                
                data.append([sample, prediction, WHO_name, 1.0, 'nextclade'])
        
        # Create DataFrame
        result_df = pd.DataFrame(data, columns=['sample', 'prediction', 'WHO_name', 'fraction', 'tool'])
        
        return result_df

    def get_kallisto_results(csv_files) -> pd.DataFrame:

        data = []

        for file_path in csv_files:
            # Extract sample name (SRR ID) from the folder path
            sample_name = os.path.basename(os.path.dirname(file_path)).split("_")[-1]

            # Read the CSV file
            df = pd.read_csv(file_path, sep="\t")  # Assuming TSV format

            # Sort by 'tpm' in descending order
            df_sorted = df.sort_values(by='tpm', ascending=False)

            # Select the top 5 rows with the highest TPM values
            top5 = df_sorted.head(5)

            # Calculate the total TPM of the top 5
            total_tpm = top5['tpm'].sum()

            # Compute the fraction (percentage) of each prediction in the top 5 mix
            top5['fraction'] = top5['tpm'] / total_tpm

            # Rename 'target_id' to 'prediction'
            top5 = top5.rename(columns={'target_id': 'prediction'})

            # Add 'sample' and 'tool' columns
            top5['sample'] = sample_name
            top5['tool'] = "kallisto"

            #set WHO names
            #check kallisto names manually to see if there are no exceptions.
            for index, row in top5.iterrows():
                #DELTA 
                if row["prediction"].startswith("AY") or row["prediction"] == "B.1.617.2":
                    # top5['WHO_name']="delta"
                    top5.at[index,'WHO_name'] = "delta"
                #OMICRON 
                elif row["prediction"].startswith("BA"):
                    top5.at[index,'WHO_name'] = "omicron"
                #RECOMBINANT 
                elif row["prediction"].startswith("X"):
                    top5.at[index,'WHO_name'] = "recombinant"
                #OTHER
                else:
                    top5.at[index,'WHO_name'] = "other"

            # Append to the results list
            data.append(top5[['sample', 'prediction', 'WHO_name' , 'fraction', 'tool']])

        # Combine all results into a single DataFrame
        return pd.concat(data, ignore_index=True)

    #read the truth
    truth_table = pd.read_csv(snakemake.input.truth, sep="\t")
    print(truth_table)

    ### ORTHANQ ###
    orthanq_predictions = get_orthanq_results(snakemake.input.orthanq)

    ### PANGOLIN ###
    pangolin_predictions = get_pangolin_results(snakemake.input.pangolin)

    ### NEXTCLADE ###
    nextclade_predictions = get_nextclade_results(snakemake.input.nextclade)

    ### KALLISTO ###
    kallisto_predictions = get_kallisto_results(snakemake.input.kallisto)

    #final_table
    final_table = pd.concat([truth_table, orthanq_predictions, pangolin_predictions, nextclade_predictions, kallisto_predictions], ignore_index=True)

    # Function to generate shades of a base color
    def generate_shades(cmap_name, num_shades):
        cmap = plt.get_cmap(cmap_name)  # Get colormap
        
        # Define how close shades should be based on count
        if num_shades <= 3:
            shade_positions = np.linspace(0.5, 0.8, num_shades)  # Closer shades
        else:
            shade_positions = np.linspace(0.3, 0.9, num_shades)  # Wider range for large numbers
        
        return [cmap(pos) for pos in shade_positions]

    # Convert RGB to HEX
    def rgb_to_hex(rgb):
        return '#%02x%02x%02x' % tuple(int(c * 255) for c in rgb[:3])  # Convert RGB tuple to HEX

    # Get unique lineages per class
    lineages_omicron = final_table[final_table['WHO_name'] == 'omicron']['prediction'].unique()
    lineages_delta = final_table[final_table['WHO_name'] == 'delta']['prediction'].unique()
    lineages_rec = final_table[final_table['WHO_name'] == 'recombinant']['prediction'].unique()
    lineages_other = final_table[final_table['WHO_name'] == 'other']['prediction'].unique()

    # Generate color shades dynamically
    shades_omicron = generate_shades('Reds', len(lineages_omicron))  # Shades of red
    shades_delta = generate_shades('Purples', len(lineages_delta))  # Shades of blue
    shades_rec = generate_shades('Greens', len(lineages_rec))  # Shades of blue
    shades_other = generate_shades('Grays', len(lineages_other))  # Shades of blue

    # Create color mapping
    color_mapping = {**dict(zip(lineages_omicron, map(rgb_to_hex, shades_omicron))), **dict(zip(lineages_delta, map(rgb_to_hex,shades_delta))), **dict(zip(lineages_rec, map(rgb_to_hex,shades_rec))), **dict(zip(lineages_other, map(rgb_to_hex,shades_other)))}

    # Create the stacked bar chart with grouped bars
    chart = alt.Chart(final_table).mark_bar().encode(
        x=alt.X('tool:N', title=None, axis=alt.Axis(labelAngle=-45),sort=['truth', 'orthanq', 'kallisto', 'pangolin', 'nextclade']),   
        # xOffset='tool:N',
        y=alt.Y('sum(fraction)', title="fraction" ,scale=alt.Scale(domain=[0.0, 1.0], clamp=True), axis=alt.Axis(grid=False)),
        column=alt.Column('sample:N', header=alt.Header(titleFontSize=24), title=None),
        # color=alt.Color('WHO_name:N'),
        color=alt.Color(
        'prediction:N',
        legend=alt.Legend(symbolLimit=100),
        scale=alt.Scale(domain=list(color_mapping.keys()), range=list(color_mapping.values())),
        # legend=alt.Legend(title='Lineage Colors')
    ),
        tooltip=['sample', 'tool', 'WHO_name', 'fraction']
    ).configure_view(
    strokeOpacity=0    
    )
    # save the chart
    chart.save(snakemake.output.plot_svg)
    chart.save(snakemake.output.plot_html)

    #write the table
    final_table.to_csv(
        snakemake.output.table, sep="\t", index=False, header=True
    )
    