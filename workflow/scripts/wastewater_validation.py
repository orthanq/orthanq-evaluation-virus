#check 2 lineages found in Baajens paper (B.1.1.7 and B.1.526), also see if we can find the additional 4 lineages that was presented in their supplementary.
#note: orthanq may output insufficient observations hence empty outputs. question: should we count them as 0.0? currently that's what we do.

import pandas as pd
import altair as alt
import os
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    #this function checks for the abundance of the given Pango lineage in orthanq results and returns the abundance, if not found, 0.0 is returned.
    def orthanq_check_abundance(lineage):
        sample_abundances={}
        for file in snakemake.input.orthanq:
            #find sample name
            sample = os.path.basename(file).split(".")[0]
            print("sample",sample)
                    
            #read table
            results = pd.read_csv(file)

            #names of haplotypes in the results
            haplotypes = results.columns.tolist()
            print("haplotypes", haplotypes)

            abundance = 0.0
            if lineage in haplotypes:
                #find the records that have the same density
                best_odds = [1, 1.0, 1.00]
                best_results = results[results['odds'].isin(best_odds)]
                print("---best_results---")
            
                #remove density and odds columns
                best_results_filtered = best_results.drop(columns=['density', 'odds'])
                for (i,row) in best_results_filtered.iterrows():
                    row_abun = row[lineage]
                    if row_abun > abundance:
                        abundance=row_abun
            sample_abundances[sample]=abundance*100
        return sample_abundances

    lineages = ["B.1.526","B.1.1.7"]
    
    #create a barplot with dates in the x axis and abundances in the y axis
    metadata=pd.read_csv(snakemake.input.metadata, sep=r'\s+', header=0)

    for l in lineages:
        orthanq_lineage_abundances=orthanq_check_abundance(l)
        print(orthanq_lineage_abundances)

        # Map the abundance values from the dictionary to the DataFrame
        metadata['Abundance'] = metadata['SRA_ID'].map(orthanq_lineage_abundances)
        print(metadata)

        chart_title =  'Abundance Over Time ({})'.format(l) 
        chart = alt.Chart(metadata).mark_bar().encode(
        x=alt.X('CollectionDate:T', title='Date', axis=alt.Axis(format='%Y-%m-%d', labelAngle=-45, tickCount='day', labelFontSize=5, labelOverlap='greedy'))   , 
        y=alt.Y('Abundance:Q',title='Abundance (%)'),   
        tooltip=['SRA_ID:N', 'CollectionDate:T', 'Abundance:Q']).properties(
        title=chart_title)
        
        if l == "B.1.1.7":
            chart.save(snakemake.output.validation_B117)
        elif l == "B.1.526":
            chart.save(snakemake.output.validation_B1526)