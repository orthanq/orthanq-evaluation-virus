import pandas as pd
import altair as alt
import numpy as np
import pyarrow
import os
import json
import sys
import collections
import scipy.stats
import pysam
from itertools import combinations

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    #create a dictionary with sample name and lineages that it contains
    simulation = snakemake.input.simulation

    sample_lineages = {}
    for file in simulation:
        sample_name = os.path.basename(file).split(".")
        results = pd.read_csv(file)
        lineages = results['lineage'].tolist()
        print("lineages",lineages)

        sample_lineages[sample_name[0]] = lineages

    print(sample_lineages)

    #calculate pairwise variant differences for each simulated sample

    vcf_path = snakemake.input.candidates 
    vcf = pysam.VariantFile(vcf_path)

    # mapvlineage -> set of variants (POS_REF_ALT)
    lineage_variants = {}

    # collect variants per lineage based on sample name
    all_lineages = {ln for v in sample_lineages.values() for ln in v}
    print("known_lineages",all_lineages)
    for record in vcf:
        for sample_name, sample_data in record.samples.items():
            for lineage in all_lineages:
                if  sample_name == lineage:
                    # Check if genotype is 1|1 or (1,1)
                    gt = sample_data.get('GT')
                    if gt == (1, 1):
                        var_id = f"{record.pos}_{record.ref}_{record.alts[0]}"
                        lineage_variants.setdefault(lineage, set()).add(var_id)
    print("lineage_variants",lineage_variants)

    # for each simulated sample, compute pairwise differences
    results = []

    for sim_sample, lineages in sample_lineages.items():
        for lin1, lin2 in combinations(lineages, 2):
            variants1 = lineage_variants.get(lin1, set())
            variants2 = lineage_variants.get(lin2, set())
            
            shared = variants1 & variants2
            unique1 = variants1 - variants2
            unique2 = variants2 - variants1
            diff = len(unique1) + len(unique2)

            results.append({
                'SimulatedSample': sim_sample,
                'Lineage1': lin1,
                'Lineage2': lin2,
                'Unique_Lineage1': len(unique1),
                'Unique_Lineage2': len(unique2),
                'Shared': len(shared),
                'Total_Differences': diff
            })
    # convert to dataframe and write to file
    df = pd.DataFrame(results)
    print(df)
    df.to_csv(snakemake.output.table, index=False)

    #bubble plot
    df["Lineage_Pair"] = df["Lineage1"] + " vs " + df["Lineage2"]
    bubble_plot = alt.Chart(df).mark_circle().encode(
        x=alt.X('SimulatedSample:N', title="Sample"),
        y=alt.Y('Lineage_Pair:N', title="Lineage Pair"),
        size=alt.Size('Total_Differences:Q', scale=alt.Scale(range=[10, 1000]), legend=alt.Legend(title="Total Differences")),
        color=alt.Color('Total_Differences:Q', scale=alt.Scale(scheme='viridis'), legend=None),
        tooltip=[
            alt.Tooltip('SimulatedSample:N', title="Sample"),
            alt.Tooltip('Lineage_Pair:N', title="Lineages"),
            alt.Tooltip('Total_Differences:Q', title="Differences")
        ]
    ).properties(
        width=500,
        height=300,
        title="Total Differences Between Lineage Pairs per Sample"
    )

    bubble_plot

    bubble_plot.save(snakemake.output.plot_svg)
    bubble_plot.save(snakemake.output.plot_html)