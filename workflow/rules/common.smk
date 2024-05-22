import pandas as pd

configfile: "config/config.yaml"

n_reads = config["n_reads"]

num_samples = config["number_of_samples"]

orthanq_sequences = ['23D', '20F', '22F', '21A', '23E', '21G', '21J', '22E', '24B', '20I', '21L', '23B', '21I', '22B', '20G', '23I', '21E', '20D', '20E', '21K', '23F', '22A', '21F', '23A', '20H', '22D', '23G', '23H', '24A', '20A', '21B', '21H', '22C', '20B', '21C', '20J', '23C', '21D', '19A', '20C', '21M', '19B']

num_list = [num+1 for num in range(num_samples)]

#modify sample tables to calculate required number of reads per lineage
def generate_numreads(lineages, n_reads):
    lineages['num_reads'] = lineages['fraction']*n_reads #assume the mixture will have 2000 reads
    lineages['num_reads'] = lineages['num_reads'].astype(int).astype('str')
    return lineages

# input function for get_fractions (executed after checkpoint 'create_sample_compositions')
# decision based on content of output file
def aggregate_input_get_fractions(wildcards):
    # Important: check if this works on a shared filesystem.
    all_samples = []
    for i in range(num_samples):
        with checkpoints.create_sample_compositions.get(**wildcards).output[i].open() as f:
            table = pd.read_csv(f)
            table_modified = generate_numreads(table, n_reads)
            all_samples.append(table_modified)
    all_samples_df = pd.concat(all_samples)
    print(all_samples_df)

    fq1=[f"results/fractions/{row["sample"]}-{row["lineage"]}-{row["num_reads"]}_1.fq" for i,row in all_samples_df.iterrows()]
    fq2=[f"results/fractions/{row["sample"]}-{row["lineage"]}-{row["num_reads"]}_2.fq" for i,row in all_samples_df.iterrows()]

    print(fq1)
    return fq1 + fq2

# input functions for concat_fractions (executed after checkpoint 'create_sample_compositions')
def aggregate_input_concat_fractions_fq1(wildcards):
    paths = []
    for i in range(num_samples):
        with checkpoints.create_sample_compositions.get(**wildcards).output[i].open() as f:
            table = pd.read_csv(f)
            table_modified = generate_numreads(table, n_reads)
            paths.append(table_modified)
    all_samples_df = pd.concat(paths)
    fq1=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_1.fq",
        zip,
        lineage=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['lineage'],
        num=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['num_reads']
        )
    return fq1

def aggregate_input_concat_fractions_fq2(wildcards):
    paths = []
    for i in range(num_samples):
        with checkpoints.create_sample_compositions.get(**wildcards).output[i].open() as f:
            table = pd.read_csv(f)
            table_modified = generate_numreads(table, n_reads)
            paths.append(table_modified)
    all_samples_df = pd.concat(paths)
    # print(all_samples_df)
    fq2=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_2.fq",
        zip,
        lineage=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['lineage'],
        num=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['num_reads']
        )
    return fq2

# input function to retrieve fastq samples
def get_fastq_input(wildcards):
    if config["simulation"] == False:
        sample = samples.loc[wildcards.sample]
        return [sample["fq1"], sample["fq2"]]
    else:
        sample = wildcards.sample
        simulated = ["results/mixed/{sample}_1.fq", "results/mixed/{sample}_2.fq"]
        return simulated