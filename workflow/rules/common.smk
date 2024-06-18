import pandas as pd

configfile: "config/config.yaml"

n_reads = config["n_reads"]

num_samples = config["number_of_samples"]

#removed 22F temporarily
orthanq_sequences = ['23D', '20F', '21A', '23E', '21G', '21J', '22E', '24B', '20I', '21L', '23B', '21I', '22B', '20G', '23I', '21E', '20D', '20E', '21K', '23F', '22A', '21F', '23A', '20H', '22D', '23G', '23H', '24A', '20A', '21B', '21H', '22C', '20B', '21C', '20J', '23C', '21D', '19A', '20C', '21M', '19B']

num_list = [num+1 for num in range(num_samples)]

#modify sample tables to calculate required number of reads per lineage
def generate_numreads(lineages, n_reads):
    lineages['num_reads'] = lineages['fraction']*n_reads #assume the mixture will have 2000 reads
    lineages['num_reads'] = lineages['num_reads'].astype(int).astype('str')
    return lineages


#following are only relevant when simulate_given is true
#for sample simulation from given lineage.tsv
lineages = pd.read_csv(config["lineages"], sep ="\t")

# the simulated sample is created below
simulated_sample = generate_numreads(lineages, n_reads)
print(simulated_sample)
#upstream are only relevant when simulate_given is true

# input function for get_fractions (executed after checkpoint 'create_sample_compositions')
# decision based on content of output file
def aggregate_input_get_fractions(wildcards):
    if config["simulate_given"] == False:
        # Important: check if this works on a shared filesystem.
        all_samples = []
        for i in range(num_samples):
            with checkpoints.create_sample_compositions.get(**wildcards).output[i].open() as f:
                table = pd.read_csv(f)
                table_modified = generate_numreads(table, n_reads)
                all_samples.append(table_modified)
        all_samples_df = pd.concat(all_samples)
        print(all_samples_df)
        fqs=[f"results/fractions/{row["sample"]}-{row["lineage"]}-{row["num_reads"]}_1.fq" for i,row in all_samples_df.iterrows()] + [f"results/fractions/{row["sample"]}-{row["lineage"]}-{row["num_reads"]}_2.fq" for i,row in all_samples_df.iterrows()]
        print(fqs)
        return fqs
    else:
        fqs=[f"results/fractions/{row['sample']}-{row['lineage']}-{row['num_reads']}_1.fq" for _, row in simulated_sample.iterrows()] + [f"results/fractions/{row['sample']}-{row['lineage']}-{row['num_reads']}_2.fq" for _, row in simulated_sample.iterrows()]
        return fqs


# input functions for concat_fractions (executed after checkpoint 'create_sample_compositions')
def aggregate_input_concat_fractions(wildcards):
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
    fq2=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}_2.fq",
        zip,
        lineage=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['lineage'],
        num=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['num_reads']
        )
    return [fq1,fq2]

# input function to retrieve fastq samples
def get_concat_fractions_input(wildcards):
    if config["simulate_given"] == False:
        fq1=aggregate_input_concat_fractions(wildcards)[0]
        fq2=aggregate_input_concat_fractions(wildcards)[1]
        return [fq1,fq2]
    else:
        fq1=expand("results/fractions/{{sample}}-{lineage}-{num}_1.fq",
        zip,
        lineage=simulated_sample.loc[simulated_sample['sample'] == wildcards.sample]['lineage'],
        num=simulated_sample.loc[simulated_sample['sample'] == wildcards.sample]['num_reads']
        )
        fq2=expand("results/fractions/{{sample}}-{lineage}-{num}_2.fq",
        zip,
        lineage=simulated_sample.loc[simulated_sample['sample'] == wildcards.sample]['lineage'],
        num=simulated_sample.loc[simulated_sample['sample'] == wildcards.sample]['num_reads']
        )
        return [fq1,fq2]

# input function to retrieve fastq samples
def get_fastq_input(wildcards):
    sample = wildcards.sample
    simulated = ["results/mixed/{sample}_1.fastq", "results/mixed/{sample}_2.fastq"]
    return simulated

def get_results(wildcards):
    if config["simulate_given"] == False:
        final_output = expand("results/orthanq/calls/SimulatedSample{num}/SimulatedSample{num}.tsv", num=num_list) + expand("results/orthanq/calls/SimulatedSample{num}/viral_solutions.html", num=num_list)
        return final_output
    else:
        final_output = expand("results/orthanq/calls/{sample}/{sample}.tsv", sample=simulated_sample["sample"].unique()) + expand("results/orthanq/calls/{sample}/viral_solutions.html", sample=simulated_sample["sample"].unique())
        return final_output

#pseudodate tobe used in uncovar workflow
DATE="13062024"
UNCOVAR_SAMPLE_SHEET="uncovar/config/pep/samples.csv"
UNCOVAR_CONFIG="uncovar/config/config.yaml"
UNCOVAR_PEP_CONFIG="uncovar/config/pep/config.yaml"