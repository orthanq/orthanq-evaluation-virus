import pandas as pd

configfile: "config/config.yaml"

n_reads_dict = {"100x": config["n_reads_100x"],"1000x": config["n_reads_1000x"]}
print(n_reads_dict)

num_samples = config["number_of_samples"]

num_list = [num+1 for num in range(num_samples)]

#modify sample tables to calculate required number of reads per lineage
def generate_numreads(lineages, n_reads):
    lineages['num_reads'] = lineages['fraction']*n_reads 
    lineages['num_reads'] = lineages['num_reads'].astype(int).astype('str')
    return lineages


#following are only relevant when simulate_given is true (only simulate a sample with given lineages and fractions in lineages.tsv)
#for sample simulation from given lineage.tsv
lineages = pd.read_csv(config["lineages_and_fractions"], sep ="\t")

# the simulated sample is created below
coverage_single_sample=config["coverage_simulate_given"]
simulated_given_lineages = generate_numreads(lineages, n_reads_dict[coverage_single_sample])
print("simulated_given_lineages",simulated_given_lineages)
#upstream are only relevant when simulate_given is true

# input function for get_fractions (executed after checkpoint 'create_sample_compositions')
# decision based on content of output file
def aggregate_input_get_fractions(wildcards):
    if config["simulate_pandemics"] and not config["simulate_given"]:
        all_samples_fqs = []
        for (coverage,n_reads) in n_reads_dict.items():
            # Important: check if this works on a shared filesystem.
            samples_per_coverage = []
            for i in range(num_samples):
                with checkpoints.create_sample_compositions.get(**wildcards).output[i].open() as f:
                    table = pd.read_csv(f)
                    table_modified = generate_numreads(table, n_reads)
                    samples_per_coverage.append(table_modified)
            samples_per_coverage_df = pd.concat(samples_per_coverage)
            fqs=[f"results/fractions/{row['sample']}-{row['lineage']}-{row['num_reads']}-{coverage}_1.fq" for i,row in samples_per_coverage_df.iterrows()] + [f"results/fractions/{row['sample']}-{row['lineage']}-{row['num_reads']}-{coverage}_2.fq" for i,row in samples_per_coverage_df.iterrows()]
            print("samples_per_coverage_df", samples_per_coverage_df)
            all_samples_fqs.extend(fqs)
        print("all_samples_fqs",all_samples_fqs)
        return all_samples_fqs
    elif config["simulate_given"] and not config["simulate_pandemics"]:
        fqs=[f"results/fractions/{row['sample']}-{row['lineage']}-{row['num_reads']}-{coverage_single_sample}_1.fq" for _, row in simulated_given_lineages.iterrows()] + [f"results/fractions/{row['sample']}-{row['lineage']}-{row['num_reads']}-{coverage_single_sample}_2.fq" for _, row in simulated_given_lineages.iterrows()]
        return fqs
    else:
        return []


# input functions for concat_fractions (executed after checkpoint 'create_sample_compositions')
def aggregate_input_concat_fractions(wildcards):
    paths = []
    for i in range(num_samples):
        with checkpoints.create_sample_compositions.get(**wildcards).output[i].open() as f:
            table = pd.read_csv(f)
            table_modified = generate_numreads(table, n_reads_dict[wildcards.coverage])
            paths.append(table_modified)
    all_samples_df = pd.concat(paths)
    print("all_samples_df", all_samples_df)
    fq1=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}-{{coverage}}_1.fq",
        zip,
        lineage=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['lineage'],
        num=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['num_reads'],
        coverage=["100x", "1000x"]
        )
    fq2=lambda wc: expand("results/fractions/{{sample}}-{lineage}-{num}-{{coverage}}_2.fq",
        zip,
        lineage=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['lineage'],
        num=all_samples_df.loc[all_samples_df['sample'] == wc.sample]['num_reads'],
        coverage=["100x", "1000x"]
        )
    return [fq1,fq2]

# input function to retrieve fastq samples, (todo: fix this function, more readable)
def get_concat_fractions_input(wildcards):
    if config["simulate_given"] == False:
        fq1=aggregate_input_concat_fractions(wildcards)[0]
        fq2=aggregate_input_concat_fractions(wildcards)[1]
        return [fq1,fq2]
    else:
        fq1=expand("results/fractions/{{sample}}-{lineage}-{num}-{{coverage}}_1.fq",
        zip,
        lineage=simulated_given_lineages.loc[simulated_given_lineages['sample'] == wildcards.sample]['lineage'],
        num=simulated_given_lineages.loc[simulated_given_lineages['sample'] == wildcards.sample]['num_reads'],
        coverage=coverage_single_sample
        )
        fq2=expand("results/fractions/{{sample}}-{lineage}-{num}-{{coverage}}_2.fq",
        zip,
        lineage=simulated_given_lineages.loc[simulated_given_lineages['sample'] == wildcards.sample]['lineage'],
        num=simulated_given_lineages.loc[simulated_given_lineages['sample'] == wildcards.sample]['num_reads'],
        coverage=coverage_single_sample
        )
        return [fq1,fq2]

samples = (
    pd.read_csv(
        config["samples"],
        sep="\t"
    )
)

# input function to retrieve fastq samples
def get_fastq_input(wildcards):
    sample = wildcards.sample
    if config["simulate_pandemics"] or config["simulate_given"]:
        files = ["results/mixed/{sample}_1.fastq", "results/mixed/{sample}_2.fastq"]
    else:
        files = ["results/sra/{sample}_1.fastq.gz", "results/sra/{sample}_2.fastq.gz"]
    return files

#get trimmed fastq input
def get_trimmed_fastq_input(wildcards):
    sample = wildcards.sample
    files = ["results/trimmed/{sample}.1.fastq", "results/trimmed/{sample}.2.fastq"]
    return files

#input function for create_sample_sheet_unicovar
def get_fastq_input_unicovar():
    if config["simulate_pandemics"] and not config["simulate_given"]: #make sure the other is not mistakenly chosen
        fqs1 = expand("results/mixed/SimulatedSample{num}-{coverage}_1.fastq", num=num_list, coverage=["100x","1000x"])
        fqs2 = expand("results/mixed/SimulatedSample{num}-{coverage}_2.fastq", num=num_list, coverage=["100x", "1000x"])
    elif config["simulate_given"] and not config["simulate_pandemics"]: #make sure the other is not mistakenly chosen:
        fqs1 = expand("results/mixed/{sample}-{coverage}_1.fastq",  sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        fqs2 = expand("results/mixed/{sample}-{coverage}_2.fastq",  sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
    else:
        fqs1 = expand("results/trimmed/{sample}.1.fastq.gz", sample=samples["sra"])
        fqs2 = expand("results/trimmed/{sample}.2.fastq.gz", sample=samples["sra"])
        
    return [fqs1, fqs2]
unicovar_inputs = get_fastq_input_unicovar()
print("unicovar_inputs[0]",unicovar_inputs[0])
print("unicovar_inputs[1]", unicovar_inputs[1])

#just a pseudodate used in uncovar workflow, hence pangolin output
DATE="13062024"
UNCOVAR_SAMPLE_SHEET="uncovar/config/pep/samples.csv"
UNCOVAR_CONFIG="uncovar/config/config.yaml"
UNCOVAR_PEP_CONFIG="uncovar/config/pep/config.yaml"

def get_results(wildcards):
    final_output=[]
    if config["simulate_pandemics"] and not config["simulate_given"]: #make sure the other is not mistakenly chosen
        orthanq_csv = expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/SimulatedSample{num}-{coverage}.csv", num=num_list, coverage=["100x", "1000x"])
        orthanq_solutions = expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/viral_solutions.html", num=num_list, coverage=["100x", "1000x"])
        pangolin = expand("results/pangolin/SimulatedSample{num}-{coverage}.csv", num=num_list, coverage=["100x", "1000x"])
        kallisto = expand("results/kallisto/quant_results_SimulatedSample{num}-{coverage}", num=num_list, coverage=["100x", "1000x"])
        nextclade = expand("results/nextstrain/results/SimulatedSample{num}-{coverage}", num=num_list, coverage=["100x", "1000x"])
    elif config["simulate_given"] and not config["simulate_pandemics"]: #make sure the other is not mistakenly chosen:
        orthanq_csv = expand("results/orthanq/calls/{sample}-{coverage}/{sample}-{coverage}.csv", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample) 
        orthanq_solutions = expand("results/orthanq/calls/{sample}-{coverage}/viral_solutions.html", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        pangolin = expand("results/pangolin/{sample}-{coverage}.csv", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        kallisto = expand("results/kallisto/quant_results_{sample}-{coverage}", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        nextclade = expand("results/nextstrain/results/{sample}", sample=simulated_given_lineages["sample"].unique())
    else:
        orthanq_csv = expand("results/orthanq/calls/{sample}/{sample}.csv", sample=samples["sra"])
        orthanq_solutions = expand("results/orthanq/calls/{sample}/viral_solutions.html", sample=samples["sra"])
        pangolin = expand("results/pangolin/{sample}.csv", sample=samples["sra"]) 
        kallisto = expand("results/kallisto/quant_results_{sample}", sample=samples["sra"])
        nextclade = expand("results/nextstrain/results/{sample}", sample=samples["sra"])
    final_output.extend(orthanq_csv + orthanq_solutions + pangolin + kallisto + nextclade)
    return final_output