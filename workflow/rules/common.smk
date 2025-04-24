import pandas as pd

configfile: "config/config.yaml"

n_reads_dict = {"100x": config["n_reads_100x"],"1000x": config["n_reads_1000x"]}
print(n_reads_dict)

num_samples = config["number_of_samples"]

num_list = [num+1 for num in range(num_samples)]

# #just a pseudodate used in uncovar workflow, hence pangolin output
DATE="13062024"

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

def get_primers(wildcards):
    primer_version = samples[samples['sra'] == wildcards.sample]['primer'].values[0]
    return f"resources/primer_{primer_version}.bed"

# input function to retrieve fastq samples
def get_raw_fastq_input(wildcards):
    sample = wildcards.sample
    if config["simulate_pandemics"] or config["simulate_given"]:
        files = ["results/mixed/{sample}_1.fastq", "results/mixed/{sample}_2.fastq"]
    else:
        files = ["results/sra/{sample}_1.fastq.gz", "results/sra/{sample}_2.fastq.gz"]
    return files

#get trimmed fastq input
def get_primer_trimmed_fastq_input(wildcards):
    sample = wildcards.sample
    files = ["results/extracted_reads/{sample}.forward.fastq", "results/extracted_reads/{sample}.reverse.fastq"]
    return files

#get adapter trimmed fastq input
def get_adapter_trimmed_fastq_input(wildcards):
    sample = wildcards.sample
    files = ["results/fastp-trimmed/pe/{sample}.1.fastq", "results/fastp-trimmed/pe/{sample}.2.fastq"]
    return files

#do primer trimming if primer column exists, if not, just do adapter trimming
def get_processed_fastq_input(wildcards):
    """Dynamically determine the input FASTQ based on whether a primer is present for the sample."""
    sample = wildcards.sample  # Get sample name from wildcards
    if sample in samples["sra"].values:  # Ensure sample exists
        row = samples[samples["sra"] == sample]  # Get the row for the sample
        if "primer" in row.columns and not row["primer"].isna().values[0]:  # Check if primer column exists (does not exist in lavmix) and has a value
            return get_primer_trimmed_fastq_input(wildcards)  # Use primer-trimmed
    return get_adapter_trimmed_fastq_input(wildcards)  # Use adapter-trimmed by default

def fastp_params():
    if "labmix" in config["samples"]:
        return "--detect_adapter_for_pe --phred64 --qualified_quality_phred 30 --cut_by_quality3 --cut_by_quality5"
    else:
        return "--detect_adapter_for_pe"

def get_tool_outputs(wildcards):
    final_output=[]
    if config["simulate_pandemics"] and not config["simulate_given"]: #make sure the other is not mistakenly chosen
        orthanq_csv = expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/SimulatedSample{num}-{coverage}.csv", num=num_list, coverage=["100x", "1000x"])
        orthanq_solutions = expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/viral_solutions.html", num=num_list, coverage=["100x", "1000x"])
        orthanq_lp_datavzrd_tsv = expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/lp_solution.tsv", num=num_list, coverage=["100x", "1000x"])
        orthanq_lp_datavzrd_report = expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/datavzrd_report", num=num_list, coverage=["100x", "1000x"])
        orthanq_final_solution=expand("results/orthanq/calls/SimulatedSample{num}-{coverage}/final_solution.html", num=num_list, coverage=["100x", "1000x"])

        pangolin = expand("results/pangolin/SimulatedSample{num}-{coverage}.csv", num=num_list, coverage=["100x", "1000x"])
        kallisto = expand("results/kallisto/quant_results_SimulatedSample{num}-{coverage}", num=num_list, coverage=["100x", "1000x"])
        nextclade = expand("results/nextstrain/results/SimulatedSample{num}-{coverage}", num=num_list, coverage=["100x", "1000x"])
    elif config["simulate_given"] and not config["simulate_pandemics"]: #make sure the other is not mistakenly chosen:
        orthanq_csv = expand("results/orthanq/calls/{sample}-{coverage}/{sample}-{coverage}.csv", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample) 
        orthanq_solutions = expand("results/orthanq/calls/{sample}-{coverage}/viral_solutions.html", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        orthanq_lp_datavzrd_tsv = expand("results/orthanq/calls/{sample}-{coverage}/lp_solution.tsv", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        orthanq_lp_datavzrd_report = expand("results/orthanq/calls/{sample}-{coverage}/datavzrd_report", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        orthanq_final_solution=expand("results/orthanq/calls/{sample}-{coverage}/final_solution.html", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)

        pangolin = expand("results/pangolin/{sample}-{coverage}.csv", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        kallisto = expand("results/kallisto/quant_results_{sample}-{coverage}", sample=simulated_given_lineages["sample"].unique(), coverage=coverage_single_sample)
        nextclade = expand("results/nextstrain/results/{sample}", sample=simulated_given_lineages["sample"].unique())
    else:
        orthanq_csv = expand("results/orthanq/calls/{sample}/{sample}.csv", sample=samples["sra"])
        orthanq_solutions = expand("results/orthanq/calls/{sample}/viral_solutions.html", sample=samples["sra"])
        orthanq_lp_datavzrd_tsv = expand("results/orthanq/calls/{sample}/lp_solution.tsv", sample=samples["sra"])
        orthanq_lp_datavzrd_report = expand("results/orthanq/calls/{sample}/datavzrd_report", sample=samples["sra"])
        orthanq_final_solution=expand("results/orthanq/calls/{sample}/final_solution.html", sample=samples["sra"])

        pangolin = expand("results/pangolin/{sample}.csv", sample=samples["sra"]) 
        kallisto = expand("results/kallisto/quant_results_{sample}", sample=samples["sra"])
        nextclade = expand("results/nextstrain/results/{sample}", sample=samples["sra"])
    #orthanq datavzrd creation takes too long and it was disabled temporarily
    final_output.extend(orthanq_csv + orthanq_solutions + orthanq_final_solution + orthanq_lp_datavzrd_tsv + orthanq_lp_datavzrd_report + kallisto + pangolin + nextclade)
    final_output.extend(orthanq_csv + orthanq_solutions + orthanq_final_solution + kallisto + pangolin + nextclade)
    return final_output

def get_results_real_data(wildcards):
    tool_outputs = list(get_tool_outputs(wildcards))
    
    if "labmix" in config["samples"]:
        scatter_plots = [
            "results/evaluation-hiv/plots/orthanq/scatter_plot.svg",
            "results/evaluation-hiv/plots/orthanq/scatter_plot.html",
            "results/evaluation-hiv/plots/kallisto/scatter_plot.svg", 
            "results/evaluation-hiv/plots/kallisto/scatter_plot.html"
        ]
        
        scatter_plots += expand("results/orthanq/calls/{sample}/{sample}.csv", sample=samples["sra"]) 
        scatter_plots += expand("results/orthanq/calls/{sample}/viral_solutions.html", sample=samples["sra"])
        # scatter_plots += expand("results/orthanq/calls/{sample}/lp_solution.tsv", sample=samples["sra"])
        # scatter_plots += expand("results/orthanq/calls/{sample}/datavzrd_report", sample=samples["sra"])
        scatter_plots += expand("results/orthanq/calls/{sample}/final_solution.html", sample=samples["sra"])

        return scatter_plots

    else:
        return tool_outputs + ["results/evaluation-real-data/plots/stacked_barchart.svg",
        "results/evaluation-real-data/plots/stacked_barchart.html",
        "results/evaluation-real-data/tables/all_tools_predictions.csv",
        "results/evaluation-real-data/datavzrd-report/all_tool_predictions"]
        

def get_viral_lineages_path():
    if "labmix" in config["samples"]:
        return "results/kallisto_index/hiv_viral_lineages.idx"
    else:
        return "results/kallisto_index/sarscov2_viral_lineages.idx"

def get_ref_seq_path():
    if "labmix" in config["samples"]:
        return "results/ref/hiv_reference_sequence.fasta"
    return "results/ref/sarscov2_reference_sequence.fasta"
