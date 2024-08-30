import pandas as pd
import re
import os

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    pseudodate_unicovar="13062024"

    unicovar_sheet_template=pd.read_csv(snakemake.input.template)

    fq1=list(snakemake.input.fq1)
    fq2=list(snakemake.input.fq2)

    def extract_substring(full_string, start_substring):
        # Find the starting position of the start_substring in the full_string
        start_pos = full_string.find(start_substring)
        
        # If the start_substring is found, extract the part of the string from that position onwards
        if start_pos != -1:
            return full_string[start_pos:]
        else:
            # If the start_substring is not found, return an empty string or an appropriate message
            return "Substring not found"

    for f in snakemake.input.fq1:
        print("f", f)
        result = extract_substring(f, "/results")
        print("result", ".." + result)

    #natural sorting (the sample names contains both string and number hence sorted() cannot sort them properly)
    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]

    fq1.sort(key=natural_keys)
    fq2.sort(key=natural_keys)

    print(unicovar_sheet_template) 

    #copy
    unicovar_sample_sheet = pd.DataFrame(columns=unicovar_sheet_template.columns)
    print(unicovar_sample_sheet)

    #loop over fastqs and fillin unicovar_sample_sheet
    for index, (f1, f2) in enumerate(zip(fq1, fq2)):
        results1 = extract_substring(f1, "/results")
        results2 = extract_substring(f2, "/results")
        print(results1)
        unicovar_sample_sheet.loc[index, "fq1"] = ".." + results1
        unicovar_sample_sheet.loc[index, "fq2"] = ".." + results2
        unicovar_sample_sheet.loc[index, "date"] = pseudodate_unicovar
        unicovar_sample_sheet.loc[index, "is_amplicon_data"] = 1
        unicovar_sample_sheet.loc[index, "technology"] = "illumina"
        unicovar_sample_sheet.loc[index, "include_in_high_genome_summary"] = 1 #?

        #sample_name
        splitted=f1.split("/")
        sample_name=splitted[-1].split("_")[0]
        unicovar_sample_sheet.loc[index, "sample_name"] = sample_name

    print("unicovar_sample_sheet",unicovar_sample_sheet)
    unicovar_sample_sheet.to_csv(snakemake.output.sample_sheet, index=False)