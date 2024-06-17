import pandas as pd
import re

pseudodate_unicovar="13062024"

unicovar_sheet_template=pd.read_csv(snakemake.input.template)

fq1=list(snakemake.input.fq1)
fq2=list(snakemake.input.fq2)

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
    # print(fq1, fq2)
    unicovar_sample_sheet.loc[index, "fq1"] = "../" + f1
    unicovar_sample_sheet.loc[index, "fq2"] = "../" + f2
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