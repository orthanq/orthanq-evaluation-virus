rule get_metadata: 
    output:
        "results/project_metadata/sra_collection_dates.tsv",
    params:
        sra_ids = samples.sra.to_list()
    log:
        "logs/get_metadata/get_metadata.log",
    conda:
        "../envs/entrez.yaml"
    shell:
        "echo -e 'SRA_ID\tCollectionDate' > {output}; "
        "for SRA_ID in {params.sra_ids}; do COLLECTION_DATE=$(esearch -db sra -query $SRA_ID | efetch -format xml | xtract -pattern SAMPLE_ATTRIBUTE -if TAG -equals collection_date -element VALUE); "
        "if [ -z '$COLLECTION_DATE' ]; then COLLECTION_DATE='NA'; fi; "
        "echo -e ${{SRA_ID}}\t${{COLLECTION_DATE}} >> {output}; done"