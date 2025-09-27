#deleted manually the entry for HXB2 from the accession list.
# the accessions below are from the 5-virus-mix paper, given in the supplementary.
#now, grep -P "U39362|K03455|M93258|AF324493|M38429" results/hiv_genomes/hiv_combined.fasta results in
# >M93258.1 Human immunodeficiency virus type 1 (YU-2) RNA
# >U39362.2 Human immunodeficiency virus type 1 proviral DNA, gag, pol, vif, vpr, vpu, env, tat, rev and nef genes, complete cds
# so, M93258.1 -> YU2 and U39362.2 -> 896, we only need to add the remaining three strains K03455, AF324493 and M38429 to the final collection of fasta sequences.
REM_ACCESSIONS = ["K03455", "AF324493", "M38429"]

rule download_ref_genomes:
    output:
        "results/hiv_genomes/{accession}.fasta"
    log:
        "logs/download_{accession}.log",
    conda:
        "../envs/entrez.yaml"
    params:
        accession = lambda wildcards: wildcards.accession
    benchmark:    
        "benchmarks/ncbi_datasets/download_{accession}.tsv" 
    shell:
        # "esearch -db nucleotide -query {params} | efetch -format fasta > {output} 2> {log}"
        "efetch -db nucleotide -id {params.accession} -format fasta > {output} 2> {log}"

rule concatenate_genomes:
    input:
        expand("results/hiv_genomes/{accession}.fasta", accession=ACCESSIONS)
    output:
        "results/hiv_genomes/hiv_combined.fasta"
    log: 
        "logs/concatenate_genomes.log"
    shell:
        "cat {input} > {output} 2> {log}"

# add the remaining refs
rule download_remaining_ref_genomes:
    output:
        "results/hiv_remaining_genomes/{accession}.fasta"
    log:
        "logs/download_{accession}.log",
    conda:
        "../envs/entrez.yaml"
    params:
        accession = lambda wildcards: wildcards.accession
    benchmark:    
        "benchmarks/ncbi_datasets/download_{accession}.tsv" 
    shell:
        "efetch -db nucleotide -id {params.accession} -format fasta > {output} 2> {log}"

rule concat_remaining:
    input:
        expand("results/hiv_remaining_genomes/{accession}.fasta", accession=REM_ACCESSIONS)
    output:
        "results/hiv_remaining_genomes/3_lineages.fasta"
    log: 
        "logs/concatenate_genomes_ref.log"
    shell:
        "cat {input} > {output} 2> {log}"    

rule add_remaining_refs:
    input:
        all_hiv_genomes="results/hiv_genomes/hiv_combined.fasta",
        remaining_refs="results/hiv_remaining_genomes/3_lineages.fasta"
    output:
        final_list="results/hiv_genomes/hiv_final.fasta"
    log:
        "logs/add_5_virus_refs.log"
    shell:
        "cat {input.all_hiv_genomes} {input.remaining_refs} > {output} 2> {log}"