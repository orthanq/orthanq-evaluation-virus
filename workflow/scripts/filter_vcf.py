import pysam

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
   
    input_vcf = snakemake.input.vcf
    output_vcf = snakemake.output.filtered_vcf

    # Open the input VCF file
    vcf_in = pysam.VariantFile(input_vcf, "r")

    # Get sample names from the header
    samples = list(vcf_in.header.samples)

    # Create dictionaries to count `1|1` occurrences and total variants for each sample
    sample_counts = {sample: 0 for sample in samples}  # Count of `1|1` genotypes per sample
    sample_totals = {sample: 0 for sample in samples}  # Total number of variants per sample

    # Step 1: Count `1|1` occurrences and total variants for each sample
    for record in vcf_in:
        for sample in samples:
            gt = record.samples[sample]["GT"]
            if gt is not None:
                gt_str = "|".join(map(str, gt))  # Convert GT tuple to "0|1" format
                sample_totals[sample] += 1
                if gt_str == "1|1":
                    sample_counts[sample] += 1

    # Step 2: Determine which samples meet the 30% `1|1` threshold
    samples_to_keep = [
        sample for sample in samples if sample_totals[sample] > 0 and 
        (sample_counts[sample] / sample_totals[sample]) >= 0.3
    ]

    print("Samples to keep based on the 30% `1|1` threshold:")
    print(samples_to_keep)

    # Step 3: Create a new header and copy the metadata from the original header
    new_header = pysam.VariantHeader()

    # Add INFO metadata from the original header
    for info_record in vcf_in.header.info.items():
        key, value = info_record
        new_header.add_meta(str(key), str(value))

    # Add FORMAT metadata from the original header
    for format_record in vcf_in.header.formats.items():
        key, value = format_record
        new_header.add_meta(str(key), str(value))

    # Add FILTER metadata from the original header
    for filter_record in vcf_in.header.filters.items():
        key, value = filter_record
        new_header.add_meta(str(key), str(value))

    # Step 4: Add the selected samples to the new header
    for sample in samples_to_keep:
        new_header.add_sample(sample)

    # Step 5: Create the output VCF file with the new header
    vcf_out = pysam.VariantFile(output_vcf, "w", header=new_header)

    # Step 6: Write the variants to the output VCF, keeping only selected samples
    vcf_in.reset()  # Reset the VCF file reader to the beginning

    # Iterate over each variant record
    for record in vcf_in:
        # Create a dictionary for the filtered samples
        filtered_samples = {sample: record.samples[sample] for sample in samples_to_keep if sample in record.samples}

        # Create a new record by copying the original record
        new_record = record.copy()  # Copy the original variant record
        new_record.samples = filtered_samples  # Assign the filtered samples to the new record

        # Write the new record to the output file
        vcf_out.write(new_record)

    # Close the VCF files
    vcf_in.close()
    vcf_out.close()
