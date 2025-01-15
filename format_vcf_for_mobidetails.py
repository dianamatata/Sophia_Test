import pandas as pd

# Load the VCF file into a pandas DataFrame
vcf_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patient_variants.vcf"
vcf_df = pd.read_csv(vcf_file, sep="\t", comment="#", header=None)

# Extract and format the columns to get the desired format
formatted_entries = vcf_df.apply(lambda row: f"hg38:{row[0]}:{row[1]}:{row[3]}:{row[4]}", axis=1)

# Save the formatted entries to a new file
formatted_entries.to_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/formatted_vcf_entries.txt", index=False, header=False)

