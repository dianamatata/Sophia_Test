# libraries -------------

import gzip
import pandas as pd

# functions -------------

def extract_gene_variants_clinvar(input_vcf, output_file, gene_name="SAMD11"):
    """
    Extracts all entries related to the specified gene from a VCF file and writes them to a text file.

    Parameters:
        input_vcf (str): Path to the input VCF file (compressed .gz format).
        output_file (str): Path to the output text file.
        gene_name (str): Gene name to filter (default is 'SAMD11').
    """
    with gzip.open(input_vcf, 'rt') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip header lines
            if line.startswith("##"):
                continue

            # Check if the line contains the desired gene
            if f"GENEINFO={gene_name}" in line:
                print(line)
                outfile.write(line)

# Function to parse the INFO column into a dictionary
def parse_info(info_str):
    info_dict = {}
    info_parts = info_str.split(';')  # Split by semicolon
    for part in info_parts:
        if '=' in part:  # Ensure it's a key-value pair
            key, value = part.split('=', 1)  # Split into key and value
            info_dict[key] = value
    return info_dict

# file paths -------------

# Specify file paths
input_vcf = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/downloaded_data/clinvar_20230527.vcf.gz"
clinvar_gene_data_output_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/SAMD11_variants.txt"

# Extract SAMD11 variants
extract_gene_variants_clinvar(input_vcf, clinvar_gene_data_output_file)
print(f"Extracted SAMD11 variants saved to {clinvar_gene_data_output_file}")

# Read the file without a header and assign column names
column_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
clinvar_gene_data = pd.read_csv(clinvar_gene_data_output_file, sep='\t', header=None, names=column_names)
# Apply the function to the INFO column
info_dicts = clinvar_gene_data['INFO'].apply(parse_info)
# Convert the list of dictionaries into a DataFrame
info_df = pd.DataFrame(list(info_dicts))
# Combine the original DataFrame with the expanded INFO DataFrame
clinvar_gene_data = pd.concat([clinvar_gene_data, info_df], axis=1)
clinvar_gene_data['MC'] = clinvar_gene_data['MC'].str.split('|').str[1]
# clinvar_gene_data.head(5)
# clinvar_gene_data.columns

# visualize columns --------------

print(clinvar_gene_data['CLNVC'].unique())
# ['single_nucleotide_variant' 'Microsatellite' 'Indel' 'Deletion'
#  'Duplication' 'Insertion' 'Variation']
print(clinvar_gene_data['CLNSIG'].unique())
# ['Uncertain_significance' 'Likely_benign' 'Benign'
#  'Conflicting_interpretations_of_pathogenicity' 'Pathogenic']
print(clinvar_gene_data['MC'].unique())
# ['missense_variant' 'synonymous_variant' 'splice_donor_variant'
#  'intron_variant' 'inframe_deletion' 'nonsense' 'frameshift_variant'
#  'inframe_insertion' ,  'no_sequence_alteration']
# TODO  and a few with errors to fix

# save data --------------

# 'GENEINFO': gene;number
# keep only interesting columns
clinvar_gene_data = clinvar_gene_data[['CHROM', 'POS', 'ID', 'REF', 'ALT',
       'ALLELEID', 'CLNDISDB',  'CLNHGVS', 'CLNSIG',
       'CLNVC', 'CLNVCSO', 'GENEINFO', 'MC']]

formatted_entries2 = clinvar_gene_data.apply(lambda row: f"hg38:{row[0]}:{row[1]}:{row[3]}:{row[4]}", axis=1)
# formatted_entries2.head()
clinvar_gene_data['hg38_vcf_clinvar'] = formatted_entries2

clinvar_gene_data.to_csv(clinvar_gene_data_output_file, sep='\t', index=False)
clinvar_gene_data = pd.read_csv(clinvar_gene_data_output_file, sep='\t')




