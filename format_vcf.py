# Libraries ---------------
import pandas as pd
import numpy as np

# Functions ---------------

# Function to split INFO column and extract relevant data
def extract_info(info):
    info_dict = {}
    # Split the INFO column by semicolon
    info_parts = info.split(';')

    for part in info_parts:
        # Split each part by '=' to get the key-value pairs
        key_value = part.split('=')
        if len(key_value) == 2:
            key, value = key_value
            info_dict[key] = value

    return info_dict


# Function to split separate INFO column from dict to new columns
def process_vcf(vcf_data):
    # Create new columns from the extracted data
    info_dicts = vcf_data['INFO'].apply(extract_info)
    vcf_data = pd.concat([vcf_data, pd.json_normalize(info_dicts)], axis=1)

    # Split AD and DP4 into separate columns
    vcf_data[['AD_1', 'AD_2']] = vcf_data['AD'].str.split(',', expand=True)
    vcf_data['AD_1'] = pd.to_numeric(vcf_data['AD_1'], errors='coerce')
    vcf_data['AD_2'] = pd.to_numeric(vcf_data['AD_2'], errors='coerce')
    vcf_data['freq_ref_allele'] = vcf_data['AD_1'] / (vcf_data['AD_1'] + vcf_data['AD_2'])
    vcf_data['freq_alt_allele'] = vcf_data['AD_2'] / (vcf_data['AD_1'] + vcf_data['AD_2'])

    vcf_data = vcf_data[['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AD', 'freq_ref_allele', 'freq_alt_allele','DP4']]
    return vcf_data

def add_zygosity_column(vcf_data):
    # Define the conditions for Zygocity
    conditions = [
        (vcf_data['freq_ref_allele'] >= 0.66),  # Zygocity = 0
        (vcf_data['freq_ref_allele'] >= 0.33) & (vcf_data['freq_ref_allele'] < 0.66),  # Zygocity = 1
        (vcf_data['freq_ref_allele'] <= 0.33)  # Zygocity = 2
    ]

    # Define the corresponding Zygocity values
    values = [0, 1, 2]

    # Create the new 'Zygocity' column using numpy's select function
    vcf_data['Zygocity'] = np.select(conditions, values, default=np.nan)

    return vcf_data

# Main ---------------

# Define file paths
vcf_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patient_variants.vcf"
vcf_output_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/patient_variants_vcf_df.txt"

# Load VCF
vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
vcf_data.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
print(vcf_data.head())
print(f"Number of rows in VCF data: {len(vcf_data)}")
vcf_data = process_vcf(vcf_data)
vcf_data = add_zygosity_column(vcf_data)

# add column with hg38 formatting to merge with the other dataframe with mobidetails
formatted_entries = vcf_data.apply(lambda row: f"hg38:{row[0]}:{row[1]}:{row[2]}:{row[3]}", axis=1)
vcf_data['hg38_vcf_mobidetails'] = formatted_entries

vcf_data.to_csv(vcf_output_file, sep='\t', index=False)
vcf_data = pd.read_csv(vcf_output_file, sep='\t')

formatted_entries.to_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/formatted_vcf_entries.txt", index=False, header=False)
