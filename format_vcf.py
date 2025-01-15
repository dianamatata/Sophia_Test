# Libraries ---------------
import pandas as pd

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

    vcf_data = vcf_data[['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AD', 'freq_ref_allele', 'DP4']]
    return vcf_data

# Main ---------------

# Define file paths
vcf_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patient_variants.vcf"
vcf_output_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patient_variants_AF.csv"

# Load VCF
vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
vcf_data.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
print(vcf_data.head())
print(f"Number of rows in VCF data: {len(vcf_data)}")
vcf_data = process_vcf(vcf_data)

# DP=669: Total of 669 reads at this position.
# AD=506,163: Of those, 506 reads support the reference allele (C), and 163 reads support the alternate allele (A).
# DP4=506,0,163,0: 506 reads supporting the reference allele are from the forward strand, and 163 reads supporting the alternate allele are from the forward strand as well (none from the reverse strand).

vcf_data.to_csv(vcf_output_file, sep='\t', index=False)
vcf_data = pd.read_csv(vcf_output_file, sep='\t')

