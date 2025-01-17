# Clearing All Variables ---------------

protected = ['os', '__name__', '__file__', '__doc__', '__package__', '__loader__', '__spec__', '__builtins__']
# Remove all user-defined variables ------------
for name in list(globals().keys()):
    if name not in protected:
        del globals()[name]

# Libraries ---------------

import pandas as pd
pd.set_option('display.max_rows', 5)  # Show all rows
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.width', None)  # Don't truncate width
pd.set_option('display.max_colwidth', None)  # Don't truncate column contents

# Get MOBI Data -------------

# mobidetails
mobi_file = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/merged_data_mobidetails.txt'
mobi_data = pd.read_csv(mobi_file, sep='\t')

mobi_data['hg38_vcf_mobidetails'] = mobi_data['pseudo VCF (hg38):'].apply(
    lambda x: 'hg38:' + x.replace('-', ':'))
# mobi_data.rename(columns={'VCF_formatted': 'hg38_vcf_mobidetails'}, inplace=True)
# mobi_data[['HGNC gene symbol (ID):', 'hg38_vcf_mobidetails']]

# Get VCF Data -------------

vcf_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/patient_variants_vcf_df.txt"
vcf_data = pd.read_csv(vcf_file, sep='\t')

# MERGE ----------
merged_data = pd.merge(vcf_data, mobi_data,
                       how='left',
                       left_on='hg38_vcf_mobidetails',
                       right_on='hg38_vcf_mobidetails')


# CHECKS --------------------------------

print(f"VCF data shape: {vcf_data.shape}")
print(f"mobi_data shape: {mobi_data.shape}")
print(f"merged_data shape: {merged_data.shape}")

# VCF data shape: (2406, 11)
# mobi_data shape: (2379, 44)
# merged_data shape: (2406, 54)

print(f"Missing variants: {vcf_data.shape[0] - mobi_data.shape[0]} ")
# Missing variants: 27
print(merged_data.columns)

failed_variants_log = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/out_mobi/failed_variants_log3.txt"
failed_variants_log = pd.read_csv(failed_variants_log, sep='\t')
print(f"failed_variants_log shape: {failed_variants_log.shape}") # failed_variants_log shape: (125, 1)

# TODO extract just the variants from there

# Add a new column 'mobi_data == True/False'
merged_data['mobi_data'] = merged_data['HGNC gene symbol (ID):'].notna()

# Count the number of False values in 'mobi_data'
false_count = (merged_data['mobi_data'] == False).sum()
print(f"Number of False values in 'mobi_data': {false_count}") # 324
pb_df = merged_data[merged_data['mobi_data'] == False]

# TODO: weird because we have  print(f"failed_variants_log shape: {failed_variants_log.shape}") # failed_variants_log shape: (125, 1)
# print(f"Missing variants: {line_count - merged_df['HGVS strict genomic (hg38):'].unique().shape[0]} ")
# # Missing variants: 27
# in merge_output_mobi.py
# print(f"Number of False values in 'mobi_data': {false_count}") # 324

# Identify rows in vcf_data that are not in mobi_data
difference = vcf_data[~vcf_data['hg38_vcf_mobidetails'].isin(mobi_data['hg38_vcf_mobidetails'])]
print(difference) # [324 rows x 11 columns]

# TODO understand why diff is 324 and before 27 .....

# SAVE DATA

# Save the 'hg38_vcf_mobidetails' column from pb_df --------
pb_df['hg38_vcf_mobidetails'].to_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/problem_vcf_variants.txt", index=False, header=False, sep='\n')

# Save merged_data
merged_data.to_csv('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/merged_mobidetails_vcf.txt', index=False, sep='\t')
output_file = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/merged_mobidetails_vcf.txt'











