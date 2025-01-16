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

# Get Data -------------

# vcf
vcf_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/patient_variants_vcf_df.txt"
vcf_data = pd.read_csv(vcf_file, sep='\t')

# all the variants entries
formatted_entries_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/formatted_vcf_entries.txt"
entries_data = pd.read_csv(formatted_entries_file, sep='\t', header=None)

# variants failed
failed_variants_log = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/out_mobi/failed_variants_log2.txt"

# mobidetails + OMIM
mobi_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_with_omim_genes.txt"
mobi_data = pd.read_csv(mobi_file, sep='\t')



# merge VCF mobidata
mobi_data.rename(columns={'VCF_formatted': 'hg38_vcf_mobidetails'}, inplace=True)
# subset
mobi_data[['HGNC gene symbol (ID):', 'hg38_vcf_mobidetails']]

merged_data = pd.merge(vcf_data, mobi_data[['HGNC gene symbol (ID):', 'hg38_vcf_mobidetails']],
                       how='left',
                       left_on='hg38_vcf_mobidetails',
                       right_on='hg38_vcf_mobidetails')
print(vcf_data.shape) # 2406
print(merged_data.shape) # 2449

# Are there duplicates ?
duplicates_in_merged = merged_data[merged_data.duplicated(subset='hg38_vcf_mobidetails', keep=False)]
print(duplicates_in_merged) # [83 rows x 12 columns]
merged_data = merged_data.drop_duplicates(subset='hg38_vcf_mobidetails', keep='first')

print(merged_data.shape) # 2406




# Add a new column 'mobi_data == True/False'
merged_data['mobi_data'] = merged_data['HGNC gene symbol (ID):'].notna()

# Count the number of False values in 'mobi_data'
false_count = (merged_data['mobi_data'] == False).sum()
# Print the result
print(f"Number of False values in 'mobi_data': {false_count}") # 390
pb_df = merged_data[merged_data['mobi_data'] == False]
# Extract the 'hg38_vcf_mobidetails' column from pb_df
pb_df['hg38_vcf_mobidetails'].to_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/problem_vcf_variants.txt", index=False, header=False, sep='\n')


# Replace NaN in 'VCF_formatted' column with 'NA'
merged_data['VCF_formatted'].fillna('NA', inplace=True)

# Display the resulting merged DataFrame
print(merged_data)










# Identify rows in merged_data that are not in vcf_data
difference = vcf_data[~vcf_data['hg38_vcf_mobidetails'].isin(merged_data['hg38_vcf_mobidetails'])]

# Display the rows that are different
print(difference)


vcf_data['hg38_vcf_mobidetails']
mobi_data['VCF_formatted']
merge these 2 dataframes according to this column that should be identical,
do not remove any row from vcf_data
if data missing in mobi_data, fill with NA and add a column mobi_data == False, otherwise ==True



# Perform a left merge to keep all rows from vcf_data and fill missing values from mobi_data
merged_data = pd.merge(vcf_data, mobi_data[['VCF_formatted', 'hg38_vcf_mobidetails']],
                       how='left',
                       left_on='hg38_vcf_mobidetails',
                       right_on='hg38_vcf_mobidetails')

# Add a new column 'mobi_data == True/False'
merged_data['mobi_data == True'] = merged_data['VCF_formatted'].notna()

# Replace NaN in 'VCF_formatted' column with 'NA'
merged_data['VCF_formatted'].fillna('NA', inplace=True)

# Display the resulting merged DataFrame
print(merged_data)


# Clinvar

# ACMG

# MPA score

# Concatenate OMIM columns ----------------













