import pandas as pd

columns_to_keep = [
    'HGNC gene symbol (ID):', 'HGVS strict genomic (hg38):', 'HGVS Protein:',
    'HGVS genomic (hg38):', 'pseudo VCF (hg38):', 'Position in transcript:', 'Position / splice site',
    'Position / domain', 'Position tolerance',
    'gnomAD exome:', 'gnomAD genome:', 'gnomAD exome (non cancer):',
    'gnomAD v4 Genome:', 'gnomAD v4 Exome:',
    'dbSNP rsid:', 'Clinvar Germline:', 'hg38 InterVar:', 'GeneBe:',
    'CADD phred:', 'MPA score:', 'MPA impact:',
    'Position / protein','dbscSNV ADA:', 'dbscSNV RF:', 'spliceAI AG:',
    'spliceAI AL:', 'spliceAI DG:', 'spliceAI DL:', 'AbSplice:', 'SIFT:', 'Polyphen 2 HumDiv:',
    'Polyphen 2 HumVar:', 'Fathmm:', 'AlphaMissense:', 'REVEL:', 'ClinPred:', 'Meta SVM:', 'Meta LR:',
    'Mistic:', 'Interpretation', 'Risk', 'LOVD Effect Reported:','LOVD Matches:'
]


output_mobidetails = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/out_mobi/output_mobidetails'
output_dir_merged = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data'

iteration_end = 58
merged_df = pd.DataFrame()

for iteration in range(1, iteration_end + 1):
    # Load the current CSV file
    file_path = f"{output_mobidetails}_{iteration}.txt"
    df = pd.read_csv(file_path)
    # Ensure we only keep the desired columns that are present in the current DataFrame
    columns_to_keep_existing = [col for col in columns_to_keep if col in df.columns]
    df = df[columns_to_keep_existing]

    # Append the current DataFrame to the merged DataFrame
    merged_df = pd.concat([merged_df, df], axis=0, ignore_index=True, sort=False)
    print(merged_df.shape)
    print(iteration)

print(merged_df['HGVS strict genomic (hg38):'].shape) # (2605,)

# Drop rows where all columns are duplicated
merged_df_no_duplicates = merged_df.drop_duplicates()
print(f"Original shape: {merged_df.shape}")
print(f"New shape after removing duplicates: {merged_df_no_duplicates.shape}")

# View all rows that are duplicates (both first and subsequent occurrences)
all_duplicates = merged_df[merged_df.duplicated(subset='HGVS strict genomic (hg38):',keep='first')]

# Drop rows where 'HGVS strict genomic (hg38):' are duplicated
merged_df_no_duplicates2 = merged_df.drop_duplicates(subset='HGVS strict genomic (hg38):', keep='first')
print(f"Number of unique 'HGVS strict genomic (hg38):' : {merged_df_no_duplicates2.shape[0]}")

# SAVE DATA
merged_df_no_duplicates2.to_csv(f"{output_dir_merged}/merged_data_mobidetails.txt", sep='\t', index=False)
print(f"Merged DataFrame shape: {merged_df.shape}")



# Count the variants/lines in the initial file formatted_entries_file
formatted_entries_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/formatted_vcf_entries.txt"
with open(formatted_entries_file, 'r') as file:
    line_count = sum(1 for line in file)

print(f"formatted_entries_file contains {line_count} lines.")
# formatted_entries_file contains 2406 lines.
print(f"Missing variants: {line_count - merged_df['HGVS strict genomic (hg38):'].unique().shape[0]} ")
# Missing variants: 27

# TODO need to extract and understand failed variants
# TODO check that all the variants are present - which ones are missing



