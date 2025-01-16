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


output_mobidetails = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/output/output_mobidetails'
output_dir_merged = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/merged_data'

iteration_end = 4
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

merged_df.to_csv(f"{output_dir_merged}/merged_data.txt", sep='\t', index=False)
print(f"Merged DataFrame shape: {merged_df.shape}")

# check if some variants are there twice
print(merged_df['HGVS strict genomic (hg38):'].unique().shape)
print(merged_df['HGVS strict genomic (hg38):'].shape)


# TODO: run again the variants with iteration=52 and observation of the driver

# need to extract and understand failed variants
# check that all the variants are present - which ones are missing



