# libraries -------------
import pandas as pd
import re
import numpy as np  # For NaN


# functions -------------
# From 'HGVS strict genomic (hg38):', id del or dup: look if inframe or frameshift
# From 'HGVS strict genomic (hg38):': check protein meaning: substitution, nonsense, Uncertain

# row = mobi_data.iloc[0] #for debug

def determine_variant_type_and_impact(row):
    variant_type, protein_impact = "Unknown", "Unknown"

    # Safely extract variant and protein impact
    variant = row.get('HGVS strict genomic (hg38):', None)
    protein = row.get('HGVS Protein:', None)

    # Handle NaN or None for variant
    if pd.isna(variant):
        return pd.Series([np.nan, np.nan], index=['variant_type', 'protein_impact'])

    variant = str(variant)  # Ensure it’s treated as a string

    # Variant type
    if "del" in variant or "dup" in variant:
        if 'del' in variant:
            try:
                start, end = variant.split(':')[1].split('del')[0].split('.')[1].split('_')
                start, end = int(start), int(end)
            except (ValueError, IndexError):
                start, end = None, None
        elif 'dup' in variant:
            try:
                coords = variant.split(':')[1].split('dup')[0].split('.')
                start, end = (coords[1].split('_') if '_' in coords[1] else (coords[1], coords[1])) if len(
                    coords) > 1 else (None, None)
            except (ValueError, IndexError):
                start, end = None, None

        if start is not None and end is not None:
            diff_length = int(end) - int(start) + 1
            if diff_length % 3 != 0:
                variant_type = "Frameshift"
            else:
                variant_type = "In-frame"
    elif ">" in variant:
        variant_type = "SNP"

    # Protein impact
    if pd.isna(protein):
        protein_impact = np.nan
    else:
        protein = str(protein)  # Ensure it's treated as a string
        if "Ter" in protein:
            protein_impact = "Nonsense"
        elif "?" in protein:
            protein_impact = "Uncertain"
        elif re.search(r'p\.\([A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)(?!Ter)', protein):
            protein_impact = "Substitution of AA"
        else:
            protein_impact = "Uncertain"

    return pd.Series([variant_type, protein_impact], index=['variant_type', 'protein_impact'])

# extract Distance from splice site
def extract_dist_ss(row):
    match = re.match(r'(\d+)', str(row))  # Match the first number before any space
    if match:
        return int(match.group(1))  # Return the number as an integer
    return None  # Return None if no match is found

def column_to_numeric(df, column_name):
    # Convert the column to numeric, invalid parsing will be set as NaN
    df[column_name] = pd.to_numeric(df[column_name], errors='coerce')
    df[column_name] = df[column_name].fillna('') # Replace NaN with an empty string
    df[column_name] = pd.to_numeric(df[column_name], errors='coerce')
    return df

def add_splice_column(df):
    # columns_of_interest = ['dbscSNV ADA:', 'dbscSNV RF:', 'spliceAI AG:', 'spliceAI AL:', 'spliceAI DG:', 'spliceAI DL:',
    #                        'AbSplice:', 'SIFT:', 'SPiP Interpretation', 'SPiP Risk', "Distance from splice site"]

    # Check the condition and create the 'splice' column
    df['Splice_prediction'] = df.apply(lambda row: 1 if (row['dbscSNV ADA:'] > 0.75 or row['dbscSNV RF:'] > 0.75 or
                                                         row['spliceAI AG:'] > 0.75 or row['spliceAI AL:'] > 0.75
                                                         or row['spliceAI DG:'] > 0.75 or row['spliceAI DL:'] > 0.75)
                                                    else (0.2 if row['Distance from splice site'] < 3 else 0), axis=1)
    return df





# PHASE 2: keep cleaning the data ------------------

mobi_data = pd.read_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_with_omim_genes.txt", sep='\t')

print(mobi_data.shape) #(2406, 60)

# Rename the specified columns
mobi_data = mobi_data.rename(columns={'HGNC gene symbol (ID):': 'HGNC gene',
                                      'Interpretation': 'SPiP Interpretation',
                                      'Risk': 'SPiP Risk'})
mobi_data['HGNC gene'] = mobi_data['HGNC gene'].str.extract(r'([A-Za-z0-9]+)').fillna('')

# Apply functions to data -------------

mobi_data[['variant_type', 'protein_impact']] = mobi_data.apply(determine_variant_type_and_impact, axis=1)
mobi_data['Distance from splice site'] = mobi_data['Position / splice site'].apply(extract_dist_ss)

# Clean the 'dbscSNV ADA:' column
mobi_data = column_to_numeric(mobi_data, 'dbscSNV ADA:')
mobi_data = column_to_numeric(mobi_data, 'dbscSNV RF:')

# SPLICING INFO --------------

for col in ['spliceAI AG:', 'spliceAI AL:', 'spliceAI DG:', 'spliceAI DL:']:
    mobi_data[col] = mobi_data[col].str.extract(r'([0-9\.]+)').astype(float)

# splice column
mobi_data = add_splice_column(mobi_data)

# to only check splicing pred
# mobi_data_splicing_info = add_splice_column(mobi_data[['dbscSNV ADA:', 'dbscSNV RF:', 'spliceAI AG:', 'spliceAI AL:', 'spliceAI DG:',
# #        'spliceAI DL:', 'AbSplice:', 'SPiP Interpretation', 'SPiP Risk','Distance from splice site']])


# Reorder the columns in the desired order
mobi_data = mobi_data[['CHROM', 'POS', 'REF', 'ALT', 'DP', 'AD', 'freq_ref_allele',
       'freq_alt_allele', 'DP4', 'Zygocity', 'hg38_vcf_mobidetails',
       'HGNC gene', 'HGVS strict genomic (hg38):', 'HGVS Protein:',
       'HGVS genomic (hg38):', 'variant_type', 'protein_impact','pseudo VCF (hg38):', 'Position in transcript:',
       'Position / splice site', 'Distance from splice site', 'Position / domain', 'Position tolerance',
       'gnomAD exome:', 'gnomAD genome:', 'gnomAD exome (non cancer):',
       'gnomAD v4 Genome:', 'gnomAD v4 Exome:', 'dbSNP rsid:',
       'Clinvar Germline:', 'hg38 InterVar:', 'GeneBe:', 'CADD phred:',
       'MPA score:', 'MPA impact:', 'Position / protein',
       'dbscSNV ADA:', 'dbscSNV RF:', 'spliceAI AG:', 'spliceAI AL:', 'spliceAI DG:',
       'spliceAI DL:', 'AbSplice:', 'SPiP Interpretation', 'SPiP Risk', 'Splice_prediction',
       'SIFT:', 'Polyphen 2 HumDiv:', 'Polyphen 2 HumVar:', 'Fathmm:', 'AlphaMissense:', 'REVEL:',
       'ClinPred:', 'Meta SVM:', 'Meta LR:', 'Mistic:',
       'LOVD Matches:', 'LOVD Effect Reported:', 'mobi_data',
       'genes', 'phenotypeInheritance_mapped', 'phenotype', 'geneMimNumber', 'phenotypeMimNumber']]


# save data -----------
mobi_data.to_csv('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_omim_splice.txt', index=False, sep='\t')
mobi_data = pd.read_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_omim_splice.txt", sep='\t')
# mobi_data2 = mobi_data


# Clinvar
mobi_data_clinvar_info = mobi_data[['CHROM', 'POS', 'HGNC gene', 'dbSNP rsid:','Clinvar Germline:']]
mobi_data['HGNC gene'] = mobi_data['HGNC gene'].str.extract(r'([A-Za-z0-9]+)').fillna('')

print(mobi_data_clinvar_info['Clinvar Germline:'].unique())

# Count occurrences of each unique value in 'Clinvar Germline:' and sort by count
clinvar_counts = mobi_data_clinvar_info['Clinvar Germline:'].value_counts().reset_index()
clinvar_counts = clinvar_counts.sort_values(by='Clinvar Germline:', ascending=False)
print(clinvar_counts)

# subset_clinvar = mobi_data_clinvar_info[mobi_data_clinvar_info['Clinvar Germline:'].str.contains(r'pathogenic|Pathogenic', case=False, na=False)]

# Filter the rows where 'Clinvar Germline:' contains 'pathogenic' or 'Pathogenic', but exclude conflicting terms
subset_clinvar_pathogenic = mobi_data_clinvar_info[
    mobi_data_clinvar_info['Clinvar Germline:'].str.contains(r'pathogenic|Pathogenic', case=False, na=False) &
    ~mobi_data_clinvar_info['Clinvar Germline:'].str.contains(r'conflicting|benign', case=False, na=False)
]

print(subset_clinvar_pathogenic)

#TODO for all the rows of subset_clinvar_pathogenic, i want to count the number of entries in clinvar




# ACMG

# MPA score

# HPO
hpo_file = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/hpo_terms_and_gene_associations.csv'
# hpo_data = pd.read_csv(hpo_file, sep='\t')
















# not working :
df['criteria'] = df['hg38 InterVar:'].str.replace(r"with the following criteria:\\n\\n.*", "", regex=True)
df['criteria'] = df['criteria'].str.replace(r"\\n", "").str.replace(r"\n", "")
# def extract_ACMG(row):
#     if pd.isna(row):  # Check for NaN
#         return row  # If NaN, return it as is
#     else:
#         row = row.replace("with the following criteria:\\n\\n", "")
#         row = row.replace("\\n", "")
#         print(row)
#         return row
#
# df['criteria'] = df['hg38 InterVar:'].apply(extract_ACMG)
# df['criteria']

# var="Likely pathogenic with the following criteria:\\n\\n  PM1\\n  PM2\\n  PP3\\n  PP5"
# parts = var.split("with the following criteria:\\n\\n", 1)
# result = parts[0] + parts[-1]
# result = "".join(result.split("\\n"))

# TODO: then check the missing rows or VCF,