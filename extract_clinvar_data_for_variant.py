# libraries -------------

import gzip
import pandas as pd

# functions -------------

def extract_first_10_lines(input_file):
    # Open the compressed file in read mode
    with gzip.open(input_file, 'rt') as f:
        count = 0
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue

            # Print the line and increment the count
            print(line.strip())
            count += 1

            # Stop after 10 lines
            if count == 10:
                break

def clinvar_variant_context(input_file, chrom, pos, window=50):
    # Initialize variables to track the matching line and surrounding lines
    lines = []
    match_index = -1

    # Open the compressed file in read mode
    with gzip.open(input_file, 'rt') as f:
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue

            # Collect all data lines for indexing
            lines.append(line.strip())
            columns = line.strip().split('\t') # Split the line by tabs
            chrom_file, pos_file = columns[0], int(columns[1])

            # Check for a match
            if chrom_file == str(chrom) and pos_file == pos:
                match_index = len(lines) - 1  # Save the line index of the match

    # If a match is found, extract the surrounding lines
    if match_index != -1:
        start = max(0, match_index - window)  # Ensure no negative index
        end = min(len(lines), match_index + window + 1)  # Ensure within bounds
        surrounding_lines = lines[start:end]
        df = pd.DataFrame(surrounding_lines, columns=["Line Content"])

        df = pd.DataFrame(
            [line.split('\t', maxsplit=5) for line in surrounding_lines],
            columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO']
        )

        return df

    # If no match is found, return an empty DataFrame
    return pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO'])

def clinvar_variant_gene(clinvar_file, gene_name, keys_to_extract=None):

    if keys_to_extract is None:
        keys_to_extract = ['CLNVC', 'CLNSIG', 'MC', 'GENEINFO']

    with gzip.open(clinvar_file, 'rt') as infile:
        lines = []
        for line in infile:
            if line.startswith("##"):  # Skip header lines
                continue

            # Check if the line contains the desired gene
            if f"GENEINFO={gene_name}" in line:
                # Collect all data lines for indexing
                lines.append(line.strip())
    df = pd.DataFrame(
        [line.split('\t', maxsplit=5) for line in lines],
        columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO']
    )
    return df

def parse_info(info_str, keys_to_extract=None):
    """
    Parse the INFO field from a VCF line and extract specific keys.

    Parameters:
    - info_str (str): The INFO field as a semicolon-delimited string.
    - keys_to_extract (list): List of keys to extract (default: None, which means extract all).

    Returns:
    - dict: A dictionary with the extracted key-value pairs.
    """
    info_dict = {}
    if keys_to_extract is None:
        keys_to_extract = ['CLNVC', 'CLNSIG', 'MC', 'GENEINFO']

    # Split the INFO field by semicolon
    info_parts = info_str.split(';')
    for part in info_parts:
        if '=' in part:  # Ensure it's a key-value pair
            key, value = part.split('=', 1)  # Split into key and value
            if key in keys_to_extract:  # Only extract specified keys
                info_dict[key] = value

    return info_dict

# info_str = ".\\t.\\tALLELEID=2470506;CLNDISDB=MeSH:D030342,MedGen:C0950123;CLNDN=Inborn_genetic_diseases;CLNHGVS=NC_000004.12:g.186708110T>C;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=FAT1:2195;MC=SO:0001583|missense_variant;ORIGIN=1\\t2470506\\tMeSH:D030342,MedGen:C0950123\\tNC_000004.12:g.186708110T>C\\tUncertain_significance\\tsingle_nucleotide_variant\\tSO:0001483\\tFAT1:2195\\tSO:0001583|missense_variant"
# keys_to_extract = ['CLNVC', 'CLNSIG', 'MC']

def format_data_after_extraction(df, keys_to_extract=None):
    if keys_to_extract is None:
        keys_to_extract = ['CLNVC', 'CLNSIG', 'MC', 'GENEINFO']

    # Parse 'INFO' column into individual components
    parsed_info = df['INFO'].dropna().apply(lambda x: parse_info(x, keys_to_extract))

    # Convert parsed_info into a DataFrame and combine with original DataFrame
    parsed_info_df = pd.DataFrame(parsed_info.tolist())
    df = pd.concat([df.reset_index(drop=True), parsed_info_df], axis=1)

    # Ensure 'MC' and 'GENEINFO' columns are strings before applying .str.split()
    df['MC'] = df['MC'].astype(str).str.split('|').str[1].str.split(',').str[0]
    df['GENEINFO'] = df['GENEINFO'].astype(str).str.split(':').str[0]

    # Select and return relevant columns
    df = df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'CLNVC', 'CLNSIG', 'MC', 'GENEINFO']]

    return df


def extract_clinvar_window_around_variant(clinvar_dir, chrom, pos, window=None, keys_to_extract=None):
    if window is None:
        window = 50
    if keys_to_extract is None:
        keys_to_extract = ['CLNVC', 'CLNSIG', 'MC', 'GENEINFO']

    clinvar_file = f"{clinvar_dir}clinvar_{chrom}.txt.gz"

    # Extract variant context
    clinvar_subset_data = clinvar_variant_context(clinvar_file, chrom, pos, window)

    # Check if the dataframe is empty
    if clinvar_subset_data.empty:
        # print("Warning: No data found for the given variant in the file.")
        return pd.DataFrame()  # Return an empty dataframe or handle accordingly

    # Convert data to string types before applying further processing
    clinvar_subset_data = clinvar_subset_data.apply(lambda x: x.astype(str) if x.dtype != 'O' else x)

    # Format the data after extraction
    clinvar_subset_data = format_data_after_extraction(clinvar_subset_data)

    return clinvar_subset_data


def extract_clinvar_gene_of_variant(clinvar_dir, chrom, gene_name=None, keys_to_extract =None):

    if keys_to_extract is None:
        keys_to_extract = ['CLNVC', 'CLNSIG', 'MC', 'GENEINFO']

    clinvar_file = f"{clinvar_dir}clinvar_{chrom}.txt.gz"
    clinvar_subset_data = clinvar_variant_gene(clinvar_file, gene_name)
    if not clinvar_subset_data.empty:
        clinvar_subset_data = format_data_after_extraction(clinvar_subset_data)

    return clinvar_subset_data



# MAIN ------------------------------

clinvar_dir = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/clinvar_chr/"

clinvar_window_df = extract_clinvar_window_around_variant(clinvar_dir, chrom=4, pos=186709159, window=50)
clinvar_gene_df = extract_clinvar_gene_of_variant(clinvar_dir, chrom=1, gene_name="SAMD11")


# try to handle duplicates in chrom and pos
duplicates = clinvar_window_df[clinvar_window_df.duplicated(subset=['CHROM', 'POS'], keep=False)]
print(duplicates)


# Observe all the possibilities -------------------

a = extract_clinvar_window_around_variant(clinvar_dir, chrom=4, pos=186709159, window=50000000000)

print(a['CLNVC'].unique())
# ['single_nucleotide_variant' 'Microsatellite' 'Indel' 'Deletion'
#  'Duplication' 'Insertion' 'Variation', 'Inversion']
print(a['CLNSIG'].unique())
# ['Likely_benign' 'Uncertain_significance' 'not_provided' 'Benign'
#  'Pathogenic/Likely_pathogenic' 'Likely_pathogenic' 'Pathogenic'
#  'Benign/Likely_benign' 'Conflicting_interpretations_of_pathogenicity'
#  'Benign/Likely_benign|other'
#  'Conflicting_interpretations_of_pathogenicity|other' nan 'association'
#  'risk_factor' 'other' 'drug_response'
#  'Uncertain_significance/Uncertain_risk_allele' 'Likely_risk_allele'
#  'Uncertain_risk_allele' 'Likely_pathogenic/Likely_risk_allele'
#  'Uncertain_risk_allele|protective' 'protective' 'association_not_found'
#  'Affects' 'Affects|association' 'Uncertain_significance|association'
#  'Pathogenic|other']

print(a['MC'].unique())
# ['missense_variant', 'frameshift_variant', 'intron_variant',
#        'synonymous_variant', '5_prime_UTR_variant', nan,
#        'initiatior_codon_variant', 'non-coding_transcript_variant',
#        'splice_acceptor_variant', 'nonsense', 'splice_donor_variant',
#        'genic_upstream_transcript_variant', 'inframe_deletion',
#        'inframe_insertion', 'no_sequence_alteration', 'stop_lost',
#        '3_prime_UTR_variant', 'inframe_indel']


# COUNT the different cases ------------------

# TODO check all the possibilities analyzing 1 chr for  print(clinvar_gene_data['CLNVC'].unique()), we have not all the combinations

df = clinvar_window_df

# Define the conditions

condition_missense = df['MC'] == 'missense_variant'
condition_LoF = df['MC'].isin(['nonsense','frameshift_variant','splice_donor_variant','splice_acceptor_variant','stop_lost'])
condition_inframe = df['MC'].isin(['inframe_insertion','inframe_deletion'])

condition_pathogenic = df['CLNSIG'].isin(['Pathogenic', 'Likely_pathogenic']) # 'Conflicting_interpretations_of_pathogenicity'?
condition_benign = df['CLNSIG'].isin(['Likely_benign' 'Benign']) # 'Conflicting_interpretations_of_pathogenicity'?

# Count
missense_count = condition_missense.sum()
missense_pathogenic_count = (condition_missense & condition_pathogenic).sum()
pathogenic_count = condition_pathogenic.sum()

LOF_pathogenic_count = (condition_LoF & condition_pathogenic).sum()


# Output the results
print(f'Number of missense variants: {missense_count}')
print(f'Number of missense variants pathogenic: {missense_pathogenic_count}')
print(f'Number of pathogenic or likely pathogenic variants: {pathogenic_count}')

# TODO: for a given variant, we check whether it is a hotspot and also in the gene





