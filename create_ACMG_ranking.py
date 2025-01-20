# libraries -------------

import pandas as pd
from extract_clinvar_data_for_variant import extract_clinvar_window_around_variant, extract_clinvar_gene_of_variant
import re

# Define Germline classes and rules -------------

# https://varsome.com/about/resources/germline-implementation/#pp5 implementation rules in varsome clinical

# Germline Rules dictionary
Germline_Rules = {
    'Very Strong': 8,
    'Strong': 4,
    'Moderate': 2,
    'Supporting': 1
}

# Functions ----------------

# Function to classify Germline class based on points
def classify_germline_class(points):
    if points >= 10:
        return 'Pathogenic'
    elif 6 <= points <= 9:
        return 'Likely pathogenic'
    elif 0 <= points <= 5:
        return 'VUS'
    elif -6 <= points <= -1:
        return 'LB'
    else:
        return 'B'


# DATA  ----------------

# MOBI data
mobi_data = pd.read_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_omim_splice_clinvarentries.txt", sep='\t')
mobi_data.head(3)
# print(mobi_data['Clinvar_entries'].value_counts())

# PVS1 -----------------------
# nonsense variant
# frameshift variant
# exon deletion variant
# intronic variant within ±2 bases of the transcript
# splice site start loss variant
# LOF is a “Known Mechanism of Disease” from either:
#    -The gene statistics: if at least 2 LOF variants in this gene have been reliably reported as pathogenic.
#    -GnomAD gene constraints LOF Observed/Expected is less than 0.7555.

def calculate_pvs1(df):
    df["PVS1"] = 0  # Initialize PVS1 column

    for index, row in df.iterrows():
        if row['variant_type'] in ['Frameshift', 'Nonsense']:
            df.at[index, "PVS1"] = 8
        # Check for splice site prediction greater than 0.75
        elif row.get('Splice_prediction', 0) > 0.75:
            df.at[index, "PVS1"] = 8
    return df


print(mobi_data['variant_type'].value_counts())
print(mobi_data['Splice_prediction'].value_counts())
mobi_data = calculate_pvs1(mobi_data)
mobi_data["PVS1"].value_counts()

def apply_extract_clinvar_gene(da, clinvar_dir):
    # Initialize an empty list to store results
    clinvar_gene_dfs = []

    # Iterate over each row in da
    for index, row in da.iterrows():
        chrom = row['CHROM']  # Extract chromosome
        gene_name = row['HGNC gene']  # Extract gene name
        clinvar_gene_df = extract_clinvar_gene_of_variant(clinvar_dir, chrom=chrom, gene_name=gene_name)
        clinvar_gene_dfs.append(clinvar_gene_df)
    clinvar_gene_combined_df = pd.concat(clinvar_gene_dfs, ignore_index=True)

    return clinvar_gene_combined_df


def count_LoF_pathogenic(df):
    condition_LoF = df['MC'].isin(
        ['nonsense', 'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant', 'stop_lost'])
    condition_pathogenic = df['CLNSIG'].isin(['Pathogenic', 'Likely_pathogenic'])
    LOF_pathogenic_count = int((condition_LoF & condition_pathogenic).sum())

    return LOF_pathogenic_count


def process_LoF_pathogenic_counts(mobi_data):
    clinvar_dir = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/clinvar_chr/"
    # Filter the rows where PVS1 == 8
    filtered_data = mobi_data[mobi_data['PVS1'] == 8]
    da = filtered_data[['CHROM', 'POS', 'HGNC gene', 'variant_type', 'Splice_prediction', 'PVS1']]

    # Extract clinvar info
    clinvar_gene_combined_df = apply_extract_clinvar_gene(da, clinvar_dir)

    # Initialize the 'LOF_pathogenic_count' column in da
    da['LOF_pathogenic_count'] = ""

    # Loop through each unique gene in the 'HGNC gene' column
    for gene_name in da['HGNC gene'].unique():
        # Filter the DataFrame for the rows containing the current gene in 'GENEINFO'
        filtered_data = clinvar_gene_combined_df[
            clinvar_gene_combined_df['GENEINFO'].str.contains(gene_name, case=False, na=False)]

        # Apply the count_LoF_pathogenic function to the filtered data
        LOF_pathogenic_count = count_LoF_pathogenic(filtered_data)

        # Add the LOF_pathogenic_count to each row that matches the current gene in da
        da.loc[da['HGNC gene'] == gene_name, 'LOF_pathogenic_count'] = LOF_pathogenic_count

    # Update the original mobi_data with the LOF_pathogenic_count column from da
    mobi_data.loc[da.index, 'LOF_pathogenic_count'] = da['LOF_pathogenic_count']

    return mobi_data


mobi_data.shape()
mobi_data = process_LoF_pathogenic_counts(mobi_data)
mobi_data['LOF_pathogenic_count'].value_counts()


# PM2 -------------------
# we do not have the gnomad coverage ...
# For dominant genes (including X-Linked and AD/AR) we check that the allele count is less than 5.
#  For recessive genes (AR): the rule will trigger if the homozygous allele count is less than 2.
#  Alternatively the rule also checks whether the allele frequency is less than 0.0001.

def assign_pm2(df, gnomad_col):
    df["PM2"] = 0
    df[gnomad_col] = pd.to_numeric(df[gnomad_col], errors='coerce')  # Convert to numeric, set invalid parsing as NaN
    df.loc[(df[gnomad_col] < 0.0001), "PM2"] = 2
    df.loc[(df[gnomad_col] == "No match in gnomADv4 Genome"), "PM2"] = 2
    return df

mobi_data = assign_pm2(mobi_data,'gnomAD v4 Exome:')

mobi_data.columns
mobi_data.shape
mobi_data.to_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_omim_splice_clinvarentries.txt", index=False, sep='\t')
# is_common_mechanism_of_disease if mobi_data_updated['LOF_pathogenic_count']>2
# If at least two distinct LOF variants for the gene are classified as pathogenic/likely pathogenic, this suggests a common mechanism of disease.
# TODO: Cross-Check with GnomAD:  Validate this against GnomAD gene constraints by checking the observed/expected LOF ratio (LOF_OE). A value below 0.7555 further supports a strong LOF constraint for the gene.





clinvar_gene_df = extract_clinvar_gene_of_variant(clinvar_dir, chrom=1, gene_name="SAMD11")




# Null variant (frame-shift) in gene GJB2, not predicted to cause NMD. Loss-of-function is a known mechanism of disease (gene has 99 reported pathogenic LOF variants). The truncated region contains 173 pathogenic variants. It removes 75.33% of the protein.
# todo: add NMD - we need to know how many exons in the protein and also how far from end of exon it is if it is the last one
# todo: removes % of protein




clinvar_gene_df = extract_clinvar_gene_of_variant(clinvar_dir, chrom=1, gene_name="SAMD11")








# PS1 -------------------
# Same amino acid change as a previously established pathogenic variant regardless of nucleotide change.
# TODO: Use a mapping tool (e.g., VEP, Ensembl, or similar tools) to annotate the variant at the protein level.
# TODO: include info from Clinvar to check, and also https://www.ncbi.nlm.nih.gov/research/litvar2/
# is there another variant in clinvar in same AA? ( we do not have the AA...

# PS3 -------------------
# Well-established in vitro or in vivo functional studies supportive of a damaging effect on the gene or gene product. (Pathogenic, Strong)
# TODO need to include literature from Clinvar, Pubmed, are there functional studies? id drug response in clinvar column, we can activate PS3
mobi_data
print(mobi_data['Clinvar Germline:'].unique())

# PM1 -------------------
# Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without benign variation. (Pathogenic, Moderate)
# Extract gene info from Clinvar
    # check how many variants in hotspot
        # using a region of 25 base-pairs on either side of the variant, the rule checks that there are at least 4 pathogenic variants (only using missense and inframe-indel variants)
    # functional domain reported by UniProt, at least 2 pathogenic clinically reported missense/in-frame variants in domain


# PM4 -------------------

# Protein length changes as a result of in-frame deletions/insertions in a non-repeat region or stop-loss variants. (Pathogenic, Moderate)
# PM4 applies when in-frame (so change multiple of 3) deletions/insertions or stop-loss variants change the protein length, but not in repeat regions, and does not require NMD.
df["PM4"]=0
if df["PVS1"]!=0, df["PM4"]=0
    # if variant_type = "Frameshift", PVS1, if "in frame", PM4


# PM5 -------------------
# Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before. (Pathogenic, Moderate)
# TODO: Need to sync with clinvar, is there another missense in the same AA?
    # todo need to get gene info in clinvar, or are these variants in another df?

# PP2 -------------------
# Missense variant in a gene that has a low rate of benign missense variation and in which missense variants are a common
# mechanism of disease. (Pathogenic, Supporting)
# In order to avoid double-counting the same evidence, rule PP2 will not be applied if rule PM1 was triggered.
# Look at MPA score, CADD and Alphamissense
df["PP2"]=0

# PP3 -------------------
# Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc.) (Pathogenic, Supporting)# very strong
df["PP3"]=0
# Rule PP3 is disabled if either rule PVS1 or PM4 are triggered, in order to avoid double-counting similar evidence.

# df_simple = df[['HGNC gene', 'HGVS strict genomic (hg38):', 'variant_type', 'HGVS Protein:', 'protein_impact',
#          'HGVS genomic (hg38):', 'pseudo VCF (hg38):', 'Position in transcript:', 'Position / splice site',
#          'Position / domain', 'Position tolerance',
#          'gnomAD exome:', 'gnomAD genome:', 'gnomAD exome (non cancer):',
#          'gnomAD v4 Genome:', 'gnomAD v4 Exome:',
#          'dbSNP rsid:', 'Clinvar Germline:', 'hg38 InterVar:', 'GeneBe:',
#          'CADD phred:', 'MPA score:', 'MPA impact:']]

# PP5 ------------------
# Reputable source recently reports variant as pathogenic, but the evidence is not available to the laboratory to perform an independent evaluation. (Pathogenic, Supporting)
# clinvar pathogenic
# supporting >3 entries, if more than 3 entries:  moderate, if more than 5: strong



# create column pathogenic_variant empty
# TODO: write code in which if phenotypeInheritance_mapped is AD and zyg= 1 and ACMG_score is > 6 : fill yes in pathogenic_variant
# or if phenotypeInheritance_mapped is AR and zyg = 2 and ACMG_score is > 6: fill yes in pathogenic_variant
# careful because sometimes we have "AD; AR" in phenotypeInheritance_mapped

# rank variants by ACMG score in dataframe
# TODO take into account NMD
Variants introducing a premature stop codon may trigger NMD
# if the codon is in an upstream exon and meets the distance criteria of 55 nucleotide
# TODO compute the frequency of protein left
# TODO look at https://varsome.com/about/resources/germline-implementation/

# CADD  >= 25.6: LP, >= 29: P
# AlphaMissense: >=0.8 LP, >=0.95: P

# TODO extract clinvar data








