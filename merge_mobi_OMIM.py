# Libraries ---------------
import pandas as pd

# DATA --------------

# mobidetails
mobi_file = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/merged_mobidetails_vcf.txt'
mobi_data = pd.read_csv(mobi_file, sep='\t')

#OMIM
omim_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/omim.txt"
omim_data = pd.read_csv(omim_file, sep='\t')


# Concatenate OMIM columns ----------------

# Function to concatenate OMIM data -----------

def concatenate_omim_data(group):
    # Ensure all columns are converted to strings before joining
    phenotypeInheritance = '; '.join(map(str, group['phenotypeInheritance_mapped']))
    phenotypes = '; '.join(map(str, group['phenotype']))
    geneMimNumber = '; '.join(map(str, group['geneMimNumber'].unique()))  # Deduplicate and convert to string
    phenotypeMimNumber = '; '.join(map(str, group['phenotypeMimNumber'].unique()))  # Deduplicate and convert to string

    return {
        'genes': group.name,  # Group name (unique 'genes')
        'phenotypeInheritance_mapped': phenotypeInheritance,
        'phenotype': phenotypes,
        'geneMimNumber': geneMimNumber,
        'phenotypeMimNumber': phenotypeMimNumber
    }

# Merge mobi_subset2 with omim2 based on matching gene names
def get_omim_data_for_gene(gene_symbol, omim2):
    # Handle cases where gene_symbol is NaN
    if pd.isna(gene_symbol):
        return None, None, None, None, None

    # Extract the gene symbol (strip '()' to match omim data gene names)
    gene_name = str(gene_symbol).split(' ')[0]  # Ensure gene_symbol is a string
    match = omim2[omim2['genes'].str.contains(rf"\b{gene_name}\b", case=False, na=False)]


    if not match.empty:
        # Return the 'genes' column along with the other relevant data
        return match.iloc[0]['genes'], match.iloc[0]['phenotypeInheritance_mapped'], match.iloc[0]['phenotype'], \
        match.iloc[0]['geneMimNumber'], match.iloc[0]['phenotypeMimNumber']
    else:
        return None, None, None, None, None



# Code to concatenate OMIM data -----------

print(f"merged_data shape: {mobi_data.shape}") # merged_data shape: (2406, 55)

# 1 Group by 'genes' and concatenate the columns into omim2 DataFrame
omim2 = omim_data.groupby('genes').apply(concatenate_omim_data).apply(pd.Series).reset_index(drop=True)
omim2.to_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/omim_grouped_by_genes.txt", sep='\t', index=False)

# 2 Apply the function to each gene in mobi_subset2

# test on small mobi
# mobi_subset2 = mobi_data[['HGNC gene symbol (ID):', 'HGVS strict genomic (hg38):']]
# mobi_subset2[['genes', 'phenotypeInheritance_mapped', 'phenotype', 'geneMimNumber', 'phenotypeMimNumber']] = mobi_subset2['HGNC gene symbol (ID):'].apply(
#     lambda x: pd.Series(get_omim_data_for_gene(x, omim2)))

omim2 = pd.read_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/omim_grouped_by_genes.txt", sep='\t')

mobi_data[['genes', 'phenotypeInheritance_mapped', 'phenotype', 'geneMimNumber', 'phenotypeMimNumber']] = mobi_data['HGNC gene'].apply(
    lambda x: pd.Series(get_omim_data_for_gene(x, omim2)))
print(f"mobi_data shape: {mobi_data.shape}")


# Save -------------------

mobi_data.to_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_with_omim_genes.txt", sep='\t', index=False)
mobi_data = pd.read_csv("/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/mobi_data_with_omim_genes.txt", sep='\t')


# DEBUG
# term="DPYD"
# mobi_subset = mobi_data[mobi_data['HGNC gene symbol (ID):'].str.contains(term)]
# omim_subset = omim_data[omim_data['genes'].str.contains(term)]
# print(mobi_subset[['HGNC gene symbol (ID):', 'HGVS strict genomic (hg38):']])
# print(omim_subset)
# omim2 = omim_subset.groupby('genes').apply(concatenate_omim_data).apply(pd.Series).reset_index(drop=True)

