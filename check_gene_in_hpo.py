import pandas as pd

# HPO ------------
# check if gene in hpo_data, get associated HPO_ID,Term_Name
def get_hpo_terms_from_gene(gene,hpo_data):

    result = hpo_data[hpo_data['Gene_Associations'].str.contains(gene, case=False, na=False)]

    # Extract relevant columns if gene is found
    if not result.empty:
        associated_hpo = result[['HPO_ID', 'Term_Name']]
        print(associated_hpo)
        return associated_hpo

    else:
        print(f"Gene {gene} not found in any HPO terms.")
        return pd.DataFrame(columns=['HPO_ID', 'Term_Name'])

gene = 'PSAP'
hpo_data = pd.read_csv('/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/hpo_terms_and_gene_associations.csv')
associated_hpo = get_hpo_terms_from_gene(gene, hpo_data)
associated_hpo_dict = associated_hpo.set_index('HPO_ID')['Term_Name'].to_dict()
print(associated_hpo_dict)