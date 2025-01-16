# Libraries ---------------
import pandas as pd
import numpy as np

# Functions ---------------

def map_inheritance(value):
    if pd.isna(value):
        return np.nan

    mappings = {
        "Autosomal dominant": "AD",
        "Autosomal recessive": "AR",
        "X-linked dominant": "XLD",
        "X-linked recessive": "XLR",
        "Mitochondrial": "MT",
        "X-linked": "XL",
        "Y-linked": "YL"
    }

    result = []

    # sorted(mappings.items(), key=lambda x: len(x[0]), reverse=True). This ensures that "X-linked dominant" is matched before "X-linked."
    # mappings.items(): This converts the dictionary mappings into a list of key-value pairs (tuples)
    # Sort mappings by key length in descending order to not have "X-linked dominant": "XLD", also mapped as "XL"
    # For the tuple ("X-linked", "XL"), len(x[0]) is 8
    for key, short_form in sorted(mappings.items(), key=lambda x: len(x[0]), reverse=True):
        if key in value and f"?{key}" not in value:
            result.append(short_form)

    return "; ".join(result) if result else np.nan

# value = '?Autosomal dominant; Autosomal recessive'
# key = "Autosomal dominant"


# Main ---------------

omim_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/downloaded_data/full_omim_table.txt"
omim_output_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/omim.txt"

omim_df = pd.read_csv(omim_file, sep='\t')
print(omim_df.columns)
print(omim_df['phenotypeInheritance'].unique())

# Test with unique values to check correctly mapped in dummy df
cols = omim_df['phenotypeInheritance'].unique()
df = pd.DataFrame(cols)
df.rename(columns={df.columns[0]: "phenotypeInheritance"}, inplace=True)

omim_df["phenotypeInheritance_mapped"] = omim_df["phenotypeInheritance"].apply(map_inheritance)
df["phenotypeInheritance_mapped"] = df["phenotypeInheritance"].apply(map_inheritance)

# simplify OMIM:
# remove columns: 'hgnc_genes', 'phenotypeInheritance','chromosome', 'comments',
columns_to_keep = ['genes', 'hgnc_synonyms',  'phenotypeInheritance_mapped', 'phenotype', 'geneMimNumber', 'phenotypeMimNumber']
omim_df = omim_df[columns_to_keep]
omim_df.to_csv(omim_output_file, sep='\t', index=False)


