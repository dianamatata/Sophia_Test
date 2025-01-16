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

# Clinvar

# ACMG

# MPA score

# Concatenate OMIM columns ----------------













