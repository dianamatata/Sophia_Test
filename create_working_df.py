# Libraries ---------------
import pandas as pd

# vcf
vcf_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/patient_variants_vcf_df.txt"
vcf_data = pd.read_csv(vcf_file, sep='\t')

# all the variants entries
formatted_entries_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/formatted_vcf_entries.txt"
# variants failed
failed_variants_log = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/out_mobi/failed_variants_log2.txt"

# mobidetails
mobi_file = '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/merged_data_mobidetails.txt'
mobi_data = pd.read_csv(mobi_file, sep='\t')

#OMIM
omim_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/omim.txt"
omim_data = pd.read_csv(omim_file, sep='\t')

# Clinvar

# ACMG

# MPA score


