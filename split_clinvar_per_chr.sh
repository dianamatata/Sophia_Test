#!/bin/bash

# Define the input file and the output directory
input_file="/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/clinvar_expanded_variants.txt.gz"
output_dir="/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/clinvar_chr/"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over each chromosome and process
for chrom in {1..22} X Y; do
  output_file="${output_dir}clinvar_${chrom}.txt.gz"
  # Extract lines where the first column matches the chromosome and compress
  echo "Processing chromosome $chrom..."
  gzcat "$input_file" | awk -v chrom="$chrom" '$1 == chrom' | gzip > "$output_file"
  echo "Chromosome $chrom data written to $output_file"
done

echo "Splitting by chromosome completed."

rm $input_file