path="/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/others"
$path
cmd="gzcat $path/clinvar_20230527.vcf.gz | grep -v '^##' > $path/clinvar_temp.txt"
echo $cmd
eval $cmd


while read -r line; do
  # Extract the 5th field from each line using 'cut' and append it to the output file
  echo "$line" | cut -d " " -f5 >> /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/failed_variants_log_hg38.txt
done < /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/failed_variants_log.txt


# TODO understand why it does not work with "$failed_variants_log"
## extract variants that did not work
#failed_variants_log="/Users/dianaavalos/Desktop/failed_variants_log.txt"
#failed_variants_log2="/Users/dianaavalos/Desktop/failed_variants_log2.txt"
#
## Clear the output file before writing to it
#> "$failed_variants_log2"
#
## Process each line and extract the 5th field
#while read -r line; do
#  # Extract the 5th field from each line using 'cut' and append it to the output file
#  echo "$line" | cut -d " " -f5 >> "$failed_variants_log2"
#done < "$failed_variants_log"
