path="/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/others"
$path
cmd="gzcat $path/clinvar_20230527.vcf.gz | grep -v '^##' > $path/clinvar_temp.txt"
echo $cmd
eval $cmd


# clone mobidetails code
git clone https://github.com/vidboda/MDUtils.git


# check clinvar size
wget --spider 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20230527.vcf.gz' # Taille : 79817563 (76M)
# get clinvar
curl -O 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20230527.vcf.gz'



gzcat /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/downloaded_data/clinvar_20230527.vcf.gz | less

# first unzip file, and at the end zip before saving
gzcat /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/downloaded_data/clinvar_20230527.vcf.gz | \
awk -F'\t' 'BEGIN {
    OFS="\t";
    print "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "ALLELEID", "CLNDISDB", "CLNHGVS", "CLNSIG", "CLNVC", "CLNVCSO", "GENEINFO", "MC"
}
{
    split($8, info, ";");  # Split INFO column by ";"
    for (i in info) {
        split(info[i], kv, "=");  # Split key-value pairs
        if (kv[1] == "ALLELEID") alleleid = kv[2];
        if (kv[1] == "CLNDISDB") clndisdb = kv[2];
        if (kv[1] == "CLNHGVS") clnhgvs = kv[2];
        if (kv[1] == "CLNSIG") clnsig = kv[2];
        if (kv[1] == "CLNVC") clnvc = kv[2];
        if (kv[1] == "CLNVCSO") clnvcso = kv[2];
        if (kv[1] == "GENEINFO") geneinfo = kv[2];
        if (kv[1] == "MC") mc = kv[2];
    }
    print $1, $2, $3, $4, $5, $6, $7, $8, alleleid, clndisdb, clnhgvs, clnsig, clnvc, clnvcso, geneinfo, mc
}' | gzip > /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/clinvar_expanded_variants.txt.gz



gzcat /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/clinvar_expanded_variants.txt.gz | less

gzcat /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/clinvar_expanded_variants.txt.gz  | \
awk -F'\t' 'BEGIN { OFS="\t" }
NR == 1 { print $0, "hg38_vcf_clinvar" }  # Add new column header
NR > 1 { print $0, "hg38:" $1 ":" $2 ":" $4 ":" $5 }'  > /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/summary_variants_clinvar.txt

#clinvar_gene_data_output_file = "/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/SAMD11_variants.txt"
#clinvar_gene_data = pd.read_csv(clinvar_gene_data_output_file, sep='\t')


# for each variant,take hotspot on +25 and -25:
# check hoy many types of variants in each
# best is to extract +100/-100 from chr, pos





#  other stuff -------------------

#while read -r line; do
#  # Extract the 5th field from each line using 'cut' and append it to the output file
#  echo "$line" | cut -d " " -f5 >> /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/failed_variants_log_hg38.txt
#done < /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/failed_variants_log.txt

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
