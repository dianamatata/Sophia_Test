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

