# Sophia_Test

SOPHiA GENETICS job application– 4-Day Assignment


Assignment Overview:
The proposed task consists of:

Annotating the variants in the attached file.
Classifying them according to their pathogenicity.
Providing a report on the findings.
The goal of this task is to assess the candidate’s knowledge in genomics and computational skills, including:

Evaluating variant pathogenicity.
Clinically interpreting the variant.

## Data Processing Workflow

### **1. [format_vcf.py](format_vcf.py)** 
**Function**: Processes the `patient_variants.vcf` file and formats it into a DataFrame for further analysis.  
**Input**:
- `patient_variants.vcf`: The raw VCF file containing patient variants.  
**Output**:
- `data/patient_variants_vcf_df.txt`: The processed VCF data stored in a tabular format.  
- `data/formatted_vcf_entries.txt`: Contains formatted entries of VCF variants, which will be used in the next script.  

---

### **2. [extract_mobidetails_info.py](extract_mobidetails_info.py) **
**Function**: Extracts information related to mobidetails from the formatted VCF entries.  
**Input**:
- `data/formatted_vcf_entries.txt`: Formatted VCF entries from `format_vcf.py`.  
**Output**:
- Extracted variants on mobidetails can be found there: `out_mobi/output_mobidetails`

---

### **3. [merge_output_mobi.py](merge_output_mobi.py) **
**Function**: Merges the outputs from `extract_mobidetails_info.py` located in `out_mobi/output_mobidetails` and creates a combined dataset.  
**Input**:
- Extracted data from `extract_mobidetails_info.py` in `out_mobi/output_mobidetails`.  
**Output**:
- `'merged_data_mobidetails.txt'`: Merged data containing mobidetails information.  

TO DO need to extract and understand failed variants

TO DO check that all the variants are present - which ones are missing

---

### **4. [create_working_df.py](create_working_df.py) **
**Function**: Merges the VCF data with the mobidetails data to create a working DataFrame.  
**Input**:
- `data/patient_variants_vcf_df.txt`: Processed VCF data from `format_vcf.py`.  
- `merged_data_mobidetails.txt`: Merged data with mobidetails from `merge_output_mobi.py`.  
**Output**:
- `merged_mobidetails_vcf.txt`: The final merged file with both VCF and mobidetails data.  

---

### **5. [sync_OMIM.py](sync_OMIM.py) **
**Function**: Processes OMIM data, formats it, and maps inheritance for all genes.  
**Input**:
- `downloaded_data/full_omim_table.txt`: Raw OMIM data file.  
**Output**:
- `omim.txt`: Formatted OMIM data with inheritance information for each gene.  

---

### **6. [merge_mobi_OMIM.py](merge_mobi_OMIM.py) **
**Function**: Merges the OMIM data with the mobidetails data.  
**Input**:
- `merged_data_mobidetails.txt`: Merged VCF and mobidetails data from `create_working_df.py`.  
- `omim.txt`: OMIM data with mapped inheritance from `sync_OMIM.py`.  
**Output**:
- `mobi_data_with_omim_genes.txt`: Final merged file containing both mobidetails and OMIM data.  

TO DO rerun and check if some OMIM data is missing

---

### **7. [get_clinvar_data.sh](get_clinvar_data.sh) **
**Function**: Get Clinvar data, select a subset of column, format and save
**Input**:
**Output**:
- `data/summary_variants_clinvar.txt`: Formatted Clinvar data

---


### **8. [split_clinvar_per_chr.sh](split_clinvar_per_chr.sh)**
**Function**:  
- extract window of variants around a given: 
**Input**: 
- clinvar_directory, chromosome, position
**Output**:
- clinvar_subset_data

**Function**:
- extract_clinvar_gene_of_variant
**Input**: 
- clinvar_directory, chromosome, gene name
- **Output**:
- clinvar_subset_data

- 
---

### **9. [extract_clinvar_data_for_variant.py](extract_clinvar_data_for_variant.py)**
**Function**:  extract_clinvar_window_around_variant to check if it is a hotspot, extract_clinvar_gene_of_variant
For ACMG crieria, we need to know for instance if "missenses" are a common mechanism of disease in the gene,
or if the variant id located in a hotspot
# TODO implement systematic check in function 
**Input**: 
**Output**:

---

### **10. [extract_clinvar_entries.py](extract_clinvar_entries.py)`extract_clinvar_entries.py`**
**Function**: Sync with clinvar website. For a given dbSNP rsid, Take the dataframe with clinvar entries P/LP, and check how many entries are found, to adapt the ACMG
criteria afterwards. 
**Input**: 
- mobi_data dataframe
**Output**:
- for instance we get 15 entries clinvar(8P/ 7LP)

---
### **11. [extract_hpo.py](extract_hpo.py)**
**Function**: for each hpo term in '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patient_phenotype.tsv', 
it browses the website 'https://hpo.jax.org/browse/term/{hpo_id}' to extract the term name and all the gene-associations
**Input**: 
**Output**:
- '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/hpo_summary.txt' with ['HPO_ID', 'Term_Name']
- '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/hpo_terms_and_gene_associations.csv' with ['HPO_ID', 'Term_Name', 'Gene_Associations':]

The code is not inside a function and needs to be run as it is.
---

### **12. [check_gene_in_hpo.py](check_gene_in_hpo.py) **
**Function**: get_hpo_terms_from_gene look up if gene name is associated to any hpo terms of the patient
**Input**: gene, hpo_data= '/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/data/hpo_terms_and_gene_associations.csv'
**Output**: associated_hpo_dict in which the hpo terms are provided in the form of a dictionnary: {'HP:0000486': 'Strabismus ', 'HP:0002119': 'Ventriculomegaly '}

TO DO implement this function to the mobi_data dataframe, right now it is run on the side for specific genes of interest.

### **13. [create_ACMG_ranking.py](create_ACMG_ranking.py)**
**Function**:  Create ACMMG ranking of all variants
Script not finished. Implemented :PVS1, PM2, PS1

---

## Explanation of the implementation

### Goal 1:
The objective is to create a file containing all annotated variants.

---

### 1. Create VCF File Formatted for Mobidetails Input
First, I create a file with VCF variants formatted in the way Mobidetails accepts as input, specifically formatted like `hg38:12:8932856:C:A`. This step ensures the data is in the correct format for Mobidetails processing.
with `format_vcf.py`
I separated in different columns:
- DP (Depth of Coverage): The total number of reads covering the position for the given variant
- AD (Allelic Depth): The number of reads supporting each allele (the reference and alternate alleles)
- DP4 (Depth of Coverage for Strand-Specific Data): The number of reads supporting the reference and alternate alleles on each strand (forward and reverse).
- And computed the Major AF and minor AF

---

### 2. Advantages of Mobidetails
Mobidetails allows us to extract a variety of data at once, including:

- **HGNC gene symbol (ID)**
- **HGVS strict genomic (hg38)**
- **HGVS Protein**: Impact on protein, such as `p.(Arg399Leu)`, `p(?)`, `p(Glu187del)`
- **HGVS genomic (hg38)**
- **Pseudo VCF (hg38)** (TODO: Need to reformat it into entries)
- **Position in transcript**: For instance, `Exon 3`, `Intron 6`
- **Position / splice site**: For example, `270 bp from acceptor`, `52 bp from donor`
- **Position / domain**: For example, `ABC transmembrane type-1 (247 - 530)`, `Region: Disordered (1275 - 1501)`, `Not in a UNIPROT-defined domain`
- **Wild type sequence** (Mutant sequence is not retained)
- **Position tolerance** according to Metadome (https://stuart.radboudumc.nl/metadome/dashboard)
- **gnomAD data**: Includes exome, genome, and non-cancer gnomAD (v2.1.1) and gnomAD v4 (v4.1.0), including global Minor Allele Frequency (MAF)
- **dbSNP rsid**: For example, `rs139044494` (Note: we removed it)
- **Clinvar Germline**: No match, Pathogenic, Benign, etc. (TODO: Add ACMG functional studies for Drug Response)
- **hg38 InterVar, GeneBe**: Some ACMG predictors, though we are building our own
- **CADD phred score**: (https://cadd.gs.washington.edu/)
- **MPA score and MPA impact**: MoBiDiC Prioritization Algorithm (https://github.com/mobidic/MPA), using input from:
    - Curated database: ClinVar
    - Biological assumption: refGene
    - Splicing prediction: SpliceAI, dbscSNV
    - Missense prediction: dbNSFP
- **Position / protein**: For example, `118 / 123`. This may need to be expressed as a percentage. For frame-shift or nonsense variants, estimate the percentage of protein lost. Also, assess if the variant causes non-mediated decay.
- **Splicing predictors**: dbscSNV ADA, dbscSNV RF, spliceAI AG, spliceAI AL, spliceAI DG, spliceAI DL, AbSplice, SPiP Interpretation, SPiP Risk
- **Missense predictors**: SIFT, Polyphen 2 HumDiv, Polyphen 2 HumVar, Fathmm, AlphaMissense, REVEL, ClinPred, Meta SVM, Meta LR, Mistic
- **LOVD database**: LOVD Matches, LOVD Effect Reported


- **TODO: Get VEP working for missing variants. or MPA? **

---

### 3. Challenges with the API
I attempted to run VEP from the command line but had some issues building it, check : [get_vep_working.sh](get_vep_working.sh)`get_vep_working.sh`
TO DO: solve this issue

I attempted to run the Mobidetails API, but the documentation was minimal and I ended up spending a lot of time trying to implement it. To address this, I designed a web scraping Python script (`extract_mobidetails_info.py`).

- **Advantages**: This approach is faster to implement.
- **Disadvantages**: Sometimes, Mobidetails takes a long time to load, possibly due to website saturation. As a result, the script may break. The script can be run and monitored in PyCharm for online debugging, but it isn't the cleanest way to handle this.
When checking manually, I saw that actually **some: quantify** entries were missing from mobidetails, like failed to load variant, or transcript. Which is annoying cause it might be rare variants..
This means I need to find other ways to retrieve the data, with MPA and Gnomad

---

### 4. Handling Missing Data

After running the script, I noticed that some variants failed to be extracted. A log file containing these failed variants was generated, and I reran the process. Some missing variants were successfully retrieved, but others were not. Upon manual inspection, I found that some entries were missing from Mobidetails, either because the variant or transcript failed to load. This is frustrating, as these could potentially be rare variants.
TO DO: understand why fo we have so many missing variants, we were expecting 27 with missing transcripts in mobidetails but actually there are 324..

---

By following this workflow, we aim to extract and annotate as much variant data as possible, while addressing the limitations of the current tools and processes. The TODOs will be addressed in future work to improve the robustness of the system.

### Goal 2:
ACMG Criteria: implemented in [create_ACMG_ranking.py](create_ACMG_ranking.py)



### Goal 3:
Processing the data.
If the variant is heterozygous with the reference allele, we discard this line.
If the variant is heterozygous and OMIM inheritance is AD, or variant is homozygous and  OMIM inheritance is AR: check variant
If splicing prediction is 1 (we are computing the splice prediction from dbscSNV ADA, dbscSNV RF, spliceAI AG, spliceAI AL, spliceAI DG, spliceAI DL, AbSplice, SPiP Interpretation and SPiP Risk ): check variant

---
## Improving perspectives:

To enhance the project, we could consider the following improvements:

- **Fix Mobidetails API**: Resolve issues with the Mobidetails API to ensure reliable data extraction and improve overall functionality.
- **Manage to integrate VEP**: Resolve the issues to run VEP on our vcf. Indeed, there are 324 variants that Mobidetails did not manage to get (some have unavailable transcript but I believe part of it is my code?)

`get_vep_working.sh`: error faced:

#Test Summary Report
#-------------------
#./t/Runner.t                                       (Wstat: 65280 Tests: 80 Failed: 0)
#Non-zero exit status: 255
#Parse errors: No plan found in TAP output
#Files=42, Tests=1808, 137 wallclock secs ( 0.29 usr  0.16 sys + 124.87 cusr  9.53 csys = 134.85 CPU)
#Result: FAIL
#Failed 1/42 test programs. 0/1808 subtests failed.
- **Integrate GERP RS**: Incorporate GERP (Genomic Evolutionary Rate Profiling) scores to provide insights into the conservation of variants across species.
- **Sync with gnomAD**: Connect with the gnomAD database to obtain ancestry-specific frequencies and counts for heterozygous and homozygous variants.
- **Integrate GTEx expression data**: GTEx data enables us to understand in which tissue the gene is expressed. If there is no OMIM phenotype associated, it can help understand the gene product
- **Integrate GenCC**: sometimes, OMIM data is not complete, and GenCC is a good option to consider
- **Integrate MoBiDiC Prioritization Algorithm**: Leverage the MoBiDiC Prioritization Algorithm (https://github.com/mobidic/MPA) for more robust variant prioritization.
- **Integrate pLoF and pLI**: from gnomAD in gene page: retrieve pLoF (Probability of Loss of Function), pLI (Probability of being Loss-of-function Intolerant) and o/e (Observed/Expected ratio).
if pLI close to 0 suggest the gene tolerates LoF mutations well, if pLI close to 1,  the gene is highly intolerant to LOF.
- **Correct OMIM integration**: some genes are associated with many OMIM phenotype entries, and I have the impression that some have gone missing
- 


 
