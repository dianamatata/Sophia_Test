# https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#download
# code to get VEP

 # Xcode and GCC
gcc -v

# Perlbrew
curl -L http://install.perlbrew.pl | bash
echo 'source $HOME/perl5/perlbrew/etc/bashrc' >> ~/.zshrc

Perlbrew
perlbrew install -j 5 --as 5.26.2 --thread --64all -Duseshrplib perl-5.26.2 --notest
perlbrew switch 5.26.2
perlbrew install-cpanm

# not working
# /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

brew install xz
brew install mysql

# get cpanm
brew install cpanminus
cpanm -l $HOME/cpanm Bio::DB::BigFile

cpanm DBI
cpanm DBD::mysql@4.050

# get BioPerl
curl -O https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.924.tar.gz
tar zxvf BioPerl-1.6.924.tar.gz
echo 'export PERL5LIB=${PERL5LIB}:##PATH_TO##/bioperl-1.6.924' >> ~/.zshrc

cpanm Test::Differences Test::Exception Test::Perl::Critic Archive::Zip PadWalker Error Devel::Cycle Role::Tiny::With Module::Build
cpanm DBD::mysql

export DYLD_LIBRARY_PATH=/usr/local/mysql/lib/:$DYLD_LIBRARY_PATH # wait

# get VEP
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git checkout release/113
perl INSTALL.pl

# manually downloading caches
# https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/
cd $HOME/.vep
# 23G
curl -O https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
tar xzf homo_sapiens_vep_113_GRCh38.tar.gz

# running VEP
cd ensembl-vep
./vep --cache --offline \
--dir_cache /Users/dianaavalos/ensembl-vep/ \
-i /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patient_variants.vcf \
-o /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patients_out.vcf

# we can add Clinvar
# Compressed VCF file
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
# Index file
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

./vep --cache --offline --custom file=/Users/dianaavalos/Desktop/Tertiary_Research_Assignment/downloaded_data/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
--dir_cache /Users/dianaavalos/ensembl-vep/ \
-i /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patients_out.vcf \
-o /Users/dianaavalos/Desktop/Tertiary_Research_Assignment/patients_out2.vcf

./vep --id "1  230710048 230710048 A/G 1" --species homo_sapiens -o /path/to/output/output.txt --cache --offline --assembly GRCh38 --custom file=/path/to/custom_files/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN



# TODO: remove in /Users/dianaavalos/ensembl-vep/ the homo sapiens file > 23Gb


