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

export DYLD_LIBRARY_PATH=/usr/local/mysql/lib/:$DYLD_LIBRARY_PATH # wait

# get VEP
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git checkout release/113
perl INSTALL.pl

# get Bio::DB::BigFile
#wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
#tar xzf v335_base.tar.gz



# errors
cpanm Bio::EnsEMBL::Utils::Logger
cpanm Bio::EnsEMBL::IO
cpanm Bio::EnsEMBL::VEP::Config
