#!/bin/bash

##########################################################################
#   I have started already, feel free to improve any part of the code.   #
##########################################################################


# Request for machine. This is just temporary, maybe we will need to change settings later.
qsub -l walltime=1h -l mem=4gb -l scratch=4gb -l nodes=1:ppn=4

# Results from each step of analysis will be in separate directory.
cd /storage/brno2/home/marek_bfu # change to your favorite storage
mkdir -p Bi5444/raw_data
cd Bi5444/raw_data

# Because the files have random names we cannot use for loop. Once each file is downloaded it is renamed to avoid confusion.
# I'm going to download them into my storage, then rename them and only after that copy them to  $scratch and work with them.
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852089/ERR852089.fastq.gz
mv ERR852089.fastq.gz control_1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852099/ERR852099.fastq.gz
mv ERR852089.fastq.gz control_2
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852097/ERR852097.fastq.gz
mv ERR852097.fastq.gz control_3
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852092/ERR852092.fastq.gz
mv ERR852092.fastq.gz control_4
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852091/ERR852091.fastq.gz
mv ERR852091.fastq.gz control_5

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852095/ERR852095.fastq.gz
mv ERR852095.fastq.gz patient_1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852093/ERR852093.fastq.gz
mv ERR852093.fastq.gz patient_2pwd
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852096/ERR852096.fastq.gz
mv ERR852096.fastq.gz patient_3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852090/ERR852090.fastq.gz
mv ERR852090.fastq.gz patient_4
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852094/ERR852094.fastq.gz
mv ERR852094.fastq.gz patient_5
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852098/ERR852098.fastq.gz
mv ERR852098.fastq.gz patient_6


# for loop moving everything to $SCRATCH
cd $SCRATCH
mkdir raw_data

# running fastqc

mkdir fastqc_before_trim
cd fastqc_before_trim
module add fastqc # and number of version?

# copying resuls of fastqc from $SCRATCH to home storage

# adapter trimming

# fastqc after trimming

# copying resuls of fastqc from $SCRATCH to home storage

# other things we are supposed to do

# cleaning $SCRATCH
rm -r $SCRATCH/*
exit # maybe not necesary

