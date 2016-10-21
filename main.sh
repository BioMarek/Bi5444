#!/bin/bash

##########################################################################
#   I have started already, feel free to improve any part of the code.   #
##########################################################################
#                             General idea                               #
# Download files to storage - only once                                  #
# Copy files to $SCRATCH and do all the work there                       #
# After each step in pipeline transfer meanigful results back to storage #
##########################################################################


# Request for machine. This is just temporary, maybe we will need to change settings later.
qsub -l walltime=2h -l mem=4gb -l scratch=40gb -l nodes=1:ppn=4

# Results from each step of analysis will be in separate directory.
cd /storage/brno2/home/marek_bfu # change to your favorite storage
mkdir -p Bi5444/raw_data
cd Bi5444/raw_data

# Because the files have random names we cannot use for loop. Once each file is downloaded it is renamed to avoid confusion.
# I'm going to download them into my storage, then rename them and only after that copy them to  $scratch and work with them.
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852089/ERR852089.fastq.gz
mv ERR852089.fastq.gz control_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852099/ERR852099.fastq.gz
mv ERR852089.fastq.gz control_2.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852097/ERR852097.fastq.gz
mv ERR852097.fastq.gz control_3.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852092/ERR852092.fastq.gz
mv ERR852092.fastq.gz control_4.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852091/ERR852091.fastq.gz
mv ERR852091.fastq.gz control_5.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852095/ERR852095.fastq.gz
mv ERR852095.fastq.gz patient_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852093/ERR852093.fastq.gz
mv ERR852093.fastq.gz patient_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852096/ERR852096.fastq.gz
mv ERR852096.fastq.gz patient_3.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852090/ERR852090.fastq.gz
mv ERR852090.fastq.gz patient_4.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852094/ERR852094.fastq.gz
mv ERR852094.fastq.gz patient_5.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852098/ERR852098.fastq.gz
mv ERR852098.fastq.gz patient_6.fastq.gz


# for loop copying everything to $SCRATCH and unpacking - NEEDS TEST!!!
# on metacentrum fasqc refuses to take *.gz file as input so that's the reason I'm unpacking it.
mkdir $SCRATCH/raw_data
for file in *
do
  cp $file $SCRATCH/raw_data/$file
  gunzip $file
done

# TODO running fastqc
module add fastQC-0.10.1
mkdir fastqc_before_trim
cd fastqc_before_trim


# TODO copying resuls of fastqc from $SCRATCH to home storage

# TODO adapter trimming

# TODO fastqc after trimming

# TODO copying resuls of fastqc from $SCRATCH to home storage

# TODO other things we are supposed to do

# cleaning $SCRATCH
rm -r $SCRATCH/*
exit # maybe not necessary

