#!/bin/bash
#
##############################################################################################################################
###INFORMATION ABOUT THE SCRIPT###
# Script designed to perform quality check using FastQC version 0.10.1. (http://www.bioinformatics.babraham.ac.uk/projects/fastqc)
#
# 1) copping sequences to $SCRATCH
# 2) quality analysis
# 3) coping compresed FastQC from $SCRATCH to our directory

##############################################################################################################################
###SPECIFY DATA VARIABLES###
# Change PROJECT_DIR variable to your favorite storage. Results from further analysis steps will be stored here.
PROJECT_DIR=/storage/brno2/home/marek_bfu/Bi5444

##############################################################################################################################
###SCRIPT BODY###
cd $PROJECT_DIR
mkdir fastqc_before_trim #creates folder for results of this analysis step
cd $PROJECT_DIR/raw_sequences

# for loop copying everything to $SCRATCH
for file in *
do
  echo copying "$file" # tells us which file is being copied, just to know where we are
  cp $file $SCRATCH/$file
done

# FastQC analysis
cd $SCRATCH
module add fastQC-0.10.1
for file in *
do
  fastqc $file
done

# copies results to our storage directory
mv *.zip $PROJECT_DIR/fastqc_before_trim/ 

rm -rf  $SCRATCH/*
