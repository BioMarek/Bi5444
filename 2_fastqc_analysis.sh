#!/bin/bash
#
##############################################################################################################################
###INFORMATION ABOUT THE SCRIPT###
# Script designed to perform quality check using FastQC version 0.10.1. (http://www.bioinformatics.babraham.ac.uk/projects/fastqc)
#
# The script performs following steps:
# 1) copping sequences to $SCRATCH
# 2) quality analysis
# 3) copping compresed FastQC results from $SCRATCH to our directory

##############################################################################################################################
###SPECIFY DATA VARIABLES###
# Change PROJECT_DIR variable to your favorite storage. Results from further analysis steps will be stored here.
PROJECT_DIR=/storage/brno2/home/marek_bfu/Bi5444
OUTPUT_DIR=$PROJECT_DIR/fastqc_before_trim

###ADD MODULES###
# add FastQC module
module add fastQC-0.10.1

##############################################################################################################################
###SCRIPT BODY###
# Make output directory with including all directories (up and down)
mkdir -p $OUTPUT_DIR
cd $PROJECT_DIR/raw_sequences

# for loop copies everything to $SCRATCH
for file in *
do
  echo copying "$file" # tells us which file is being copied, just to know where we are
  cp $file $SCRATCH/$file
done

cd $SCRATCH

# FastQC analysis
for file in *
do
  fastqc $file
done

# copies results to our storage directory
mv *.zip $OUTPUT_DIR

rm -rf  $SCRATCH/*
