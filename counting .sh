#!/bin/bash
#
#######################################################################################################################
###INFORMATION ABOUT THE SCRIPT###
# The script counts length distribution of sequences so we can make histogram.

#######################################################################################################################
###SPECIFY DATA VARIABLES###
# Change PROJECT_DIR variable to your favorite storage. Results from further analysis
#  steps will be stored here.
PROJECT_DIR=/storage/brno2/home/marek_bfu/Bi5444 # path to project dir
DATASET_DIR=$PROJECT_DIR/trimming # path to input filtered sequences
OUTPUT_DIR=$PROJECT_DIR/counts # path to output sequences

#######################################################################################################################
###SCRIPT BODY###
mkdir $OUTPUT_DIR
cd $DATASET_DIR

for file in *mirna.fastq.gz
do
   # leaves only sequences in fastq file, removes everything else. Gunzip unpacks the file and sends output on STDOUT; 
   # sed removes all lines except nucleotide sequence; awk prints lengths of sequencs; uniq counts the sequencesl;
   # sort is there because uniq works on adjacent matching lines (see manual)
   gunzip -c $file | sed -n 'n;p;n;n;' | awk '{print length}' | sort | uniq -c > $OUTPUT_DIR/${file:0:9}_counts.txt
done
