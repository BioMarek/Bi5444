#!/bin/bash
#
# Following script performs quality check using FastQC.
# for loop copying everything to $SCRATCH
PROJECT_DIR=/storage/brno2/home/marek_bfu/Bi5444
cd $PROJECT_DIR
mkdir fastqc_before_trim
cd $PROJECT_DIR/raw_sequences

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

mv *.zip $PROJECT_DIR/fastqc_before_trim/ # copies results to our storage directory
rm -rf *
