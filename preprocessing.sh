#!/bin/bash
#
##############################################################################################################################
###INFORMATION ABOUT THE SCRIPT###
# Pre-processing script designed for microRNAs analysis using Cutadapt (http://cutadapt.readthedocs.io/en/stable/index.html) 
# and Fastx-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/) It is focused for later use of Chimira 
# (http://www.ebi.ac.uk/research/enright/software/chimira) and DESeq2 
# (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) later on
#
# The script does following steps in order to preprocess data for mapping:
# 1) Adapter trimming
# 2) Quality trimming
# 3) Size filtering
# 4) Quality filtering 

##############################################################################################################################
##SPECIFY DATA VARIABLES###
PROJECT_DIR=/storage/brno2/home/marek_bfu/Bi5444
INPUT_SUFFIX=".fastq.gz" # Suffix of files to launch the analysis on
DATASET_DIR=$PROJECT_DIR/raw_sequences # path to input raw sequences
OUTPUT_DIR=$PROJECT_DIR/trimming # path to output sequences

## Following line of code detects quality encoding (Phred+33, Phred+64, Solexa+64), uncomment if you don't know your encoding
# head -n 40 $datafile | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding\!";}'

FILE_FORMAT=fastq # File format
QUALITY=33 # Phred coding of input files, change if you have something else

###CUTADAPT VARIABLES###
QT_THRESHOLD=5 # Threshold for quality trimming; we filter by number of mismatches so we need high quality reads

# selecting sequences that have lenght of miRNAs 
DISC_SHORT=15 # Discard too short sequences after the pre-processing
DISC_LONG=26 # Discard too long after the pre-processing

###FASTX - QUALITY FILTERING VARIABLES###
QF_THRESHOLD=10 # Threshold for quality filtering
QF_PERC=85 # Minimal percentage of bases with $QF_THRESHOLD

###ADD MODULES###
# add Cudadapt module
module add python27-modules-gcc  
module add python27-modules-intel

# add FastX toolkit
module add fastx-0.0.13 

##############################################################################################################################
###SCRIPT BODY###
# Make output directory with including all directories (up and down)
mkdir -p $OUTPUT_DIR

cd $DATASET_DIR

for sample in *$INPUT_SUFFIX # For each file with specified suffix in the directory do the pre-processing loop
do
	# Cutadapt adapter removal, quality trimming, N bases removal and length filtering
	# for quality filtering we need Phred coding +33, otherwise --quality-base=$QUALITY
	# --untrimmed-output=$OUTPUT_DIR/${sample%.fastq*}.ad3untrimmed.fastq.gz \
	# we do not have adapters, we do only size and quality filtering, so we cannot use --untrimmed-output =Write all reads without adapters to FILE (in FASTA/FASTQ format) instead of writing them to the regular output file.
	cutadapt --quality-cutoff $QT_THRESHOLD,$QT_THRESHOLD --trim-n --max-n=0 \
	--minimum-length $DISC_SHORT --maximum-length $DISC_LONG -o $OUTPUT_DIR/${sample%.fastq*}.ad3trim.fastq.gz $sample #cutadapt detects format automatically

	# Fastx-toolkit quality filtering; to use gz as input/output https://www.biostars.org/p/83237/
	# Cutadapt can trim only ends of the reads. To filter sequences with low qualitys in the middle, we need to use FastX
	gunzip -c $OUTPUT_DIR/${sample%.fastq*}.ad3trim.fastq.gz | fastq_quality_filter -Q $QUALITY \
	-q $QF_THRESHOLD -p $QF_PERC -z -o $OUTPUT_DIR/${sample%.fastq*}.mirna.fastq.gz
done
