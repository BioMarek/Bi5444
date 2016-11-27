

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


##########################################################################
#                       DOWNLOADING SCRIPT                               #
##########################################################################
#!/bin/bash
#
# Don't forget dos2unix if working on windows
# Results from each step of analysis will be in separate directory.
cd /storage/brno2/home/marek_bfu # (☞ﾟヮﾟ)☞ change to your favorite storage ☜(ﾟヮﾟ☜)
mkdir -p Bi5444/raw_sequences
cd Bi5444/raw_sequences

# Because the files have random names we cannot use for loop. ¯\_(ツ)_/¯ Once each file is downloaded it is renamed to avoid confusion.
# I'm going to download them into my storage, then rename them and only after that copy them to  $scratch and work with them.
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852089/ERR852089.fastq.gz
mv ERR852089.fastq.gz control_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852099/ERR852099.fastq.gz
mv ERR852099.fastq.gz control_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852097/ERR852097.fastq.gz
mv ERR852097.fastq.gz control_3.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852092/ERR852092.fastq.gz
mv ERR852092.fastq.gz control_4.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR852/ERR852091/ERR852091.fastq.gz
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


##########################################################################
#                    FASTQC BEFORE TRIMMING SCRIPT                       #
##########################################################################

#!/bin/bash
#
# for loop copying everything to $SCRATCH and unpacking
# on metacentrum fasqc refuses to take *.gz file as input so that's the reason I'm unpacking it.
cd /storage/brno2/home/marek_bfu/Bi5444 # (☞ﾟヮﾟ)☞ change to your favorite storage ☜(ﾟヮﾟ☜)
mkdir fastqc_before_trim
cd /storage/brno2/home/marek_bfu/Bi5444/raw_sequences

for file in *
do
  echo copying "$file" # just to know where we are
  cp $file $SCRATCH/$file
  gunzip $file
done

cd $SCRATCH
module add fastQC-0.10.1
for file in *
do
  fastqc $file
done
mv *.zip /storage/brno2/home/marek_bfu/Bi5444/fastqc_before_trim/ # don't forget to copy file to your computer

rm -rf $SCRATCH/*


##########################################################################
#                            ADAPTER SEARCHING                           #
##########################################################################

#!/bin/bash
#
# Simple script for minion adapter search
# within the folder  (DATASET_DIR) it takes all files with
# specified suffix (SUFFIX) and looks for the adapters and compares
# them with adapter file (ADAPTERS) using Swan
# Minion can also take .gz input files as itself
# Minion and Swan are part of the Kraken pipeline http://www.ebi.ac.uk/research/enright/software/kraken
#

# Set variables - input folder, output folder and suffix of files to check
DATASET_DIR=/storage/brno2/home/marek_bfu/Bi5444/raw_sequences
OUTPUT_DIR=/storage/brno2/home/marek_bfu/Bi5444/minion

ADAPTERS=/storage/brno2/home/marek_bfu/adapters_merge.txt # List of possible adapters

mkdir -p $OUTPUT_DIR # Make output directory with including all directories (up and down)

cd $DATASET_DIR # Go to the input folder

# Start minion and swan
for i in *
do	
	# I used absoulte paths to folder where I compiled swan and minion.
	/storage/brno2/home/marek_bfu/reaper-15-065/src/minion search-adapter -i $i -show 5 -write-fasta $OUTPUT_DIR/${i%.*}.minion.fasta # Identify top 5 over-represented sequences
	/storage/brno2/home/marek_bfu/reaper-15-065/src/swan -r $ADAPTERS -q $OUTPUT_DIR/${i%.*}.minion.fasta > $OUTPUT_DIR/${i%.*}.minion.compare # Compare them with list of adapters
done








##########################################################################
#                       PREPROCESSING SCRIPT                             #
##########################################################################

#!/bin/bash

###################################################################################################
###INFORMATION ABOUT THE SCRIPT###
# Pre-processing script designed for microRNAs analysis using Cutadapt (http://cutadapt.readthedocs.io/en/stable/index.html) and Fastx-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/)
# It is focused for later use of Chimira (http://www.ebi.ac.uk/research/enright/software/chimira) and DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) later on
#
# The script does following stepts in order to preprocess data for mapping
# 1) Adapter trimming
# 2) Quality trimming
# 3) Size filtering
# 4) Quality filtering 

###################################################################################################
##SPECIFY DATA VARIABLES###
INPUT_SUFFIX=".fastq.gz" # Suffix of files to launch the analysis on
DATASET_DIR=/storage/brno2/home/marek_bfu/Bi5444/raw_sequences #path to input raw sequences
OUTPUT_DIR=/storage/brno2/home/marek_bfu/Bi5444/trimming #path to output seqquences

FILE_FORMAT=fastq # File format
QUALITY=33 # Phred coding of input files

### What kind of coding is used for Phred? ##
##head -n 40 $datafile | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding\!";}'


##CUTADAPT VARIABLES#
QT_THRESHOLD=5 # Threshold for quality trimming; we filter by number of mismatches so we need high quality reads

DISC_SHORT=14 # Discard too short sequences after the pre-processing
DISC_LONG=30 # Discard too long after the pre-processing



##FASTX - QUALITY FILTERING VARIABLES###
QF_THRESHOLD=10 # Threshold for quality filtering
QF_PERC=85 # Minimal percentage of bases with $QF_THRESHOLD

##Add modules

#add Cudadapt module
module add python27-modules-gcc  
module add python27-modules-intel

#add FastX toolkit
module add fastx-0.0.13 

###################################################################################################
###SCRIPT BODY###
mkdir -p $OUTPUT_DIR # Create output directory including tmp folder for additional results

cd $DATASET_DIR # Go to the input directory

for sample in *$INPUT_SUFFIX # For each file with specified suffix in the directory do the pre-processing loop
do
	# Cutadapt adapter removal, quality trimming, N bases removal  and length filtering
	cutadapt --quality-cutoff $QT_THRESHOLD,$QT_THRESHOLD --trim-n --max-n=0 \
	--minimum-length $DISC_SHORT --maximum-length $DISC_LONG -o $OUTPUT_DIR/${sample%.fastq*}.ad3trim.fastq.gz \
	# --untrimmed-output=$OUTPUT_DIR/${sample%.fastq*}.ad3untrimmed.fastq.gz \
	#we do not have adapters, we do only size and quality filtering, so we cannot use --untrimmed-output =Write all reads without adapters to FILE (in FASTA/FASTQ format) instead of writing them to the regular output file.\
	$sample #cuadapt detects format automatically

	# Fastx-toolkig quality filtering; to use gz as input/output https://www.biostars.org/p/83237/
	gunzip -c $OUTPUT_DIR/${sample%.fastq*}.ad3trim.fastq.gz | fastq_quality_filter -Q $QUALITY \
	-q $QF_THRESHOLD -p $QF_PERC -z -o $OUTPUT_DIR/${sample%.fastq*}.mirna.fastq.gz
done

###################################################################################################
###ANALYSIS DONE###








# TODO copying resuls of fastqc from $SCRATCH to home storage

# TODO other things we are supposed to do

# cleaning $SCRATCH
rm -r $SCRATCH/*
exit # maybe not necessary

