#! /bin/bash
#
##############################################################################################################################
###INFORMATION ABOUT THE SCRIPT###
# Simple script for minion adapter search using Minion and Swan which are part of the Kraken pipeline 
# (http://www.ebi.ac.uk/research/enright/software/kraken) 
#
# The script performs following steps:
# 1) downloads file containing sequences of possible adapters
# 2) dowloads and installs minion and swan
# 3) takes all files and looks for the adapters and compares them with adapter file using Swan

##############################################################################################################################
###SPECIFY DATA VARIABLES###
# Change PROJECT_DIR variable to your favorite storage. Results from further analysis steps will be stored here.
PROJECT_DIR=/storage/brno2/home/marek_bfu/Bi5444
DATASET_DIR=$PROJECT_DIR/raw_sequences
OUTPUT_DIR=$PROJECT_DIR/minion

##############################################################################################################################
###SCRIPT BODY###
# Downloading of adapters, list of possible adapters is stored in study materials
cd $PROJECT_DIR
wget https://is.muni.cz/el/1431/podzim2016/Bi5444/um/65638858/adapters_merge.txt
ADAPTERS=$PROJECT_DIR/adapters_merge.txt

# Downloading and instalation of minion nad swan
wget http://wwwdev.ebi.ac.uk/enright-dev/kraken/reaper/src/reaper-15-065.tgz 
tar zxvf reaper-15-065.tgz
cd reaper-15-065/src
make # command that compiles minion and swan

# Make output directory with including all directories (up and down)
mkdir -p $OUTPUT_DIR 
cd $DATASET_DIR

# Start minion and swan. I used absoulte paths to folder where I compiled swan and minion.
for i in *
do	
	# Identify top 5 over-represented sequences
	$PROJECT_DIR/reaper-15-065/src/minion search-adapter -i $i -show 5 -write-fasta $OUTPUT_DIR/${i%.*}.minion.fasta
	# Compare them with list of adapters
	$PROJECT_DIR/reaper-15-065/src/swan -r $ADAPTERS -q $OUTPUT_DIR/${i%.*}.minion.fasta > $OUTPUT_DIR/${i%.*}.minion.compare 
done
