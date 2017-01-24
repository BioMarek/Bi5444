### TODO do it in scratch and check whether it works

# Simple script for minion adapter search
# within the folder  (DATASET_DIR) it takes all files with
# specified suffix (SUFFIX) and looks for the adapters and compares
# them with adapter file (ADAPTERS) using Swan
# Minion can also take .gz input files as itself
# Minion and Swan are part of the Kraken pipeline http://www.ebi.ac.uk/research/enright/software/kraken
#

# Set variables - input folder, output folder and suffix of files to check
DATASET_DIR=$PROJECT_DIR/raw_sequences
OUTPUT_DIR=$PROJECT_DIR/minion

# Downloading of adapters, list of possible adapters is stored in study materials
cd $SCRATCH
wget https://is.muni.cz/el/1431/podzim2016/Bi5444/um/65638858/adapters_merge.txt
ADAPTERS=$SCRATCH/adapters_merge.txt

# Downloading and instalation of minion nad swan
wget http://wwwdev.ebi.ac.uk/enright-dev/kraken/reaper/src/reaper-15-065.tgz 
tar zxvf reaper-15-065.tgz
cd reaper-15-065/src
make

# Make output directory with including all directories (up and down)
mkdir -p $OUTPUT_DIR 
cd $DATASET_DIR

# Start minion and swan. I used absoulte paths to folder where I compiled swan and minion.
for i in *
do	
	# Identify top 5 over-represented sequences
	$SCRATCH/reaper-15-065/src/minion search-adapter -i $i -show 5 -write-fasta $OUTPUT_DIR/${i%.*}.minion.fasta
	# Compare them with list of adapters
	$SCRATCH/reaper-15-065/src/swan -r $ADAPTERS -q $OUTPUT_DIR/${i%.*}.minion.fasta > $OUTPUT_DIR/${i%.*}.minion.compare 
done
