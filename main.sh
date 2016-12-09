

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
PROJECT_DIR=/storage/brno2/home/marek_bfu/Bi5444 # (☞ﾟヮﾟ)☞ change to your favorite storage ☜(ﾟヮﾟ☜)
mkdir -p $PROJECT_DIR/raw_sequences
cd $PROJECT_DIR/raw_sequences

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
cd $PROJECT_DIR
mkdir fastqc_before_trim
cd $PROJECT_DIR/raw_sequences

for file in *
do
  echo copying "$file" # just to know where we are
  cp $file $SCRATCH/$file
  #gunzip $file
done

cd $SCRATCH
module add fastQC-0.10.1
for file in *
do
  fastqc $file
done

mv *.zip $PROJECT_DIR/fastqc_before_trim/ # don't forget to copy file to your computer


##########################################################################
#                            ADAPTER SEARCHING                           #
##########################################################################

### TODO do it in scratch and check whether it works


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
DATASET_DIR=$PROJECT/raw_sequences
OUTPUT_DIR=$PROJECT/minion

cd $PROJECT # List of possible adapters is stored in study materials 
wget https://is.muni.cz/auth/el/1431/podzim2016/Bi5444/um/65638858/adapters_merge.txt?studium=694611
ADAPTERS=$PROJECT/adapters_merge.txt

mkdir -p $OUTPUT_DIR # Make output directory with including all directories (up and down)

cd $DATASET_DIR # Go to the input folder

### TODO downloading and compilation of swan and minion automated

# Start minion and swan
for i in *
do	
	# I used absoulte paths to folder where I compiled swan and minion.
	/storage/brno2/home/marek_bfu/reaper-15-065/src/minion search-adapter -i $i -show 5 -write-fasta $OUTPUT_DIR/${i%.*}.minion.fasta # Identify top 5 over-represented sequences
	/storage/brno2/home/marek_bfu/reaper-15-065/src/swan -r $ADAPTERS -q $OUTPUT_DIR/${i%.*}.minion.fasta > $OUTPUT_DIR/${i%.*}.minion.compare # Compare them with list of adapters
done


###################################################################################################
# Check result from swan and minion by eye, if everything is OK you can procedd to the next step.
###################################################################################################


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

####################################################################################################
####################################################################################################
# Script to calculate differential gene expression using DESeq2 and edgeR package
# Designed for miRNA differential gene expression analysis based on miRBase results#add FastX toolkit
module add fastx-0.0.13 

###################################################################################################
###SCRIPT BODY###
mkdir -p $OUTPUT_DIR # Create output directory including tmp folder for additional results

cd $DATASET_DIR # Go to the input directory

for sample in *$INPUT_SUFFIX # For each file with specified suffix in the directory do the pre-processing loop
do
	# Cutadapt adapter removal, quality trimming, N bases removal  and length filtering
	# --untrimmed-output=$OUTPUT_DIR/${sample%.fastq*}.ad3untrimmed.fastq.gz \
	# we do not have adapters, we do only size and quality filtering, so we cannot use --untrimmed-output =Write all reads without adapters to FILE (in FASTA/FASTQ format) instead of writing them to the regular output file.
	cutadapt --quality-cutoff $QT_THRESHOLD,$QT_THRESHOLD --trim-n --max-n=0 \
	--minimum-length $DISC_SHORT --maximum-length $DISC_LONG -o $OUTPUT_DIR/${sample%.fastq*}.ad3trim.fastq.gz $sample #cuadapt detects format automatically

	# Fastx-toolkit quality filtering; to use gz as input/output https://www.biostars.org/p/83237/
	gunzip -c $OUTPUT_DIR/${sample%.fastq*}.ad3trim.fastq.gz | fastq_quality_filter -Q $QUALITY \
	-q $QF_THRESHOLD -p $QF_PERC -z -o $OUTPUT_DIR/${sample%.fastq*}.mirna.fastq.gz
done

###################################################################################################
###ANALYSIS DONE###


##########################################################################
#                       	DESeq2 analysis                          #
##########################################################################
## DESeq2 is an R package, so it is necessarry to work in R, either in Metacentrum (described bellow) or install R and RStudio at home computer##

ssh -X skirit.metacentrum.cz #if we want to use graphical programs, they will open in our computer
qsub -l walltime=2h -l mem=4gb -l scratch=40gb -l nodes=1:ppn=4 -I -X #ask for interactive job, ask for graphical output to run in our window


#go to $SCRATCH and perform the analysis there
#cd $SCRATCH/

#copy data to $SCRATCH
#cp my/favourite/metacentrum/storage/file $SCRATCH/

module add rstudio
rstudio

# in rstudio File -> New File -> R Script
####################################################################################################
####################################################################################################
# Script to calculate differential gene expression using DESeq2 and edgeR package
# Designed for miRNA differential gene expression analysis based on miRBase results
#
# INPUT_COUNTS - table with raw genes counts -> first column = miRBase gene id, first row = patient id
####################################################################################################
####################################################################################################
#
#rm(list=ls(all=TRUE)) ##removes all variables that were stored in rstudio from older projects

## Install packages
install.packages("pheatmap")
install.packages("rgl")
install.packages("gplots")

# install bioconductor, run with ctrl + enter, when asked "update all/some/none" pres n
source("http://bioconductor.org/biocLite.R")
biocLite()

#install DESeq2
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")


####################################################################################################
# General variables
INPUT_COUNTS<-"/storage/brno2/home/indrakovaa/Bi5444/All_plain_counts.counts" 
OUTPUT_DIR<-"/storage/brno2/home/indrakovaa/Bi5444/Mirna_results/"

####################################################################################################
# Custom variables
P_THRESHOLD<-0.05
LFC_THRESHOLD<-log(1.5, 2) #log2 fold change, important results will be in range given in brackets
TOP<-20 # How many top DE miRNAs should be plotted

####################################################################################################
dir.create(OUTPUT_DIR, recursive = T)
setwd(OUTPUT_DIR)

# Processing
mrcounts<-read.table(INPUT_COUNTS, header=TRUE, row.names=1)

# if we need to reorder the column order we use following line
# mrcounts<-mrcounts[,c("header_column_4", "header_column_2", "header_column_1", "header_column_3")] # Select and reorder columns, in R we read tables table[rows,columns], so if I want to see first column: table [,1]

conds<-factor(c(rep("control", 5), rep("patient", 6)), levels=c("control", "patient")) # Set conditions and make sure they are in correct order, factor-in R define levels

mrcounts<-mrcounts[rowSums(mrcounts)!=0,] # Remove not expressed genes in any sample



####################################################################################################
# create table for DESeq2, 
coldata<-as.data.frame(t(t(conds))) #create table, use number of rows according to conds
colnames(coldata)<-"condition" #name header of the column "condition"
rownames(coldata)<-colnames(mrcounts) #assign rows the rownames from table mrcounts
coldata<-as.data.frame(coldata) #make sure it is a data frame (table)

sink("deseq2_design_control.txt") #save output from R into a file (it saves the output that would be in terminal into a file)
  coldata
sink()
#otveru si "deseq2_design_control.txt" a zkontroluju, jak to dopadlo
####################################################################################################
####################################################################################################
# DESeq2 part
# Make the count object, normalise, dispersion and testing, read manual DESeq2
library("DESeq2")

dds<-DESeqDataSetFromMatrix(countData = mrcounts, colData = coldata, design = ~condition)

# Make sure the levels are correct
dds$condition<-factor(dds$condition, levels=levels(coldata$condition)) # Make sure conditions are levels 

# The same thing which follows be calculated by >DESeq(dds) instead of three separate commands
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds)
dds<-nbinomWaldTest(dds) 
cds<-dds

####################################################################################################
# DESeq2 results visualization
res<-results(cds, contrast=c("condition", levels(conds)[2], levels(conds)[1])) # Extract contrasts of interest
res<-res[order(res$padj),]
sig<-rownames(res[(abs(res$log2FoldChange) >= LFC_THRESHOLD) & (res$padj < P_THRESHOLD) & !is.na(res$padj),]) # Extract "significant" results

res2<-res
res2<-merge(res2, counts(dds, normalized=TRUE), by="row.names") # Add norm. counts
rownames(res2)<-res2$Row.names
res2$Row.names<-NULL

res2<-merge(res2, counts(dds, normalized=FALSE), by="row.names") # Add raw counts
rownames(res2)<-res2$Row.names
res2$Row.names<-NULL
res2$gene_id<-NULL
res2<-res2[with(res2, order(padj)), ]

colnames(res2)<-gsub(pattern = ".x", replacement = "_normCounts", fixed = T, colnames(res2))
colnames(res2)<-gsub(pattern = ".y", replacement = "_rawCounts", fixed = T, colnames(res2))

write.table(x = res2, file = "deseq2.csv", sep = "\t", col.names = NA)

####################################################################################################
# Visualization and other
pdf(file="deseq2_disperison_plot.pdf")
  plotDispEsts(cds,main="Dispersion Plot")
dev.off()

rawcounts<-counts(cds, normalized=FALSE) # Save raw counts
normcounts<-counts(cds, normalized=TRUE) # Save normalized counts
log2counts<-log2(normcounts+1) # Save log2 of normalized counts

vsd<-varianceStabilizingTransformation(cds) # Save counts tranformed with variance Stabilizing Transformation
vstcounts<-assay(vsd)

# Normalization check
cond_colours<-brewer.pal(length(unique(conds)), "Accent")[as.factor(conds)]
names(cond_colours)<-conds

pdf(file="deseq2_pre_post_norm_counts.pdf")
  par(mfrow=c(3,1))
  barplot(colSums(rawcounts), col=cond_colours, las=2,cex.names=1.3,main="Pre Normalised Counts", ylim=c(0,(max(apply(rawcounts,2,sum)))*1.1))
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend("center", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=0.6, horiz=TRUE)
  barplot(colSums(normcounts), col=cond_colours, las=2, cex.names=1.3, main="Post Normalised Counts", ylim=c(0,(max(apply(normcounts,2,sum)))*1.1))
dev.off()

# Heatmaps
library("gplots")

pdf(file="heatmap_raw.pdf")
  heatmap.2(cor(rawcounts), trace="none", col=hmcol, main="Sample to Sample Correlation (Raw Counts)", RowSideColors=cond_colours, margins=c(9,7))
dev.off()
pdf(file="heatmap_log2.pdf")
  heatmap.2(cor(log2counts), trace="none", col=hmcol, main="Sample to Sample Correlation (Log2)", RowSideColors=cond_colours, margins=c(9,7))
dev.off()
pdf(file="heatmap_vst.pdf")
  heatmap.2(cor(vstcounts), trace="none", col=hmcol, main="Sample to Sample Correlation (VST)", RowSideColors=cond_colours, margins=c(9,7))
dev.off()

vstcounts<-vstcounts[apply(vstcounts, 1, max) != 0,]
pca<-princomp(vstcounts)

pdf(file="sample_to_samples_PCA.pdf")
  par(mfrow=c(1,1), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
  text(pca$loadings, as.vector(colnames(mrcounts)), pos=3)
  legend("topright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
dev.off()

pca2<-prcomp(t(vstcounts),scale=TRUE,center=TRUE)

pdf(file="contributions_PCA.pdf")
  plot(pca2, type = "l", main="Principal Component Contributions")
dev.off()

pdf(file="samples_to_sample2_PCA.pdf")
  par(mfrow=c(1,3), oma=c(2,0,0,0), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
  plot(pca2$x, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
  text(pca2$x, as.vector(colnames(mrcounts)), pos=3, cex=1)
  plot(pca2$x[,1],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC1")
  text(pca2$x[,1],pca2$x[,3], as.vector(colnames(mrcounts)), pos=3, cex=1)
  plot(pca2$x[,2],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC2")
  text(pca2$x[,2],pca2$x[,3], as.vector(colnames(mrcounts)), pos=3, cex=1)
  legend(-75,-38, levels(conds), fill=cond_colours[levels(conds)], cex=1) 
dev.off()

library("rgl")

plot3d(pca$loadings, col=cond_colours,  pch=19, size=20, main="Sample to Sample PCA (VST)")
rgl.postscript("samples_to_sample_3D_PCA.pdf", "pdf")

# Quick check of DE genes
tmpMatrix<-matrix(ncol=1, nrow=5)
rownames(tmpMatrix)<-c("total genes", paste("LFC >= ", round(LFC_THRESHOLD, 2), " (up)", sep=""), paste("LFC <= ", -(round(LFC_THRESHOLD, 2)), " (down)", sep=""), "not de", "low counts")
tmpMatrix[1,1]<-nrow(res2)
tmpMatrix[2,1]<-nrow(res2[res2$log2FoldChange >= (LFC_THRESHOLD) & (res2$padj < P_THRESHOLD) & !is.na(res2$padj),])
tmpMatrix[3,1]<-nrow(res2[(res2$log2FoldChange <= (-LFC_THRESHOLD)) & (res2$padj < P_THRESHOLD) & !is.na(res2$padj),])
tmpMatrix[4,1]<-nrow(res2[((res2$padj >= P_THRESHOLD) | ((res2$log2FoldChange > (-LFC_THRESHOLD)) & (res2$log2FoldChange < (LFC_THRESHOLD)))) & !is.na(res2$padj),])
tmpMatrix[5,1]<-sum(is.na(res2$padj))
tmpMatrix[,1]<-paste(": ", tmpMatrix[,1], ", ",round(tmpMatrix[,1]/((tmpMatrix[1,1]/100)), 1), "%", sep="")

sink("deseq2_de_genes_check.txt")
  print(paste("Number of DE Genes With adj.pval < ", P_THRESHOLD, " Without LogFC Cut-off", sep=""))
  summary(res, alpha=P_THRESHOLD)
  print(paste("Number of DE Genes With adj.pval < ", P_THRESHOLD, " and LogFC >= ", round(LFC_THRESHOLD, 2), sep=""))
  noquote(tmpMatrix)
sink()

####################################################################################################
#
TOP_BCKP<-TOP
  
if(length(sig)<TOP){ # Avoid error by ploting more TOP genes than significantly DE genes
   TOP<-length(sig)
}
  
pdf(file=paste("volcanoplot_", levels(conds)[1], "_v_", levels(conds)[2],".pdf", sep=""))
    par(mfrow=c(1,1), xpd=NA) # http://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r
    plot(res$log2FoldChange, -log(res$padj, 10), main=paste("Volcanoplot ",levels(conds)[2], " v ", levels(conds)[1], " top ", TOP, " genes", sep=""), cex=0.4, pch=19)
    #   text(res[1:length(sig),]$log2FoldChange, -log(res[1:length(sig),]$padj,10), 
    #        labels=parsedEnsembl[parsedEnsembl$gene_id %in% rownames(res[1:length(sig),]), "gene_name"], cex=0.3, pos=3)
    text(res[1:TOP,]$log2FoldChange, -log(res[1:TOP,]$padj,10),
         labels=parsedEnsembl[rownames(res)[1:TOP], "gene_name"], cex=0.3, pos=3)
    legend("bottomleft", levels(conds)[1], cex=0.5)
    legend("bottomright", levels(conds)[2], cex=0.5)
    abline(h=-log(P_THRESHOLD, 10), col="red", lty=2)
    abline(v=c(-LFC_THRESHOLD, LFC_THRESHOLD), col="darkblue", lty=2)
dev.off()
 
TOP<-TOP_BCKP

# Heatmaps of selected genes
library("pheatmap")

res_tmp<-res[order(res$padj, -abs(res$log2FoldChange), -(res$baseMean)),]
select<-rownames(res_tmp[(abs(res_tmp$log2FoldChange) >= LFC_THRESHOLD) & (res_tmp$padj < P_THRESHOLD) & !is.na(res_tmp$padj),])

if(length(select)>TOP){
  select<-select[1:TOP]
}

nt<-normTransform(dds) # defaults to log2(x+1)
# nt<-rlog(dds, blind=FALSE)
# nt<-varianceStabilizingTransformation(dds, blind=FALSE)
log2.norm.counts<-assay(nt)[select,]
log2.norm.counts<-log2.norm.counts[order(rowMeans(log2.norm.counts)),]

df<-as.data.frame(colData(dds)[,c("condition","patient")])

TOP_BCKP<-TOP

if(TOP>length(select)){
  TOP<-length(select)
}

pdf(file="deseq2_heatmaps_selected_orderBaseMeanCluster.pdf")
  pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, 
           main = paste("Top ", TOP, " significantly DE genes (log2norm)", sep=""))
dev.off()

pdf(file="deseq2_heatmaps_selected_orderBaseMean.pdf")
  pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=FALSE, annotation_col=df, main = paste("Top ", TOP, " significantly DE genes (log2norm)", sep=""))
dev.off()

TOP<-TOP_BCKP

graphics.off()

system("for i in deseq2_heatmaps_selected_*; do pdftk $i cat 2-end output tmp.pdf; mv tmp.pdf $i; done") # Cut first empty page from heatmap plot







# cleaning $SCRATCH
rm -r $SCRATCH/*
exit # maybe not necessary

