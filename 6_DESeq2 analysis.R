
############################################################################################################################
### INFORMATION ABOUT THE SCRIPT ###
# R-Script to calculate differential gene expression using DESeq2 package
# Designed for miRNA differential gene expression analysis based on miRBase results
#
# The script performs following:
# 1) Processing of the input counts
# 2) DESeq2 analysis
# 3) DESeq2 results visualization
#
# INPUT_COUNTS - table with raw genes counts -> first column = miRBase gene id, first row = patient id
#############################################################################################################################

# removes all variables that were stored in rstudio from older projects
# rm(list=ls(all=TRUE))

### Install packages ###
install.packages("rgl")
install.packages("gplots")
install.packages("RColorBrewer")

# install bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite()

#install DESeq2
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")

###################################################################################################
###SPECIFY DATA VARIABLES###
# General variables. Change pathway to INPUT_COUNTS to where you store your alignment counts.Custom your pahtway to directory for result storage.
INPUT_COUNTS<-"d:/Dokumenty/Adela/Bi5444/Chimira_results/All_plain_counts.counts" 
OUTPUT_DIR<-"d:/Dokumenty/Adela/Bi5444/DESeq2_results/"

# Custom variables
P_THRESHOLD<-0.05 # p-value threshold
LFC_THRESHOLD<-log(1.5, 2) # Log2 fold change, important results will be in range given in brackets
TOP<-20 # How many top DE miRNAs should be plotted

####################################################################################################
###SCRIPT BODY###
# Set working directory, create output directory.
dir.create(OUTPUT_DIR, recursive = T)
setwd(OUTPUT_DIR)

# 1) Processing

mrcounts<-read.table(INPUT_COUNTS, header=TRUE, row.names=1)

# If we need to reorder the column order we use following line
# mrcounts<-mrcounts[,c("header_column_4", "header_column_2", "header_column_1", "header_column_3")]

# Set conditions and make sure they are in correct order
conds<-factor(c(rep("control", 5), rep("patient", 6)), levels=c("control", "patient")) 

# Remove not expressed genes in any sample
mrcounts<-mrcounts[rowSums(mrcounts)!=0,] 

# Create table for DESeq2
coldata<-as.data.frame(t(t(conds))) # create table, use number of rows according to conds
colnames(coldata)<-"condition" # name header of the column "condition"
rownames(coldata)<-colnames(mrcounts) #a ssign rows the rownames from table mrcounts
coldata<-as.data.frame(coldata) # make sure it is a data frame (table)

# Save output from R into a file 
sink("deseq2_design_control.txt")
coldata
sink()

####################################################################################################
####################################################################################################
# 2) DESeq2 analysis
# Make the count object, normalise, dispersion and testing, read manual DESeq2

library("DESeq2")

dds<-DESeqDataSetFromMatrix(countData = mrcounts, colData = coldata, design = ~condition)

# Make sure the levels are correct
dds$condition<-factor(dds$condition, levels=levels(coldata$condition)) # Make sure conditions are levels 

# The same thing which follows be calculated by command >DESeq(dds) instead of three separate commands
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds)
dds<-nbinomWaldTest(dds) 
cds<-dds
####################################################################################################
####################################################################################################
# 3) DESeq2 results visualization

### Create table (.csv) with normalized and raw counts ###
# Extract contrasts of interest
res<-results(cds, contrast=c("condition", levels(conds)[2], levels(conds)[1])) 
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

write.table(x = res2, file = "deseq2.csv", sep = "\t", col.names = NA) # Table with normalized and raw counts

rawcounts<-counts(cds, normalized=FALSE) # Save raw counts
normcounts<-counts(cds, normalized=TRUE) # Save normalized counts
log2counts<-log2(normcounts+1) # Save log2 of normalized counts

vsd<-varianceStabilizingTransformation(cds) # Save counts tranformed with variance Stabilizing Transformation (DESeq2 function)
vstcounts<-assay(vsd) # Stores data for PCA, data have to be rlog or varianceStabilizingTransformation

## Asign colors according conditions for data visualization
library("RColorBrewer") 
cond_colours<-brewer.pal(length(unique(conds)), "Accent")[as.factor(conds)]
names(cond_colours)<-conds

### Normalization check ###
png(file="deseq2_pre_post_norm_counts.png") # Create a file, it can be png, pdf..

# set layout of the graphs
layout(matrix(c(1,2,3),3 ,1, byrow = TRUE), # data, nrow, ncol, byrow=fill matrix by rows
       heights=c(1,0.69,1)) # respective heights of each row
# plot number of raw and normalized reads
barplot(colSums(rawcounts), col=cond_colours, las=2,cex.names=1.1,cex.axis=0.8, main="Pre Normalised Counts", ylim=c(0,(max(apply(rawcounts,2,sum)))*1.1))
plot(1, type="n", axes=F, xlab="", ylab="")
legend("center", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1, horiz=TRUE)
barplot(colSums(normcounts), col=cond_colours, las=2, cex.names=1.1, cex.axis=0.8, main="Post Normalised Counts", ylim=c(0,(max(apply(normcounts,2,sum)))*1.1))
dev.off() # Close the created files with plots

## Dispersion plot ###

png(file="deseq2_dispersion_plot.png") # Create a file, it can be png, pdf..
plotDispEsts(cds,main="Dispersion Plot")
dev.off()


## Heatmaps ###

library("gplots")

# Add colours for hmcol, colours for visualization of heatmaps, different then default
hmcol<-brewer.pal(11,"RdBu") 

# Plot heatmap with raw counts
pdf(file="heatmap_raw.pdf") # Create a file, it can be png, pdf..
heatmap.2(cor(rawcounts), trace="none", col=hmcol,
          main="Sample to Sample Correlation (Raw Counts)", RowSideColors=cond_colours, margins=c(9,7))
dev.off()

# Plot heatmap with log2 counts
png(file="heatmap_log2.png")# Create a file, it can be png, pdf..
heatmap.2(cor(log2counts), trace="none", col=hmcol, main="Sample to Sample Correlation (Log2)", RowSideColors=cond_colours, margins=c(9,7))
dev.off()

# Plot heatmap with vst counts
pdf(file="heatmap_vst.pdf") # Create a file, it can be png, pdf..
heatmap.2(cor(vstcounts), trace="none", col=hmcol, main="Sample to Sample Correlation (VST)", RowSideColors=cond_colours, margins=c(9,7))
dev.off()

### Principal component analysis ###

vstcounts<-vstcounts[apply(vstcounts, 1, max) != 0,] # Gene counts variance stabilizing transformation
pca<-princomp(vstcounts) # Principal component analysis

# plot PCA, 
png(file="sample_to_samples_PCA.png") # Create a file, it can be png, pdf..
par(mfrow=c(1,1), xpd=NA) # center legend below two plots
plot(pca$loadings, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)") # plot PCA
text(pca$loadings, as.vector(colnames(mrcounts)), pos=1, cex=0.6)
legend("bottomright", levels(unique(conds)), fill=cond_colours[levels(unique(conds))], cex=1)
dev.off()

# Dispaly relative contribution of PCA components
load<-with(pca, unclass(loadings))
aload<-abs(load) # save absolute values (=without-)

pca2<-prcomp(t(vstcounts),scale=TRUE,center=TRUE)

png(file="contributions_PCA.png") # Create a file, it can be png, pdf..
plot(pca2, type = "l", main="Principal Component Contributions") # check for batch effect
dev.off()

# Plot three graphs with principal component 1 vs 2, 1 vs 3, and 2 vs 3
png(file="samples_to_sample2_PCA.png") # Create a file, it can be png, pdf..
par(mfrow=c(1,3), oma=c(2,0,0,0), xpd=NA) # center legend below two plots
plot(pca2$x, col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)")
plot(pca2$x[,1],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC1")
plot(pca2$x[,2],pca2$x[,3], col=cond_colours,  pch=19, cex=2, main="Sample to Sample PCA (VST)",ylab="PC3",xlab="PC2")
legend(-65,-45, levels(conds), fill=cond_colours[levels(conds)], cex=1) 
dev.off()

## PLOT 3D PCA
library("rgl")
# open 3d window
open3d()
par3d(windowRect = c(100, 100, 612, 612)) # resize window
plot3d(pca$loadings, col=cond_colours, pch=100, size=2, type="s") # plot points
title3d(main = "Sample to Sample PCA (VST)", line=4, cex=1.5) # add title
legend3d("topleft", legend = paste(levels(conds)), pch = 16, col = cond_colours[levels(conds)], cex=1, inset=c(0.02)) # add legend
rgl.postscript("samples_to_sample_3D_PCA.pdf", "pdf") # write to pdf
snapshot3d(file="samples_to_sample_3D_PCA.png", top = TRUE ) # write to png

####################################################################################################
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
TOP_BCKP<-TOP

if(length(sig)<TOP){ # Avoid error by ploting more TOP genes than significantly DE genes
  TOP<-length(sig)
}

## Volcano plot
png(file=paste("volcanoplot_", levels(conds)[1], "_v_", levels(conds)[2],".png", sep="")) # Create a file, it can be png, pdf..
par(mfrow=c(1,1), xpd=NA) # center legend below two plots

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, cex=0.4, main=paste("Volcanoplot ",levels(conds)[2], " vs ", levels(conds)[1], sep="")))
# Add colored points: red if padj<0.01, orange of log2FC>2.5, green if both)
with(subset(res, res$padj<0.01 ), points(log2FoldChange, -log10(padj), pch=20, cex=0.4, col="red"))
with(subset(res, abs(res$log2FoldChange)>2.5), points(log2FoldChange, -log10(padj), pch=20, cex=0.4, col="orange"))
with(subset(res, res$padj<0.01 & abs(res$log2FoldChange)>2.5), points(log2FoldChange, -log10(padj), pch=20, cex=0.4, col="darkgreen"))
#Add labels points with the textxy function from the calibrate plot library
library(calibrate)
with(subset(res[1:9,]), textxy(log2FoldChange, -log10(padj), labs=rownames(res[1:9,]), cex=0.8, pos=1))
# Add legend
legend("bottomleft", levels(conds)[1], cex=0.5)
legend("bottomright", levels(conds)[2], cex=0.5)
# Add lines denoting area of interest
abline(h=-log(P_THRESHOLD, 10), col="red", lty=2)
abline(v=c(-LFC_THRESHOLD, LFC_THRESHOLD), col="darkblue", lty=2)
dev.off()

TOP<-TOP_BCKP

## Heatmaps of selected DE genes
library("DESeq2")
library("gplots")

res_tmp<-res[order(res$padj, -abs(res$log2FoldChange), -(res$baseMean)),]
select<-rownames(res_tmp[(abs(res_tmp$log2FoldChange) >= LFC_THRESHOLD) & (res_tmp$padj < P_THRESHOLD) & !is.na(res_tmp$padj),])

if(length(select)>TOP){
  select<-select[1:TOP]
}

nt<-normTransform(dds) # defaults to log2(x+1)
log2.norm.counts<-assay(nt)[select,]
log2.norm.counts<-log2.norm.counts[order(rowMeans(log2.norm.counts)),]

df<-as.data.frame(colData(dds) [,c("condition")]) # it only have one column named condition, so we cannot select two columns, one named condition and the other named patient

TOP_BCKP<-TOP

if(TOP>length(select)){
  TOP<-length(select)
}

# specification of the position of the color key (below) and the sizes of the graphs
lmat = rbind(c(0,3,3),c(2,1,1),c(4,0,0))
lwid = c(1.5,4,2)
lhei = c(1.5,4,1.5)

# Heatmap of DE genes, hiearchically clustered
png(file="deseq2_heatmaps_selected_orderBaseMeanCluster.png") # Create a file, it can be png, pdf..
heatmap.2(log2.norm.counts, Rowv=TRUE, Colv=TRUE, margins=c(5,10),
          lmat = lmat, lwid = lwid, lhei = lhei,
          main = paste("Top ", TOP, " significantly DE genes (log2norm)", sep=""), na.color="Green")
dev.off()

# Heatmap of DE genes, not hiearchically clustered
png(file="deseq2_heatmaps_selected_orderBaseMean.png")  # Create a file, it can be png, pdf..
heatmap.2(log2.norm.counts, dendrogram="none", Rowv=FALSE, Colv=FALSE, margins=c(5,10),
          lmat = lmat, lwid = lwid, lhei = lhei,
          main = paste("Top ", TOP, " significantly DE genes (log2norm)", sep=""))
dev.off()

TOP<-TOP_BCKP

graphics.off()

system("for i in deseq2_heatmaps_selected_*; do pdftk $i cat 2-end output tmp.pdf; mv tmp.pdf $i; done") # Cut first empty page from heatmap plot

