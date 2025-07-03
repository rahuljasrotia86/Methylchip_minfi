# Methylchip_minfi
Epigenetic-mediated transcriptome analysis in primary Sjogren’s disease reveals new cellular pathways

# Load Required Libraries
library("minfi")
library("minfiData")
library("limma")
library("RColorBrewer")
library("missMethyl")
library("Gviz")
library("DMRcate")
library("DMRcatedata")
library("stringr")
library("mCSEA")
library("tidyverse")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(doParallel)
registerDoParallel(cores = 4)

# Load Annotation Data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

# Load Raw IDAT Files and Sample Sheet
dataDirectory <- "extdata"
list.files(dataDirectory, recursive = TRUE)

targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
RGSet <- read.metharray.exp(targets=targets)

# Quality Control (QC)
MSet <- preprocessRaw(RGSet)
qc <- getQC(MSet)
plotQC(qc)

detP <- detectionP(RGSet)
barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05, col="red")

densityPlot(MSet, sampGroups = pData(MSet)$Sample_Group)
controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")


# Normalization
mSetSq <- preprocessQuantile(RGSet)

# Multi-Dimensional Scaling (MDS) Plots
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, col=brewer.pal(8,"Dark2")[factor(targets$Sample_Group)])


# Probe Filtering
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)


# Differential Methylation: CpG-Level (DMPs)
mVals <- getM(mSetSqFlt)
design <- model.matrix(~0+factor(targets$Sample_Group)+factor(targets$Sample_Source))
fit <- lmFit(mVals, design)
contMatrix <- makeContrasts(Tre-Con, levels=design)
fit2 <- eBayes(contrasts.fit(fit, contMatrix))
DMPs <- topTable(fit2, num=Inf, coef=1)


# Export DMP Results
write.table(DMPs, file="new_DMPs.csv", sep=",", row.names=FALSE)

# Visualize Top DMPs
bVals <- getBeta(mSetSqFlt)
plotCpg(bVals, cpg="cg07499259", pheno=targets$Sample_Group, ylab = "Beta values")


# Regional Differential Methylation (DMRs)
myAnnotation <- cpg.annotate(object = mVals, ...)
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
write.table(results.ranges, file="new_results_ranges.csv", sep=",", row.names=FALSE)


