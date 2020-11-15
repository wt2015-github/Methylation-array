#!/usr/bin/env Rscript

# devtools::install_github("metamaden/methyPre")
library(methyPre)
library(minfi)
library(ggplot2)

######
# load batch 1 raw idat files
rg1 <- read.metharray.exp(base='path-to-ImageData-IDATs-1', recursive=T, force=T, verbose = T)
# add phenotype pData
pheno <- read.csv('phenotyoe1.csv', header = T, check.names = F)
rownames(pheno) <- pheno$`EPIC Array ID`
pheno <- pheno[colnames(rg1),]
pData(rg1) <- cbind(pData(rg1), pheno)
pData(rg1)$Batch <- 'Batch1'
# save(rg1, file='data_rgSet_batch1.rda')

######
# load batch 2 raw idat files
rg2 <- read.metharray.exp(base='path-to-ImageData-IDATs-2', recursive=T, force=T, verbose = T)
# add phenotype pData
pheno <- read.csv('phenotyoe2.csv', header = T, check.names = F)
rownames(pheno) <- pheno$`EPIC Array ID`
pheno <- pheno[colnames(rg2),]
pData(rg2) <- cbind(pData(rg2), pheno)
pData(rg2)$Batch <- 'Batch2'
# save(rg2, file='data_rgSet_batch2.rda')

############
# combine two batches
# load('data_rgSet_batch1.rda')
# load('data_rgSet_batch2.rda')
identical(colnames(pData(rg1)),colnames(pData(rg2))) #T
rg.all <- combineArrays(rg1, rg2, outType = 'IlluminaHumanMethylationEPIC') #or HM450
dim(rg.all) #1051815 88

# PCA of snp probes
snpinfo <- getSnpInfo(rg.all)
snps <- getSnpBeta(rg.all)
respca2 <- prcomp(t(snps))$x
# pdf('figure_pca_snpMethy.pdf', height = 5, width = 5)
ggplot(data.frame(PC1=respca2[,1], PC2=respca2[,2], Batch=rg.all$Batch), aes(x = PC1, y = PC2, color=Batch))+
  ggtitle("PCA on DNAm of SNP probes")+
  geom_point()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
# dev.off()
# save(snps,file='data_snpMethy.rda')

# QC
mset <- preprocessRaw(rg.all)
qc <- getQC(mset)
pdf('figure_qc.pdf', height = 5, width = 10)
par(mfrow=c(1,2))
plotQC(qc)
densityPlot(mset)
dev.off()
# pdf('figure_qc_bs.pdf', height = 9, width = 7)
# par(mfrow=c(1,2))
controlStripPlot(rg.all, controls="BISULFITE CONVERSION I")
controlStripPlot(rg.all, controls="BISULFITE CONVERSION II")
# dev.off()

# Illumina and SWAN normalization
datIlmn <- preprocessIllumina(rg.all)
msetSWAN <- preprocessSWAN(rg.all, mSet = datIlmn)
dim(msetSWAN) #

detP <- detectionP(rg.all)
msetSWAN <- msetSWAN[rownames(detP),]

# filter CpGs by detectable p-values
identical(rownames(msetSWAN), rownames(detP))
identical(colnames(msetSWAN), colnames(detP))
detPcutoff <- 0.01
detPpercent <- 0.1
max(colMeans(detP > detPcutoff)) #
failedCpG <- rownames(detP)[rowMeans(detP > detPcutoff) > detPpercent]
cat(paste("Remove",length(failedCpG), "undetected CpGs\n")) #
gfilt <- msetSWAN[setdiff(rownames(msetSWAN),failedCpG),]
dim(gfilt) #

# map to genome
gfilt <- mapToGenome(gfilt)

# apply minfi filters
cat("Removing cg probes containing SNPs...\n")
gfilt <- dropLociWithSnps(gfilt, snps = c("SBE", "CpG"))
cat("Removing CH and SNP-assoc. probes...\n")
gfilt <- dropMethylationLoci(gfilt, dropRS = T, dropCH = T)
cat("Removing chrX and chrY-assoc. probes...\n")
gfilt <- gfilt[!getAnnotation(gfilt)$chr %in% c("chrX","chrY"),] 
cat("After applying minfi filters,", nrow(gfilt), "CpGs remain\n") #812017
cat("Filtering cross-reactive EPIC probes\n")
data(illumina_crxcg)
data(pidsley_crxcg)
gfilt <- gfilt[!rownames(gfilt) %in% c(illumina.crxcg, pidsley.crxcg),]
# if it's HM450 data
# data(chen_crxcg)
# gfilt.all <- gfilt.all[!rownames(gfilt.all) %in% chen.crxcg,]
dim(gfilt) #
save(gfilt, file='data_gmSet_combined_filtered.rda')



