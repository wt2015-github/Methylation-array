#!/usr/bin/env Rscript
#Ting

rm(list=ls())
library(cgageR)
library(minfi)
library(plyr)
library(gplots)
library(ggplot2)
library(EpiDISH)
data(centEpiFibIC.m)

load('data_grset.rda')

pheno <- as.data.frame(pData(gfilt))
beta <- getBeta(gfilt)

# estimate  cell-type fractions of total epithelial cells, total fibroblasts and total immune cells using EpiDISH with the 716 EpiFibIC DNAm reference
res.epifibic <- epidish(beta.m = beta, ref.m = centEpiFibIC.m, method = 'RPC')
res.celltype <- res.epifibic$estF

# epigenetics ages
epiAges <- getAgeR(beta, epitoc=TRUE, horvath=TRUE, hannum=TRUE, drift=FALSE, showStatusHannum=TRUE, keepcpgs.epitoc=TRUE, keepcpgs.hannum=TRUE, keepres=FALSE, chrage=NULL)
tmpTable <- read.csv('PhenoAge_513CpGmodel.csv', header = T, row.names = 1, check.names = F)
phenoAgeCpG <- tmpTable[intersect(rownames(tmpTable), rownames(beta)),]
beta.phenoAge <- beta[rownames(phenoAgeCpG),]
phenoAge <- list(PhenoAge.output = apply(beta.phenoAge, 2, function(z){sum(phenoAgeCpG$Weight * z) + tmpTable$Weight[1]}), PhenoAge.CpGs.used = rownames(phenoAgeCpG))
pmdCpG <- as.character(read.table('hm450.comPMD.probes.tsv', sep='\t', header = T, row.names = 1, check.names = F)[,1])
pmdAge <- colMeans(beta[intersect(rownames(beta), pmdCpG),])
res.ages <- data.frame(Hannum_age = epiAges$HannumClock.output$Hannum.Clock.Est[,1],
                       Horvath_age = epiAges$HorvathClock.output$Horvath.Est[,1],
                       Pheno_age = phenoAge$PhenoAge.output,
                       EpiTOC_age = epiAges$EpiTOC.output$EpiTOC.Est[,1],
                       PMD_age = pmdAge,
                       Patient = gfilt$Patient_ID, Batch=gfilt$Batch,
                       Age = gfilt$Age, Gender = gfilt$Gender,
                       Patient_DX = gfilt$Patient_DX, Sample_DX = gfilt$Sample_DX,
                       Location = gfilt$Location, Distance = gfilt$Distance,
                       Tissue_ID = gfilt$Tissue_ID, Sample_Location = gfilt$Sample_Location)
rownames(res.ages) <- colnames(gfilt)

identical(rownames(res.ages), rownames(res.celltype)) #T
res.ages <- cbind(res.ages, res.celltype, pheno)
write.table(res.ages, file='result_ages.txt', sep='\t', quote = F)
