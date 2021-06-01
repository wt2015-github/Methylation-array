#!/usr/bin/env Rscript
#Ting, 2020.10.2

rm(list=ls())
# devtools::install_github('metamaden/cgageR')
library(cgageR)
library(minfi)
library(plyr)
library(gplots)
library(ggplot2)
library(EpiDISH)
data(centEpiFibIC.m)

load('../prep/data_gmset.swan.noFilt.rda')
pheno <- data.frame(pData(gmsetSWAN), check.names = F)
beta <- getBeta(gmsetSWAN)

# estimate  cell-type fractions of total epithelial cells, total fibroblasts and total immune cells using EpiDISH with the 716 EpiFibIC DNAm reference
res.epifibic <- epidish(beta.m = beta, ref.m = centEpiFibIC.m, method = 'RPC')
res.celltype <- res.epifibic$estF

# epigenetics ages
epiAges <- getAgeR(beta, epitoc=TRUE, horvath=TRUE, hannum=TRUE, drift=FALSE, showStatusHannum=TRUE, keepcpgs.epitoc=TRUE, keepcpgs.hannum=TRUE, keepres=FALSE, chrage=NULL)
tmpTable <- read.csv('./doc/PhenoAge_513CpGmodel.csv', header = T, row.names = 1, check.names = F)
phenoAgeCpG <- tmpTable[intersect(rownames(tmpTable), rownames(beta)),]
beta.phenoAge <- beta[rownames(phenoAgeCpG),]
phenoAge <- list(PhenoAge.output = apply(beta.phenoAge, 2, function(z){sum(phenoAgeCpG$Weight * z) + tmpTable$Weight[1]}), PhenoAge.CpGs.used = rownames(phenoAgeCpG))
source('./doc/epiTOC2.R')
outEpiTOC2 <- epiTOC2('./doc/epiTOC2_data.Rd', data.m = beta, ages.v=pheno$Age)
res.ages <- data.frame(Hannum_age = epiAges$HannumClock.output$Hannum.Clock.Est[,1],
                       Horvath_age = epiAges$HorvathClock.output$Horvath.Est[,1],
                       Pheno_age = phenoAge$PhenoAge.output,
                       EpiTOC_age = epiAges$EpiTOC.output$EpiTOC.Est[,1],
                       EpiTOC2_age = outEpiTOC2$tnsc,
                       HypoClock_age = outEpiTOC2$hypoSC)
rownames(res.ages) <- colnames(beta)

identical(rownames(res.ages), rownames(res.celltype)) #T
res.ages <- cbind(res.ages, res.celltype, pheno)
write.table(res.ages, file='result_ages.txt', sep='\t', quote = F)

