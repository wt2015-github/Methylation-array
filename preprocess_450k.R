#!/usr/bin/env Rscript
# 2016.6.20, Ting Wang

library('limma')
library('minfi')
library('IlluminaHumanMethylation450kanno.ilmn12.hg19')
library('sva')
library('IlluminaHumanMethylation450kmanifest')
library('ENmix')
library('lumi')
library('knitr')

file.Xprobes <-'/PathTo/48639-non-specific-probes-Illumina450k.csv'
path.idat <- 'PathToIDATs'

cat('load idat files.\n')
rgSet <- read.450k.exp(base=path.idat, extended = TRUE)

cat('detection p-values.\n')
detP <- detectionP(rgSet)

cat('QC.\n')
#qc <- QCinfo(rgSet, detPthre=0.01, CpGthre=0.05, samplethre=0.05, outlier=FALSE, distplot=FALSE)
failed <- detP > 0.01
failed.probe <- rowMeans(failed) > 0.1
failed.sample <- apply(failed, 2, mean) > 0.1
if(sum(failed.sample) > 0) {
  rgSet = rgSet[,!failed.sample]
  cat(sum(failed.sample), "samples were excluded before ENmix correction\n") #1
}
Samples <- data.frame(Step='Samples', Remained=ncol(rgSet), Removed=sum(failed.sample))

cat('ENmix multimodal probe identification.\n')
beta <- getBeta(preprocessRaw(rgSet))
nmode <- nmode.mc(beta, minN = 3,modedist = 0.2, nCores = 10)
outCpG <- names(nmode)[nmode > 1]
# table(nmode)
# length(outCpG)

cat('ENmix functional normalization')
GRset <- preprocessENmix(rgSet, bgParaEst="oob", dyeCorr=TRUE, nCores=10)
beta <- getBeta(GRset)
save(beta,file="beta.FuncNorm.RData")

png('Figure_BetaDistribution_FuncNorm.png')
densityPlot(beta, main="Beta after functional normalization", xlab="Beta")
dev.off()

## remove poor quality probes
cat(sum(failed.probe), 'undetectable probes.\n')
GRset.Flt <- GRset[!failed.probe,]
QC.table <- data.frame(Step="Raw probes", Remained=dim(GRset)[1], stringsAsFactors=FALSE)
QC.table <- rbind(QC.table, c("Undetectable probes", dim(GRset.Flt)[1]))

# remove SNP probes
gmGRset.Flt <- mapToGenome(GRset.Flt)
gmGRset.Flt <- dropLociWithSnps(gmGRset.Flt, snps = c("CpG", "SBE"))
QC.table <- rbind(QC.table,c( "SNP containing probes",dim(gmGRset.Flt)[1]))

# remove cross-reactive probes
Xreact <- read.csv(file=file.Xprobes, stringsAsFactors=FALSE)
noXreact <- !(featureNames(gmGRset.Flt) %in% Xreact$TargetID) 
gmGRset.Flt <- gmGRset.Flt[noXreact,] 
QC.table <- rbind(QC.table, c("Cross-reactive probes",dim(gmGRset.Flt)[1]))

# remove sexual probes
ann.450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#write.table(ann.450k, 'Annotation_Probe.450k.txt', sep='\t', quote=F)
autosomes = !(featureNames(gmGRset.Flt) %in% ann.450k$Name[ann.450k$chr %in% c("chrX","chrY")])
gmGRset.Flt = gmGRset.Flt[autosomes,]
QC.table <- rbind(QC.table, c("Sex-chromosome probes",dim(gmGRset.Flt)[1]))

## remove any remaining multi-modal probes
notMultimodal <- !(featureNames(gmGRset.Flt) %in% outCpG)
gmGRset.Flt <- gmGRset.Flt[notMultimodal,]
QC.table <- rbind(QC.table, c("Multi-modal probes",dim(gmGRset.Flt)[1]))

## summary of number of probes removed at each quality control step:
QC.table$Removed <- c(0,-1*diff(as.numeric(QC.table$Remained),1))
QC.table <- rbind(Samples, QC.table)
#kable(QC.table, caption="Number of probes removed by the quality control steps")
write.table(QC.table, 'QC.filtering.txt', sep='\t', quote=F, row.names = F)

beta <- getBeta(gmGRset.Flt)
cat(sum(is.na(beta)),'NA betas.\n')
save(beta,file="beta.FuncNorm.filtered.RData")

png('Figure_BetaDistribution_FuncNorm.filtered.png')
densityPlot(beta, main="Beta after functional normalization and filtering", xlab="Beta")
dev.off()

