# generate data for grimAge

rm(list=ls())
options(stringsAsFactors=F)
# setwd('C:/Users/akelu/Dropbox/Juno_work/DNAmOline/global')
#
#your beta methyarray file
#
# mycpg='beta/dat0NoobShortMini3.csv'
load('../../../../Archived_epiAging_Risk_201808/result/prep_combine/data_gmset.swan.noFilt.rda')
beta <- getBeta(gmsetSWAN)

#
#your output filename after adding missing CpGs
#
out.csv='dat_swan.csv'
#
#PGM
#
# beta=read.csv(mycpg)
# colnames(beta)[1]='ProbeID'
ann=read.csv('grimAge_datMiniAnnotation3.csv')
cpgs=rownames(beta)
check=is.element(ann$Name,cpgs)
table(check)
#FALSE  TRUE 
# 2996 27088
miss.cpg=ann$Name[!check]
#add NA
# beta.subject=colnames(beta)[-1]
# nsubject=dim(beta)[2]-1
nmiss.cpg=length(miss.cpg)
#add subjects into the columns
add=data.frame(matrix(data=NA,nrow=nmiss.cpg,ncol=dim(beta)[2]))
colnames(add)=colnames(beta)
rownames(add)=miss.cpg
# add[,1]=miss.cpg
#
#combine
#
output=rbind(beta,add)
output=output[ann$Name,]
# output=subset(output,ProbeID %in% ann$Name)
cat('check my new beta dimension\n')
print(dim(output))
# 30084   334
dim(ann)[1]==dim(output)[1]
#T
write.csv(output, out.csv, row.names=T, quote=F)

pheno = as.data.frame(pData(gmsetSWAN))
pheno$Female = 1
pheno$Female[pheno$Gender %in% 'M'] = 0
pheno$Tissue <- 'Colon'
write.csv(pheno[,c('Tissue','Age',"Female")], 'dat_pheon.csv', quote = F)

######
# upload to http://dnamage.genetics.ucla.edu/new to calculate GrimAge






