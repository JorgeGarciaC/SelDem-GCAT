## Addapted From Matt Hartfield https://doi.org/10.1002/evl3.263
#####
# SDS values standardization
#####     

rm(list = ls())

library(qvalue)
library(tidyverse)

SDSres <- read.table(paste0("~/PROJECTS.JORGE/Selection.Signals/raw_SDS_acen.txt.gz"),header=F)
names(SDSres) <- c("CHROM","ID", "AA", "DA", "POS", "DAF", "nG0", "nG1", "nG2", "rSDS", "SuggestedInitPoint")
SDSres <- SDSres[,-11]

##################################################################
SDSres <- cbind(SDSres,as.numeric(cut(SDSres[,6],seq(0.05,0.95,0.01),include.lowest = T)))
names(SDSres)[11] <- "Bin"

# Mean and SD without chr2 and 6 as in Field et al 2016
head(SDSres)
binM <- tapply(SDSres[(SDSres$CHR!="chr2" & SDSres$CHR!="chr6"),]$rSDS, SDSres[(SDSres$CHR!="chr2" & SDSres$CHR!="chr6"),]$Bin,mean) 
binSD <- tapply(SDSres[(SDSres$CHR!="chr2" & SDSres$CHR!="chr6"),]$rSDS, SDSres[(SDSres$CHR!="chr2" & SDSres$CHR!="chr6"),]$Bin,sd)

# Create column to insert values later
SDSres <- cbind(SDSres,vector(mode="numeric",length=dim(SDSres)[1]),vector(mode="numeric",length=dim(SDSres)[1]))
names(SDSres)[12] <- "sSDS"
names(SDSres)[13] <- "pSDS"

# We divide SDSraw between meand and SD by it's bin  
for(i in 1:length(binM)){
  SDSres[SDSres$Bin==i,12] <- (SDSres[SDSres$Bin==i,10]-binM[i])/binSD[i]
  SDSres[SDSres$Bin==i,13] <- pnorm(SDSres[SDSres$Bin==i,12],lower.tail=F)
}

# Multiple testing correction
qobj <- qvalue(p = SDSres$pSDS) # If pvalues does not cover whole distribution, you can use pi = 1 to use bejamini bonferroni
SDSres[,14] <- qobj$qvalues
names(SDSres)[14] <- "qSDS" 

write.csv(SDSres, "sds_norm_GCAT_acen.csv")
