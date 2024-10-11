###### Plot xpehh #####


rm(list=ls())

library(qvalue)
library(tidyverse)
library(GenomicRanges)
library(ggrepel)
library(httr)
library(jsonlite)
library(xml2)

args = commandArgs(trailingOnly=TRUE)

xp_in <- args[1]
csv_out <- args[2]

setwd("~/xpehh/")
xpehh <- read.table(xp_in, header = T)

xpehh$peak <- FALSE
percentiles_above <- quantile(xpehh$normxpehh, c(0.999, 0.9995)) # delimit the genome-wide 95 and 99%
print(percentiles_above)
percentiles_beneath <- quantile(xpehh$normxpehh, c(0.001, 0.0005))
print(percentiles_beneath)
new_df <- as.data.frame(NULL)

pos_bound <- 500000

xpehh$n_peak <- F

# We have to look not only for the top ones, also for the bottom ones

for (chrom in 1:22){ # I extract each chromosome
  n <- 0
  last_peak <- F
  
  xp_aux <- xpehh[xpehh$chr==chrom,]
  print(paste0("Chrom ", chrom, " ", dim(xp_aux)[1]))
  
  for (i in 1:dim(xp_aux)[1]){ # and iterate over each SNP
    
    if (xp_aux$normxpehh[i] > percentiles_above[[1]]) { # if the SNP is above the 95%
      lowb <- xp_aux[i,"pos"]-pos_bound # set lower bound
      highb <- xp_aux[i,"pos"]+ pos_bound # set upper bound
      num_top <- sum(xp_aux[xp_aux$pos >= lowb & 
                              xp_aux$pos <= highb,]$normxpehh 
                     >= percentiles_above[[2]]) # count how many SNPs in the defined region are above the 99%

      if (num_top>=10) { # I only consider the SNP to be in a peak in case there are at least 3 SNPs in the top 0.05%
        xp_aux$peak[i] <- TRUE
        if (last_peak == F) {
          n <- n+1}
        last_peak <- T
        xp_aux$n_peak[i] <- n
      } else {
        last_peak <- F
      }
   }
    else if (xp_aux$normxpehh[i] < percentiles_beneath[[1]]) {
      lowb <- xp_aux[i,"pos"]-pos_bound
      highb <- xp_aux[i,"pos"]+pos_bound
      num_top <- sum(xp_aux[xp_aux$pos >= lowb &
                              xp_aux$pos <= highb,]$normxpehh
                     <= percentiles_beneath[[2]])
      if (num_top>=10) {
        xp_aux$peak[i] <- TRUE
        if (last_peak == F) {
          n <- n+1}
        last_peak <- T
        xp_aux$n_peak[i] <- n
      } else {
        last_peak <- F
      }
      }
    }   
  new_df <- rbind(xp_aux,new_df)
}

query <- new_df[new_df$peak == T,]$id

write.csv(new_df, paste(csv_out, "_peaks.csv", sep = ''))
write.table(query, paste(csv_out, "_vep.txt", sep = ''), row.names = F, col.names = F,)
