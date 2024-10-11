##### How to run it: Rscript filter_and_plot.R SDS_file

rm(list=ls())
library("ggplot2")
library(tidyverse)

##############
# Read files
#############

#setwd("~/Jorge/SDS/")
plot_name <- "prueba"
args <- commandArgs(TRUE)
#SDS <- read.delim(args[1], header=TRUE) # read the SDS file
SDS <- read.csv(args[1], header=TRUE) # read the SDS file
#SDS <- read.csv("~/Jorge/SDS/sds_norm_GCAT_acen.csv", header=TRUE) # read the SDS file
plot_name <- strsplit(args[1], ".txt")[[1]] # keep the input filename to generate the output filename

SDS <- SDS[,2:15]
######################
# Peak definition
#####################
percentiles <- quantile(SDS$pSDS, c(0.001, 0.0005)) # delimit the genome-wide 95 and 99%
print(percentiles)

SDS$CHROM <- as.numeric(substr(as.character(SDS$CHR),4,5))
SDS$peak <- FALSE # initialize the peak variable -- TRUE will indicate those SNPs that belong to a peak

new_df <- as.data.frame(NULL)
SDS$n_peak <- F
bp_bound <- 500000

for (chrom in 1:22){ # I extract each chromosome

	SDS_aux <- SDS[SDS$CHROM==chrom,]
	print(paste0("Chrom ", chrom, " ", dim(SDS_aux)[1]))
  n <- 0
  last_peak <- F
	for (i in 1:dim(SDS_aux)[1]){ # and iterate over each SNP
 
		if (SDS_aux$pSDS[i] < percentiles[[1]]) { # if the SNP is above the 95%
			lowb <- SDS_aux[i,"POS"]-bp_bound # set lower bound
			highb <- SDS_aux[i,"POS"]+ bp_bound # set upper bound
			num_top <- sum(SDS_aux[SDS_aux$POS >= lowb & SDS_aux$POS <= highb,]$pSDS <= percentiles[[2]]) # count how many SNPs in the defined region are above the 99%
      
			if (num_top>=10) { # I only consider the SNP to be in a peak in case there are at least 3 SNPs in the top 0.05%
				SDS_aux$peak[i] <- TRUE
				if (last_peak == F) {
				  n <- n+1}
				SDS_aux$n_peak[i] <- n
				last_peak <- T
			} else {
			  last_peak <- F
			}
		}
	}
  
  new_df <- rbind(SDS_aux,new_df)
}

write_csv(new_df, "SDS_ispeak.csv")


#### PLOTTING
don <- new_df %>% 
  
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% # Positions in BP
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(new_df, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate(BPcum=POS+tot)


# Then we need to prepare the X axis. 
# Indeed we do not want to display the cumulative position of SNP in bp, 
# but just show the chromosome name instead

axisdf = don %>% 
  group_by(CHROM) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


# plot 
p <- ggplot(don, aes(x=BPcum, y=-log10(pSDS))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  geom_point(data=subset(don, peak==T), color="lightsalmon", size=1)  

plot_name <- "SDS"
png(paste0(plot_name, "_manhattan.png"), width = 1280, height = 720)
print(p)
dev.off()