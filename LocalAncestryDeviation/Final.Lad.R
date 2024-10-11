rm(list=ls())
library(tidyverse)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Gviz)
library(ggplotify)
library(org.Hs.eg.db)
library(Homo.sapiens)
library(MoMAColors)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
standard.d <- 3 #4.42
color.palette <- c("#744940","#e6bcac","#3d3d47")

##### Functions

plot.percentage.deviation <- function(df,POP1,POP2,POP3,Momma_scale = "Picabia", target) {
  # Function used to plot the genome wide percentage deviation for component
  df <- pivot_longer(df,cols = c(POP1,POP2,POP3))
  
  df <- df %>% group_by(name) %>% mutate(standard_deviation = mean(value) + standard.d * sd(value), mean = mean(value))
  
  don <- df %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(epos)) %>% # Positions in BP
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, epos) %>%
    mutate(BPcum=epos+tot)
  
  # Then we need to prepare the X axis. 
  axisdf = don %>% 
    group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  p2 <- ggplot(don, aes(BPcum,value,fill=name)) + geom_area() +
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) + 
    scale_fill_manual(values = color.palette, labels =c(POP1,POP2,POP3)) +
    theme_bw() +
    theme(
      #panel.border = element_blank(),
      axis.text = element_text(size=14),
      axis.title = element_text(size=18),
      strip.text = element_text(size=18),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(), 
      legend.position = "none"
    ) +
    labs(x = "Chromosome", y = paste0("Proportion in ",target), fill = "Ancestry")  
  p3 <- p2 + facet_grid( name ~ .) + 
    geom_hline(aes(yintercept=standard_deviation), colour="red", lwd=0.5)  +
    geom_hline(aes(yintercept=mean), colour="black", lwd=0.5)
  return(p3)
}

make.viz.plot <- function(df, seqn, starting.point,ending.point, txObject) {
  # Function to make viz plots
  df$seqnames <- seqn
  df$seqnames <- as.factor(seqn)
  
  itrack <- IdeogramTrack(genome = "hg19", chromosome = seqn)
  txTr <- GeneRegionTrack(txObject, chromosome = seqn, start = starting.point, end = ending.point, geneSymbol = T, showId = T, gene = T, name = NULL)
  z <- ranges(txTr)
  z$symbol <- mapIds(Homo.sapiens, z$symbol, "SYMBOL", "TXNAME")
  ranges(txTr) <- z
  axisTrack <- GenomeAxisTrack()
  dTrack <- DataTrack(df, genome = "hg19", name = "-log(qSDS)", data = -log(df$qSDS), chromosome = seqn, baseline=-log(0.05), col.baseline="red")
  
  plot.gviz <-
    plotTracks(
      c(dTrack, txTr, axisTrack),
      from = starting.point, to = ending.point,
      collapseTranscripts = "longest", shape = "arrow",
      transcriptAnnotation = "symbol", sizes = c(3, 3, 1),
      background.title = "white",
      background.panel = "white",
      col.axis = "black",
      col.title = "black",
      cex.title = 0.9,       # Reducir el tamaño del texto del título
      col.frame = NA,        # Eliminar los bordes del marco
      margin = 1,            # Reducir el margen alrededor de las pistas
      title.width = 0.65     # Reducir el ancho del panel del título
    )
  
  if (is.null(plot.gviz)) {
    stop("Error while creating the graph. Check")
  }
  
  return(plot.gviz)
}

# Functions to make GRanges objects
sd2GRanges <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = df[,2], end = df[,3]))
  return(gr)
}

sd2GRanges.metadata <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = df[,2], end = df[,3]),
    peak = df[,4],
    qSDS = df[,5])
  return(gr)
}

# Function to calculate ancestry population from vtb files per position for a population
calculate.ancestry.population <- function(vtb.file,P1,P2,P3) {
  hap_1 <- cbind(vtb.file[,1:6],vtb.file[,seq(7,dim(vtb.file)[2],2)])
  hap_2 <- cbind(vtb.file[,1:6],vtb.file[,seq(8,dim(vtb.file)[2],2)])
  colnames(hap_1) <- c("CHR","spos","epos","sgpos","egpos","n_snps")
  colnames(hap_2) <- c("CHR","spos","epos","sgpos","egpos","n_snps")
  
  POP1_1 <- hap_1[7:dim(hap_1)[2]] %>% apply(1, function(x) length(which(x=="0"))/dim(.)[2])
  POP2_1 <- hap_1[7:dim(hap_1)[2]] %>% apply(1, function(x) length(which(x=="1"))/dim(.)[2])
  POP3_1 <- hap_1[7:dim(hap_1)[2]] %>% apply(1, function(x) length(which(x=="2"))/dim(.)[2])
  
  POP1_2 <- hap_2[7:dim(hap_2)[2]] %>% apply(1, function(x) length(which(x=="0"))/dim(.)[2])
  POP2_2 <- hap_2[7:dim(hap_2)[2]] %>% apply(1, function(x) length(which(x=="1"))/dim(.)[2])
  POP3_2 <- hap_2[7:dim(hap_2)[2]] %>% apply(1, function(x) length(which(x=="2"))/dim(.)[2])
  
  POP1 <- rowMeans(cbind(POP1_1,POP1_2))
  POP2 <- rowMeans(cbind(POP2_1,POP2_2))
  POP3 <- rowMeans(cbind(POP3_1,POP3_2))
  
  df_for_plot <- cbind(hap_2[1:3],POP1,POP2,POP3)
  colnames(df_for_plot)[c(4,5,6)] <- c(P1,P2,P3)
  return(df_for_plot)
}

# Function to make the zooms of ancestry for a chromosome
plot.chromosome.deviation <- function(df,POP1,POP2,POP3,Momma_scale = "Picabia", chromosome, standard = standard.d ) {
  df <- pivot_longer(df,cols = c(POP1,POP2,POP3))
  
  df <- df %>% group_by(name) %>% mutate(standard_deviation = mean(value) + standard * sd(value), mean = mean(value))
  
  p2 <- ggplot(df[df$CHR == chromosome,], aes(spos,value,fill=name)) + geom_area() +
    #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    #scale_y_continuous(expand = c(0, 0) ) + 
    scale_fill_moma_d(Momma_scale, labels = c(POP1,POP2,POP3)) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(), 
      legend.position = "none"
    ) +
    labs(x = "Chromosome", y = "Proportion", fill = "Ancestry")  
  p3 <- p2 + facet_grid( name ~ .) + 
    geom_hline(aes(yintercept=standard_deviation), colour="red", lwd=0.5) +
    geom_hline(aes(yintercept=mean), colour="black", lwd=0.5)
  return(p3)
}

# Function to plot only one population
plot.chromosome.onepop.deviation <- function(df,POP, color.used, chromosome, target,standard = standard.d) {
  df <- pivot_longer(df,cols = POP)
  
  df <- df %>% group_by(name) %>% mutate(standard_deviation = mean(value) + standard * sd(value), mean = mean(value))
  
  p2 <- ggplot(df[df$CHR == chromosome,], aes(spos,value,fill=name)) + geom_area() +
    scale_x_continuous(
      labels = scales::label_scientific(digits = 2)
      ,breaks= seq(0,max(df$spos),by= max(df$spos)/10) 
      ) +
    #scale_y_continuous(expand = c(0, 0) ) + 
    scale_fill_manual(values = color.used) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      legend.position = "none"
    ) +
    labs(x = paste0("Chromosome ", chromosome), y = paste0("Proportion of ", POP, " ancestry in ",target), fill = "Ancestry")  +
    geom_hline(aes(yintercept=standard_deviation), colour="red", lwd=0.5) +
    geom_hline(aes(yintercept=mean), colour="black", lwd=0.5)
  return(p2)
}

# Get transcript for genes to plot for the MHC because 
# it has too many genes and too many transcripts in
# different patches to be nicely plotted
get.transcript.for.gene <- function(gene.symbol) {
  key <- select(org.Hs.eg.db, keys = gene.symbol, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  gene <- select(TxDb.Hsapiens.UCSC.hg19.knownGene, keys = key[,2], keytype = "GENEID", columns = c("CDSCHROM","CDSEND","CDSSTART"))
  colnames(gene) <- c("symbol", "chromosome","start","end")
  seqn <- "chr6"
  gene <- gene[gene$chromosome == seqn,]
  txTr <- GeneRegionTrack(gene, chromosome = seqn, start = 32e6, end = 32.3e6, geneSymbol = T, showId = T, gene = T, name = NULL)
  z <- ranges(txTr)
  z$symbol <- mapIds(Homo.sapiens, z$symbol, "SYMBOL", "ENTREZID")
  ranges(txTr) <- z
  return(txTr)
}

#### Plot lad
setwd("~/PROJECTS.JORGE/LAD/Zoom.Lad/")

SDS <- read_csv("~/PROJECTS.JORGE/Selection.Signals/SDS_ispeak.csv")

vtb <- read.table(file = "All.msp.ancient.tsv") #EEF 0, WHG 1, ENS 2

to_plot <- calculate.ancestry.population(vtb.file = vtb, "EEF", "WHG","ENS")
p.ancient.gcat <- plot.percentage.deviation(to_plot,"WHG","EEF","ENS", target = "GCAT")

WHG_d <- to_plot[to_plot$WHG > mean(to_plot$WHG) + standard.d *sd(to_plot$WHG),] 
EEF_d <- to_plot[to_plot$EEF > mean(to_plot$EEF) + standard.d *sd(to_plot$EEF),] 
ENS_d <- to_plot[to_plot$ENS > mean(to_plot$ENS) + standard.d *sd(to_plot$ENS),]

EEF.sd <- GenomicRanges::reduce(sd2GRanges(EEF_d))
WHG.sd <- GenomicRanges::reduce(sd2GRanges(WHG_d))
ENS.sd <- GenomicRanges::reduce(sd2GRanges(ENS_d))

sdsGR <- sd2GRanges.metadata(as.data.frame(SDS[,c(16,6,6,17,15)]))

hts <- findOverlaps(ENS.sd, sdsGR)
sds.ENS <- sdsGR[subjectHits(hts)]

LCT <- as.data.frame(sds.ENS[sds.ENS@seqnames == 2,])
MHC <- as.data.frame(sds.ENS[sds.ENS@seqnames == 6,])

hts <- findOverlaps(WHG.sd, sdsGR)
sds.WHG <- sdsGR[subjectHits(hts)]

WHG <- as.data.frame(sds.WHG[sds.WHG@seqnames == 19,])

#### Plot tracks

WHG.plot <- as.ggplot(~make.viz.plot(WHG,"chr19",starting.point = 21.97e6, ending.point = 22.35e6,txdb))
LCT.plot <- as.ggplot(~make.viz.plot(LCT,"chr2",starting.point = 135e6, ending.point = 138e6,txdb))

seqn <- "chr6"
MHC$seqnames <- as.factor("chr6")
itrack <- IdeogramTrack(genome = "hg19", chromosome = seqn)
axisTrack <- GenomeAxisTrack()
dTrack <- DataTrack(MHC, genome = "hg19", name = "-log(qSDS)", data = -log(MHC$qSDS), chromosome = seqn, baseline=-log(0.05), col.baseline="red")

PBX2.t <- get.transcript.for.gene("PBX2")
RNF5.t <- get.transcript.for.gene("RNF5")
NOTCH4.t <- get.transcript.for.gene("NOTCH4")
AGER.t <- get.transcript.for.gene("AGER")
AGPAT1.t <- get.transcript.for.gene("AGPAT1")
GPSM3.t <- get.transcript.for.gene("GPSM3")
TSBP1.t <- get.transcript.for.gene("TSBP1")
PRRT1.t <- get.transcript.for.gene("PRRT1")
PPT2.t <- get.transcript.for.gene("PPT2")
EGFL8.t <- get.transcript.for.gene("EGFL8")
FKBPL.t <- get.transcript.for.gene("FKBPL")
ATF6B.t <- get.transcript.for.gene("ATF6B")
TNXB.t <- get.transcript.for.gene("TNXB")

#
MHC.plot <- as.ggplot(~plotTracks(
    c(dTrack, PBX2.t, RNF5.t, NOTCH4.t, AGER.t, AGPAT1.t, GPSM3.t,TSBP1.t,PRRT1.t,PPT2.t,EGFL8.t,FKBPL.t,ATF6B.t,TNXB.t, axisTrack),
    from = 32e6, to = 32.3e6,
    collapseTranscripts = "gene", shape = "arrow",
    transcriptAnnotation = "symbol",
    background.title = "white",
    background.panel = "white",
    col.axis = "black",
    col.title = "black",
    cex.title = 0.9,       # Reducir el tamaño del texto del título
    col.frame = NA,        # Eliminar los bordes del marco
    margin = 1,            # Reducir el margen alrededor de las pistas
    title.width = 0.65,     # Reducir el ancho del panel del título
    sizes= c(9,1,1,1,1,1,1,1,1,1,1,1,1,1,2)  
))

# MODERN LAD

vtb.downsample <- read.table(file = "all.gcat.downsamplesud.msp.tsv") #ENA=0 PAL=1 SUD=2
to_plot.2 <- calculate.ancestry.population(vtb.downsample,"ENA","FRANCE","WNA")

# PAL ancient LAD
vtb.pal.ancient <- read.table("paint.ALL.PAL.msp.tsv") #EEF 0, WHG 1, ENS 2

PAL.ancient.df <- calculate.ancestry.population(vtb.pal.ancient,"EEF","WHG","ENS")
PAL.ancient.df$CHR <- as.numeric(gsub("chr","",PAL.ancient.df$CHR)) # as is in GRCH38 needs to change chr names

# NA ancient LAD
vtb.naf.ancient <- read.table("NorthAfrica.Ancient.msp.tsv") #EEF, WHG,1 ENS2
naf.to_plot <- calculate.ancestry.population(vtb.naf.ancient,"EEF","WHG","ENS")

# PAL and North Africa modern on GCAT
vtb.palnaf.modern <- read.table("PAL_NAF.modern.msp.tsv") #ENA, PAL, SUD
pal.naf.to_plot <- calculate.ancestry.population(vtb.palnaf.modern, "NAF","PAL","FRANCE")

#### Individual chromosomes plots
deviated.region.ENS.MHC <- plot.chromosome.onepop.deviation(to_plot,"ENS", chromosome = 6, color.used =  color.palette[2], "GCAT" )
deviated.region.ENS.LCT <- plot.chromosome.onepop.deviation(to_plot,"ENS", chromosome = 2, color.used =  color.palette[2], "GCAT")
deviated.region.EEF  <- plot.chromosome.onepop.deviation(to_plot,"EEF", chromosome = 11, color.used =  color.palette[1], "GCAT" )
deviated.region.WHG.11  <- plot.chromosome.onepop.deviation(to_plot,"WHG", chromosome = 11, color.used =    color.palette[3], "GCAT")
deviated.region.WHG.14 <-  plot.chromosome.onepop.deviation(to_plot,"WHG", chromosome = 14, color.used =    color.palette[3], "GCAT")
deviated.region.WHG.19 <-  plot.chromosome.onepop.deviation(to_plot,"WHG", chromosome = 19, color.used =    color.palette[3], "GCAT")
deviated.region.WHG.21 <-  plot.chromosome.onepop.deviation(to_plot,"WHG", chromosome = 21, color.used =    color.palette[3], "GCAT")

deviated.regions.NAF.modern <- plot.chromosome.onepop.deviation(pal.naf.to_plot, chromosome = 6, "NAF", color.used = "lightgrey", "GCAT", standard = 4.42)
deviated.regions.PAL.modern <- plot.chromosome.onepop.deviation(pal.naf.to_plot, chromosome = 6, "PAL", color.used = "darkgrey", "GCAT", standard = 4.42)
deviated.region.EEF.NA.ANCIENT <- plot.chromosome.onepop.deviation(naf.to_plot, chromosome = 6, "ENS", color.used = color.palette[2], "NA")
deviated.region.EEF.PAL.ANCIENT <- plot.chromosome.onepop.deviation(PAL.ancient.df, chromosome = 6, "ENS", color.used = color.palette[2], "PAL")

# Figure 1 Ancestral proportions in GCAT

p.ancient.gcat
ggsave("~/PROJECTS.JORGE/GCAT.SUPPTABLES/FINAL FIGURES/ancient.LAD.GCAT.pdf", device = "pdf", plot =  last_plot(), width = 24, height = 12)


# Figure 2, zoom LCT y ZNF region WHG
cowplot::plot_grid(deviated.region.ENS.LCT, LCT.plot, deviated.region.WHG.19, WHG.plot, rel_widths = c(1,1.5,1,1.5), labels = c("A","C","B","D"))
ggsave("~/PROJECTS.JORGE/GCAT.SUPPTABLES/FINAL FIGURES/LCT_WHG.figure.pdf", device = "pdf", plot =  last_plot(), width = 18, height = 9)

# Figure 3
p.MODERN <- cowplot::plot_grid(deviated.regions.NAF.modern,deviated.regions.PAL.modern, ncol = 1, labels = c("A","B"))
p.ANCIENT <- cowplot::plot_grid(deviated.region.EEF.NA.ANCIENT,deviated.region.EEF.PAL.ANCIENT,deviated.region.ENS.MHC,nrow = 1, labels = c("C","D","E"))
p.COMPOSE <- cowplot::plot_grid(p.ANCIENT, MHC.plot,ncol = 1, labels = c("","F"), rel_heights = c(1,1.75))
cowplot::plot_grid(p.MODERN,p.COMPOSE, nrow=1, rel_widths = c(1,1.5))
ggsave("~/PROJECTS.JORGE/GCAT.SUPPTABLES/FINAL FIGURES/MHC.figure.pdf", device = "pdf", plot =  last_plot(), width = 18, height = 9)

# Supp figure with all ampliations
cowplot::plot_grid(deviated.region.ENS.MHC,deviated.region.ENS.LCT,deviated.region.EEF,
                   deviated.region.WHG.11,deviated.region.WHG.14,deviated.region.WHG.19,
                   deviated.region.WHG.21)
ggsave("~/PROJECTS.JORGE/GCAT.SUPPTABLES/FINAL FIGURES/ALL.zooms.figure.pdf", device = "pdf", plot =  last_plot(), width = 18, height = 9)

# More supp
colnames(to_plot.2)[5] <- "SUD"
standard.d <- 4.42
plot.percentage.deviation(to_plot.2, "ENA","SUD","WNA", target = "GCAT")
ggsave("~/PROJECTS.JORGE/GCAT.SUPPTABLES/FINAL FIGURES/downsample_SUD.LAD.GCAT.pdf", device = "pdf", plot =  last_plot(), width = 24, height = 12)

#"#53362e""#744940" "#9f7064" "#c99582" "#e6bcac" "#e2d8d6" "#a5a6ae" "#858794" "#666879" "#515260" "#3d3d47"
