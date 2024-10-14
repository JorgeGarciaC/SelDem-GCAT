###### PCA 07/06/24 snps = 141849

rm(list=ls())

setwd("~/PROJECTS.JORGE/PCA.DEMOGRAPHY/General.Plots.Lazar/")

library(tidyverse)
library(ggalt)

#### Define plot function 
do.changes.to.data <- function (eigenvec, eigenval,labs) {
  colnames(eigenvec)[1] <- "IID"
  colnames(labs) <- c("IID","POP")
  
  to.plot <- merge(eigenvec,labs,by="IID", sort = F)
  colnames(to.plot)[c(2,3)] <- c("PC1","PC2")
  
  to.plot$SUPERPOP <- "EUROPE"
  
  F_s <- c("French_Basque","French","PAR","BR","IeV","PdD","NO","BdR","HG")
  to.plot[to.plot$POP %in% F_s,]$SUPERPOP <- "France"
  
  to.plot[to.plot$POP == "BR",]$POP <- "Alsace"
  to.plot[to.plot$POP == "BdR",]$POP <- "Provence"
  to.plot[to.plot$POP == "IeV",]$POP <- "Britanny"
  to.plot[to.plot$POP == "PdD",]$POP <- "Dordogne"
  to.plot[to.plot$POP == "NO",]$POP <- "Nord"
  to.plot[to.plot$POP == "PAR",]$POP <- "Paris"
  to.plot[to.plot$POP == "HG",]$POP <- "Pyrenees"
  
  GCAT_s <- c("ANDALUCIA","CATALUÑA","MADRID","EXTREMADURA","ARAGON","CYL","CLM","MURCIA","VALENCIANA","GALICIA","MADRID")
  to.plot[to.plot$POP %in% GCAT_s,]$SUPERPOP <- "GCAT"
  
  CAT_S <- c("València", "Catalonia", "Balearic_Is.")
  to.plot[to.plot$POP %in% CAT_S,]$SUPERPOP <- "People from Ibiza"
  return(to.plot)
} 


##### CHANGES 16 Septiembre Figure s4
shapes <- c(1,16,16,16,16,16,16,16,16,2,16,16,16,3,4,5,6,7,16,16)
plot.pca <- function (full.df, subseted.df, color.column, polygon.column) {
  p <- ggplot(full.df, aes(PC1, PC2, col={{color.column}}, shape = {{color.column}})) + geom_point(size = 1.5) + 
    coord_equal() + theme_minimal() + theme(panel.grid.minor = element_blank()) + 
    theme(panel.grid.major = element_blank()) +  scale_shape_manual(values= shapes) + 
    scale_color_discrete(type = c_palette) +
    xlab(paste0("PC1 (", signif(eigenval[1]*100, 3), "%)")) + 
    ylab(paste0("PC2 (", signif(eigenval[2]*100, 3), "%)")) 
  p <- p +
    geom_encircle(data = subseted.df, aes(group = {{polygon.column}}, fill={{polygon.column}}), s_shape = 1, expand = 0,
                  alpha = 0.2, color = "black", show.legend = FALSE, na.rm = T)
  return(p)
}


do.changes.to.data <- function (eigenvec, eigenval,labs) {
  colnames(eigenvec)[1] <- "IID"
  colnames(labs) <- c("IID","POP")
  
  to.plot <- merge(eigenvec,labs,by="IID", sort = F)
  colnames(to.plot)[c(2,3)] <- c("PC1","PC2")
  
  to.plot$SUPERPOP <- "EUROPE"
  
  F_s <- c("French_Basque","French","PAR","BR","IeV","PdD","NO","BdR","HG")
  to.plot[to.plot$POP %in% F_s,]$SUPERPOP <- "France"
  
  to.plot[to.plot$POP == "BR",]$POP <- "Alsace"
  to.plot[to.plot$POP == "BdR",]$POP <- "Provence"
  to.plot[to.plot$POP == "IeV",]$POP <- "Britanny"
  to.plot[to.plot$POP == "PdD",]$POP <- "Dordogne"
  to.plot[to.plot$POP == "NO",]$POP <- "Nord"
  to.plot[to.plot$POP == "PAR",]$POP <- "Paris"
  to.plot[to.plot$POP == "HG",]$POP <- "Pyrenees"
  
  GCAT_s <- c("ANDALUCIA","CATALUÑA","MADRID","EXTREMADURA","ARAGON","CYL","CLM","MURCIA","VALENCIANA","GALICIA","MADRID")
  to.plot[to.plot$POP %in% GCAT_s,]$SUPERPOP <- "GCAT"
  
  CAT_S <- c("València", "Catalonia", "Balearic_Is.")
  to.plot[to.plot$POP %in% CAT_S,]$SUPERPOP <- "People from Ibiza"
  return(to.plot)
} 

setwd("~/PROJECTS.JORGE/PCA.DEMOGRAPHY/General.Plots.Lazar/")

eigenvec <- read_table("GCAT.lazar.F.C.HGDP.euronly.evec", col_names = T)
eigenval <- scan("GCAT.lazar.F.C.HGDP.euronly.eval")
eigenval <- sapply(eigenval, function(x) x/length(eigenval))
labs <- read_table("../idfile.grouped.txt", col_names = F)
colnames(eigenvec)[1] <- "IID"
colnames(labs) <- c("IID","POP")

to.plot <- merge(eigenvec,labs,by="IID", sort = F)
to.plot$SUPERPOP <- "EUROPE"
colnames(to.plot)[c(2,3)] <- c("PC1","PC2")
F_s <- c("French_Basque","French","PAR","BR","IeV","PdD","NO","BdR","HG")
to.plot[to.plot$POP %in% F_s,]$POP <- "France"
GCAT_s <- c("ANDALUCIA","CATALUÑA","MADRID","EXTREMADURA","ARAGON","CYL","CLM","MURCIA","VALENCIANA","GALICIA","MADRID")
to.plot[to.plot$POP %in% GCAT_s,]$SUPERPOP <- "GCAT"
CAT_S <- c("València", "Catalonia", "Balearic_Is.")
to.plot[to.plot$POP %in% CAT_S,]$SUPERPOP <- "People from Ibiza"

change.names.aacc <- function(to.plot) {
  to.plot$POP <- gsub("ANDALUCIA","Andalusia", x= to.plot$POP)
  to.plot$POP <- gsub("Catalonia","Catalonia_B", x= to.plot$POP)
  to.plot$POP <- gsub("CATALUÑA","Catalonia_G", x= to.plot$POP)
  to.plot$POP <- gsub("ARAGON","Aragon", x= to.plot$POP)
  to.plot$POP <- gsub("MADRID","Madrid", x= to.plot$POP)
  to.plot$POP <- gsub("GALICIA","Galicia", x= to.plot$POP)
  to.plot$POP <- gsub("EXTREMADURA","Extremadura", x= to.plot$POP)
  to.plot$POP <- gsub("MURCIA","Murcia", x= to.plot$POP)
  to.plot$POP <- gsub("València","València_B", x= to.plot$POP)
  to.plot$POP <- gsub("VALENCIANA","València_G", x= to.plot$POP)
  return(to.plot)
}

to.plot <- change.names.aacc(to.plot)

ss <- subset(to.plot,to.plot$SUPERPOP == "GCAT" |  to.plot$SUPERPOP == "People from Ibiza")
to.plot$POP <- as.factor(to.plot$POP)
levels(as.factor(to.plot$POP))


shapes <- c(1,16,16,16,16,16,16,16,16,2,16,16,16,3,4,5,6,7,16,16)
c_palette <- c("grey","#00FF7B", "#FF3500", "#C100FF", "#2300FF", "brown", "#F6FF00",
               "#007BFF", "#F600FF", "grey", "#00FFE5", "black", "darkgreen","grey","grey","grey","grey","grey" , "#FF9E00", 
               "#8B864E", "#FF6A00", 
               "#FF006A","#FF0035","#0012FF","#FFD300","#58FF00","#5800FF","#8D00FF",
               "darkgrey","#C1FF00","#0046FF","#8DFF00","#FF009E","#00B0FF")

plot.pca(to.plot, ss, POP ,SUPERPOP) 

ggsave("PLOT.EUR.OCTUBRE.",
       plot = last_plot(),
       device = "png",
       scale = 1,
       width = 24,
       height = 24,
       units = "cm",
       dpi = 300)

##### CHANGES 16 Septiembre Figure s7
setwd("~/PROJECTS.JORGE/PCA.DEMOGRAPHY/")
eigenvec <- read_table("SPAIN.FRANCE.evec", col_names = T)
eigenval <- scan("SPAIN.FRANCE.eval")
eigenval <- sapply(eigenval, function(x) x/length(eigenval))
labs <- read_table("idfile.grouped.txt", col_names = F)

to.plot <- do.changes.to.data(eigenvec,eigenval,labs)
to.plot <- change.names.aacc(to.plot)

# MASK

shapes <- c(rep(1:25),15,16,17)
c_palette <- c("#00FF7B", "#FF3500", "#C100FF", "#2300FF", "#23FF00", "#F6FF00",
               "#007BFF", "#F600FF", "#00FFE5", "black", "#00FF46", "#FF9E00", 
               "#00FF12", "#FF6A00", "#FF0000",
               "#FF006A","#FF0035","#0012FF","#FFD300","#58FF00","#5800FF","#8D00FF",
               "darkgrey","#C1FF00","#0046FF","#8DFF00","#FF00D3","#FF009E","#00B0FF")



ss <- subset(to.plot,to.plot$SUPERPOP == "GCAT" |  to.plot$SUPERPOP == "France")
ss <- ss[!ss$POP == "French_Basque",]
ss <- ss[!ss$POP == "French",]
to.plot <- to.plot[!to.plot$POP == "French",]

p <- plot.pca(to.plot, ss, POP,SUPERPOP) 
p + geom_point(data = to.plot[to.plot$SUPERPOP == "People from Ibiza",], aes(x=PC1,y=PC2), colour= "grey")

ggsave("PLOT.IBE+FRA.SEPTIEMBRE.",
       plot = last_plot(),
       device = "png",
       scale = 1,
       width = 24,
       height = 24,
       units = "cm",
       dpi = 300)

####### Ancient DNA plot 18 sept

rm(eigenvec,eigenval,labs,ss,to.plot)
dev.off()

eigenvec <- read_table("aDNA.DEMO.GRCh38.EUR.evec", col_names = T)
eigenval <- scan("aDNA.DEMO.GRCh38.EUR.eval")
eigenval <- sapply(eigenval, function(x) x/length(eigenval))

labs <- read_table("idfile.grouped.txt", col_names = F)
alabs <- read_table("a.poplabs.txt", col_names = F)
labs <- rbind(labs,alabs)

to.plot <- do.changes.to.data(eigenvec, eigenval, labs  )
F_s <- c("Alsace","Pyrenees","Britanny","Nord","Dordogne", "Paris")
to.plot[to.plot$POP %in% F_s,]$POP <- "French"
to.plot$POP[to.plot$POP == "North_Italian"] <- "Italian"
to.plot$POP[to.plot$POP == "Tuscan"] <- "Italian"

plot.pca <- function (full.df, color.column, polygon.column) {
  p <- ggplot(full.df, aes(PC1, PC2, col={{color.column}}, shape = {{color.column}})) + geom_point(size = 1.5) + 
    coord_equal() + theme_light() +  scale_shape_manual(values= shapes) + 
    scale_color_discrete(type = c_palette) +
    xlab(paste0("PC1 (", signif(eigenval[1]*100, 3), "%)")) + 
    ylab(paste0("PC2 (", signif(eigenval[2]*100, 3), "%)")) 
  return(p)
}
to.plot <- change.names.aacc(to.plot)

to.plot$POP <- as.factor(to.plot$POP)
levels(to.plot$POP)

shapes <- c(1,16,16,16,16,
            16,16,16,16,2,
            3,4,16,5,6,
            16,16,7,8,9,
            11,12,16,16
            )
c_palette <- c("grey","#00FF7B", "#FF3500", "#C100FF", "#2300FF", 
               "brown", "#F6FF00","#007BFF", "#F600FF", "grey",
               "grey","grey", "#00FFE5","grey","grey", 
               "black", "darkgreen","grey","grey","grey",
               "grey","grey", "#FF9E00","#8B864E", "#FF6A00", 
               "#FF006A","#FF0035","#0012FF","#FFD300","#58FF00","#5800FF","#8D00FF",
               "darkgrey","#C1FF00","#0046FF","#8DFF00","#FF009E","#00B0FF")

c_palette <- c("grey","#00FF7B", "#FF3500", "#C100FF", "#2300FF",
               "darkgreen", "#F6FF00", "#007BFF", "#8B864E", "grey", 
               "grey", "grey", "#00FFE5", "grey","grey", 
               "black","brown","grey","grey","grey",
               "grey", "grey", "#FF9E00", "#00FF12", "#FF6A00", "#FF0000",
               "#FF006A","#FF0035","#0012FF","#FFD300","#58FF00","#5800FF","#8D00FF",
               "darkgrey","#C1FF00","#0046FF","#8DFF00","#FF00D3","#FF009E","#00B0FF")

plot.pca(to.plot, POP, POP)  

ggsave("PLOT.EUROPE.ANCIENT.OCTB",
       plot = last_plot(),
       device = "svg",
       scale = 1,
       width = 24,
       height = 24,
       units = "cm",
       dpi = 300)

