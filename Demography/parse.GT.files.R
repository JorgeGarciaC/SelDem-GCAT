##### CI for fastGL ####
rm(list=ls())
library(tidyverse)
library(cowplot)
library(ggridges)

#Function to automatically process the output of fastGT
process_main <- function(main, ft = T) {
  str.m <- strsplit(main, " ")
  if (ft == T) {
    best.fit <- "one-date)"
  } else {best.fit <- str.m[[grep("conclusion",str.m)]][8]}
  
###### Multiple dates #####
  if (best.fit == "multiple-dates)") {
  
    fit_date <- 2
    
    ##### first date
    D1.S1 <- as.data.frame(str.m[24:25])
    D1.P.S1 <- as.numeric(D1.S1[1,2])
    D1.S1 <- D1.S1[-1,] 
    colnames(D1.S1) <- c("POP","P1")
    D1.S1$FIT <- 1
    D1.S1$S <- 1
    D1.S1$P1 <- as.numeric(D1.S1$P1)
    
    D1.S2 <- as.data.frame(str.m[26:27])
    D1.P.S2 <- D1.S2[1,2]
    D1.S2 <- D1.S2[-1,] 
    colnames(D1.S2) <- c("POP","P1")
    D1.S2$FIT <- 1
    D1.S2$S <- 2
    D1.S2$P1 <- as.numeric(D1.S2$P1)
    
    ##### Second date
    D2.S1 <- as.data.frame(str.m[30:31])
    D2.P.S1 <- D2.S1[1,2]
    D2.S1 <- D2.S1[-1,] 
    colnames(D2.S1) <- c("POP","P1")
    D2.S1$FIT <- 2
    D2.S1$S <- 1
    D2.S1$P1 <- as.numeric(D2.S1$P1)
    
    D2.S2 <- as.data.frame(str.m[32:33])
    D2.P.S2 <- D2.S2[1,2]
    D2.S2 <- D2.S2[-1,] 
    colnames(D2.S2) <- c("POP","P1")
    D2.S2$FIT <- 2
    D2.S2$S <- 2
    D2.S2$P1 <- as.numeric(D2.S2$P1)
    
    
    D1.S1$Pt <- D1.S1$P1 * as.numeric(D1.P.S1)
    D1.S2$Pt <- D1.S2$P1 * as.numeric(D1.P.S2)
    D2.S1$Pt <- D2.S1$P1 * as.numeric(D2.P.S1)
    D2.S2$Pt <- D2.S2$P1 * as.numeric(D2.P.S2)
    
    new_df <- rbind(D1.S1,D1.S2,D2.S1,D2.S2)
    new_df$name <- paste0(new_df$POP,"_", new_df$S)
    return(new_df)
    } 

##### One date ####
  if (best.fit == "one-date)") {
    D1.S1 <- as.data.frame(str.m[12:13])
    D1.P.S1 <- as.numeric(D1.S1[1,2])
    D1.S1 <- D1.S1[-1,] 
    colnames(D1.S1) <- c("POP","P1")
    D1.S1$FIT <- 1
    D1.S1$S <- 1
    D1.S1$P1 <- as.numeric(D1.S1$P1)
    
    D1.S2 <- as.data.frame(str.m[14:15])
    D1.P.S2 <- D1.S2[1,2]
    D1.S2 <- D1.S2[-1,] 
    colnames(D1.S2) <- c("POP","P1")
    D1.S2$FIT <- 1
    D1.S2$S <- 2
    D1.S2$P1 <- as.numeric(D1.S2$P1)
    
    D1.S1$Pt <- D1.S1$P1 * as.numeric(D1.P.S1)
    D1.S2$Pt <- D1.S2$P1 * as.numeric(D1.P.S2)
    
    new_df <- rbind(D1.S1,D1.S2)
    new_df$name <- paste0(new_df$POP,"_", new_df$S)
    return(new_df)
  } 
  
  if (best.fit == "one-date-multiway)") {
    D1.S1 <- as.data.frame(str.m[12:13])
    D1.P.S1 <- as.numeric(D1.S1[1,2])
    D1.S1 <- D1.S1[-1,] 
    colnames(D1.S1) <- c("POP","P1")
    D1.S1$FIT <- 1
    D1.S1$S <- 1
    D1.S1$P1 <- as.numeric(D1.S1$P1)
    
    D1.S2 <- as.data.frame(str.m[14:15])
    D1.P.S2 <- D1.S2[1,2]
    D1.S2 <- D1.S2[-1,] 
    colnames(D1.S2) <- c("POP","P1")
    D1.S2$FIT <- 1
    D1.S2$S <- 2
    D1.S2$P1 <- as.numeric(D1.S2$P1)
    
    D1.S1$Pt <- D1.S1$P1 * as.numeric(D1.P.S1)
    D1.S2$Pt <- D1.S2$P1 * as.numeric(D1.P.S2)
    
    new_df <- rbind(D1.S1,D1.S2)
    new_df$name <- paste0(new_df$POP,"_", new_df$S)
    
    # Sacamos la minor contribution
    Minor.D1.S1 <- as.data.frame(str.m[18:19])
    Minor.D1.P.S1 <- as.numeric(Minor.D1.S1[1,2])
    Minor.D1.S1 <- Minor.D1.S1[-1,] 
    colnames(Minor.D1.S1) <- c("POP","P1")
    Minor.D1.S1$FIT <- 1
    Minor.D1.S1$S <- 1
    Minor.D1.S1$P1 <- as.numeric(Minor.D1.S1$P1)
    
    Minor.D1.S2 <- as.data.frame(str.m[20:21])
    Minor.D1.P.S2 <- Minor.D1.S2[1,2]
    Minor.D1.S2 <- Minor.D1.S2[-1,] 
    colnames(Minor.D1.S2) <- c("POP","P1")
    Minor.D1.S2$FIT <- 1
    Minor.D1.S2$S <- 2
    Minor.D1.S2$P1 <- as.numeric(Minor.D1.S2$P1)
    
    Minor.D1.S1$Pt <- Minor.D1.S1$P1 * as.numeric(Minor.D1.P.S1)
    Minor.D1.S2$Pt <- Minor.D1.S2$P1 * as.numeric(Minor.D1.P.S2)
    
    minor <- rbind(Minor.D1.S1,Minor.D1.S2)
    minor$name <- paste0(minor$POP,"_", minor$S)
    
    return(list(minor,new_df))
  }
  
}

main <- readLines(con = "~/PROJECTS.JORGE/Pasar al portatil-20240307T164307Z-001/Pasar al portatil/FINALL.GT/globetrotter.1.EASTE.main.txt")
main.west <- readLines(con = "~/PROJECTS.JORGE/Pasar al portatil-20240307T164307Z-001/Pasar al portatil/FINALL.GT/globetrotter.1.WESTE.main.txt")

main <- readLines(con = "~/PROJECTS/GCAT.Demography.Actual.Project/CHROMO/01-20-2023/fGT/12-19-23/globetrotter.1.EASTE.main.txt")
main.west <- readLines(con = "~/PROJECTS/GCAT.Demography.Actual.Project/CHROMO/01-20-2023/fGT/12-20-23/globetrotter.1.WESTE.main.txt")

eas.esp <- process_main(main)
west.esp <- process_main(main.west)

eas.esp <- eas.esp[,c("POP","P1","S")]
west.esp <- west.esp[,c("POP","P1","S")]

eas.esp$TARGET <- "EAST"
west.esp$TARGET <- "WEST"

#### Plot sources proportions
to.plot <- rbind(eas.esp,west.esp)

p <- ggplot(to.plot[to.plot$S == 1,], aes(x = POP, y = TARGET, fill = P1)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(name="% Contribution", 
                       colours = c("lightblue","blue","#132B43")) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

b <- ggplot(to.plot[to.plot$S == 2,], aes(x = POP, y = TARGET, fill = P1)) +
  geom_tile(color = "black") +
  scale_fill_gradientn(name="% Contribution", 
                       colours = c("lightblue","blue","#132B43")) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

j <- plot_grid(b,p,
          labels = c("Mayor Source", "Minor Source"),
          ncol = 1, nrow = 2,label_y = 1.01)
j

##### Save plots
setwd("~/PROJECTS/GCAT.Demography.Actual.Project/CHROMO/01-12-2023/fGT/12-18-23/")
ggsave("Sources.guess.fGT.ALL.pdf",
       plot = j,
       device = "pdf",
       scale = 1,
       width = 18,
       height = 24,
       units = "cm",
       dpi = 300)

ggsave("Sources.guess.fGT.ALL.jpeg",
       plot = j,
       device = "jpeg",
       scale = 1,
       width = 18,
       height = 24,
       units = "cm",
       dpi = 300)

###### 
# bootstrap stimates
#####

rm(list=ls())
setwd("~/PROJECTS.JORGE/Pasar al portatil-20240307T164307Z-001/Pasar al portatil/FINALL.GT/")
setwd("/home/bioevo/PROJECTS/GCAT.Demography.Actual.Project/CHROMO/01-12-2023/fGT/12-18-23/")

boost.easte.ALL  <- read_table("globetrotter.1.EASTE.boot.txt")
boost.weste.ALL  <- read_table("globetrotter.1.WESTE.boot.txt")

boost.easte.ALL$experiment <- "EAST"
boost.weste.ALL$experiment <- "WEST"

to.plot <- rbind(boost.easte.ALL,boost.weste.ALL)
to.plot$date <- as.integer(1960-(to.plot$date1.est.boot * 29))

qtls <- to.plot %>% group_by(experiment) %>% select(date) %>% summarise(quantile(date, probs = c(0.025,0.5,0.975))) 
qtls

ggplot(to.plot, aes(x = date, y = experiment, fill= experiment)) +
  geom_density_ridges(scale=0.95, 
                      quantile_lines=TRUE,
                      quantile_fun=function(date,...)quantile(date, probs = c(0.025,0.5,0.975))) +
  theme_ridges() +
  labs() + 
  theme(legend.position = "none") +
  ggtitle( "Bootstrap on date inference of the admixture event") 

ggsave("Boostrap.fGT.ALL.pdf",
       plot = last_plot(),
       device = "pdf",
       scale = 1,
       width = 18,
       height = 24,
       units = "cm",
       dpi = 300)
