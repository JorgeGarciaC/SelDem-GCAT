library(qvalue)
library(tidyverse)
library(biomaRt)
library(GenomicRanges)
library(ggrepel)
library(httr)
library(jsonlite)
library(xml2)

args = commandArgs(trailingOnly=TRUE)

############
# Annotation

new_df_in <- args[1]
vep_in <- args[2]
name_pop_1 <- args[3]
name_pop_2 <- args[4]
file_out <- args[5]


vep_xp <- read_tsv(vep_in, col_names = T)
new_df <- read_csv(new_df_in)

vep_xp <- vep_xp[c("#Uploaded_variation", "Allele", "Consequence", "IMPACT","SYMBOL", 
                   "BIOTYPE","CADD_PHRED","CADD_RAW","PHENOTYPES")]
names(vep_xp)[1] <- "refsnp_id" 

# Reformat VEP file
print("Reformat VEP file")
summ_vep <- vep_xp %>% group_by(refsnp_id) %>%
  summarise((across(everything(), ~(toString(na.omit(.))))))

for (j in 2:dim(summ_vep)[2]) { # Columna
  for (i in 1:dim(summ_vep)[1]){ # fila
    aux_str <- str_split(str_remove(summ_vep[i,j], "-"), ",") # por comodidad parto la celda en lista
    summ_vep[i,j] <- paste(unique(str_remove(aux_str[[1]], " ")), collapse = ",") 
    # extraigo la lista, quito los espacios y miro los valores únicos. Estos los pego y los reintroduzco en la posición
  }
}

server <- "https://grch37.rest.ensembl.org/"
endp <- "/variation/human/"
filt <- "?pops=1"

summ_vep$GBR_f <- NA
summ_vep$YRI_f <- NA
summ_vep$EUR_f <- NA

print("Adding population frequencies")

for (i in 1:dim(summ_vep)[1]) {
  rsid <- summ_vep$refsnp_id[i] 
  ext <- paste(endp,rsid,filt, sep = "")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json")) # Llamada
  warn_for_status(r)
  
  json_r <- fromJSON(toJSON(content(r)))$populations
  # Miramos si hay frecuencias para la pob deseada. 
  
  if ("1000GENOMES:phase_3:EUR" %in% json_r$population) {
    aux_json <- json_r[json_r$population == "1000GENOMES:phase_3:EUR",]
    # extraemos las f y las juntamos todas para meterlas en una casilla
    summ_vep$EUR_f[i] <- paste(paste(aux_json$allele, aux_json$frequency, sep = ":"), collapse = ";")
  }  
  
  if ("1000GENOMES:phase_3:GBR" %in% json_r$population) {
    aux_json <- json_r[json_r$population == "1000GENOMES:phase_3:GBR",]
    # extraemos las f y las juntamos todas para meterlas en una casilla
    summ_vep$GBR_f[i] <- paste(paste(aux_json$allele, aux_json$frequency, sep = ":"), collapse = ";")
    
  }
  if ("1000GENOMES:phase_3:YRI" %in% json_r$population) {
    aux_json <- json_r[json_r$population == "1000GENOMES:phase_3:YRI",]
    # extraemos las f y las juntamos todas para meterlas en una casilla
    summ_vep$YRI_f[i] <- paste(paste(aux_json$allele, aux_json$frequency, sep = ":"), collapse = ";")
    
  }
}

annotation <- new_df %>% 
  dplyr::select(id,chr,normxpehh,n_peak,p1,p2)

f_xp_anno <- merge(annotation, summ_vep, by.x = "id", by.y = "refsnp_id")

names(f_xp_anno)[5] <- paste("freq_", name_pop_1)
names(f_xp_anno)[6] <- paste("freq_", name_pop_2)

write_csv(f_xp_anno, paste(file_out,".csv", sep = ''))


if (args[6] == "PLOT_TRUE")
{
  print("Plotting")
  
  ########
  # Manhattan plotting
  
  don <- new_df[, c("chr","pos", "normxpehh", "peak", "n_peak")] %>% 
    
    # Compute chromosome size
    group_by(chr) %>% 
    summarise(chr_len=max(pos)) %>% 
    
    # Calculate cumulative BPition of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(new_df, ., by=c("chr"="chr")) %>%
    
    # Add a cumulative BPition of each SNP
    arrange(chr, pos) %>%
    mutate( BPcum=pos+tot)
  
  
  # Then we need to prepare the X axis. 
  # Indeed we do not want to display the cumulative BPition of SNP in bp, 
  # but just show the chromosome name instead
  
  axisdf = don %>% 
    group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  
  # plot 
  p <- ggplot(don, aes(x=BPcum, y=normxpehh)) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    geom_hline(yintercept = 0, color = "black", size =1 ) +
    annotate(geom = "text", x = -0.5, y = c(0+1, 0-1), 
             label = c(name_pop_1, name_pop_2), size = 6, angle =90) +
    #ggtitle(label = "Gcat Vs Gbr") +
    geom_point(data=subset(don, peak==TRUE), color="lightsalmon", size=1)  
    #geom_text_repel(aes(label=ifelse(peak == T,as.character(n_peak)
     #                                ,'')),hjust=0,vjust=0, size=3, max.overlaps = 50) # cambiar p por IHS
  
  png(paste0(file_out, "_manhattan.png"), width = 1280, height = 720)
  print(p)
  dev.off()
  
}