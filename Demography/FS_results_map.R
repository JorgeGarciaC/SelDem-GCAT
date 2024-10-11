##### See map of FS results

rm(list=ls())

library(tidyverse)
library(MoMAColors)
library(mapSpain)
library(sf)
library(grid)
library(scatterpie)
ibs <- read_table("~/PROJECTS.JORGE/CLUSTERS/labs.in.Iberian.Peninsula2.txt", col_names = F)

fulllabs <- read_table("~/PROJECTS.JORGE/CLUSTERS/NEW.full.labs")

colnames(ibs) <- c("fullorder","force.clust")
ibs <- ibs[!ibs$force.clust == "OUTLIER",]

new_df <- merge(ibs,fulllabs)

ONLY.GCAT <- new_df[grep("EGA",new_df$fullorder),]

cts <- ONLY.GCAT %>% count(fulllabs,force.clust)

cts$force.clust <- factor(cts$force.clust)

cts <- cts %>% mutate(proportion = n/sum(n)*100)

a <- esp_get_ccaa()
#a <- esp_get_ccaa(ccaa = c("Andalucía", "Aragón","La Rioja", "Cantabria","Madrid","Castilla y León","Castilla-La mancha","Asturias","País Vasco","Navarra", "Islas Baleares", "Galicia","Extremadura","Murcia","Comunidad Valenciana", "Cataluña"), moveCAN = F)

centroids <- st_centroid(a)
centroids <- st_transform(centroids, crs = st_crs(a)) # We ensure that we are using same coordinate system
centroids_df  <- data.frame(
  fulllabs = centroids$ccaa.shortname.es,
  latitude = st_coordinates(centroids)[,2],
  longitude = st_coordinates(centroids)[,1]
)

# Create centroids to plot

cts[cts$fulllabs == "ANDALUCIA",1] <- "Andalucía"
cts[cts$fulllabs == "ARAGON",1] <- "Aragón"
cts[cts$fulllabs == "CYL",1] <- "Castilla y León"
cts[cts$fulllabs == "CLM",1] <- "Castilla-La Mancha"
cts[cts$fulllabs == "CATALUÑA",1] <- "Cataluña"
cts[cts$fulllabs == "EXTREMADURA",1] <- "Extremadura"
cts[cts$fulllabs == "GALICIA",1] <- "Galicia"
cts[cts$fulllabs == "VALENCIANA",1] <- "Comunidad Valenciana"
cts[cts$fulllabs == "MADRID",1] <- "Madrid"
cts[cts$fulllabs == "MURCIA",1] <- "Murcia"

cts_coords <- merge(cts,centroids_df)

# Transform data into wider format

wider_cts <- cts_coords %>% select(fulllabs,force.clust,proportion) %>% pivot_wider(names_from = force.clust, values_from = proportion ) %>% replace_na(list(ARA = 0, VAL = 0, WEST = 0, CAT = 0))
wider_cts_part2 <- cts_coords %>% group_by(fulllabs) %>%  #summarise(n = sum(n), latitude, longitude) deprecated
  reframe(n = sum(n), longitude,latitude) %>% distinct()

cts_final <- merge(wider_cts,wider_cts_part2)

clnames <- c("fulllabs","ARA","CAT","WEST","MED","n","longitude","latitude") 
colnames(cts_final) <- clnames
cts_final$n <- as.numeric(cts_final$n)

can_box <- esp_get_can_box()
can_prov <- esp_get_can_provinces()

pal_color=c("#365C83", "#384351", "#4D8F8B", "#CDD6AD")

ggplot(a) +
  geom_sf() +
  #geom_sf(data = can_prov) +
  #geom_sf(data = can_box) + 
  theme_void() +
  geom_scatterpie(aes(x=longitude, y = latitude, group = fulllabs, r=log(n+1)/5 ), data=cts_final, cols = c("ARA","WEST","MED","CAT"), label_show_ratio = T,
                  legend_name = "Cluster", color=NA) + scale_colour_moma_d("Picabia") + scale_fill_manual(values = pal_color) +
  coord_sf(crs = st_crs(a))


ggsave(filename= "FS_pie.png",
       plot= last_plot(),
       device = "png",
       bg = "transparent",
       width = 15, 
       height = 15, 
       units = "cm", 
       dpi = 800)
