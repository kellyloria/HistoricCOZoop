## ---------------------------
## Exploring some community comp between 4 modern zoop. surverys in the the Southern Rockies 
##
## Author: Kelly A. Loria
## Date Created: 2020-04-15
## Email: kelly.loria@colorado.edu
##
## ---------------------------
## Load packages:
library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyverse)
library(zoo)
library(vegan)
library(ape)
library(devtools)


## ---------------------------
# File path setup:
#if (dir.exists('~/Users/kellyloria/Documents/Publications/Historical_Limno')){
#  inputDir<- '~/Users/kellyloria/Documents/Publications/Historical_Limno'
#  outputDir<- '~/Users/kellyloria/Documents/Publications/Historical_Limno/DataShare' 
#}

## ---------------------------
# I. For C7 output corrections 
# read in chl-a extractions: water_qualityCHLA.csv
d <- read.csv("/Users/kellyloria/Documents/Publications/Historical_Limno/PresAbs_RockMtsZoopTaxa.csv", header=T)
names(d)
summary(d)

# Clean up data to aggregate the RMNP15 and CL16 to get all observations of occurance in one row
d1 <- d %>% 
  group_by(Dataset, Lake) %>% 
  summarise(
            "Alonella"= sum(Alonella.sp., na.rm = T),
            "D.ambigua"= sum(Daphnia.ambigua, na.rm = T),
            "D.catawba"= sum(Daphnia.catawba, na.rm = T),
            "D.longispinus"= sum(Daphnia.longispinus, na.rm = T),
            "D.middendorfiana"= sum(Daphnia.middendorfiana, na.rm = T),
            "D.pulicaria"= sum(Daphnia.pulex.pulicaria, na.rm = T),
            "D.rosea"= sum(Daphnia.rosea, na.rm = T),
            "D.schodleri"= sum(Daphnia.schodleri, na.rm = T),
            "H.gibberum"= sum(Holopedium.gibberum, na.rm = T),
            "D.unidentifiable"= sum(Daphnia..unidentifiable., na.rm = T),
            "Diaphanosoma"= sum(Diaphanosoma, na.rm = T),
            "S.mucronata"= sum(Scapholeberis.mucronata, na.rm = T),
            "Macrothrix"= sum(Macrothrix.spp, na.rm = T),
            "H.shoshone"= sum(Diaptomus.shoshone., na.rm = T),
            "Epischura"= sum(Epischura.sp., na.rm = T),
            "D.thomasi"= sum(Diacyclops.thomasi, na.rm = T),
            "Macrocyclops"= sum(Macrocyclops.sp., na.rm = T),
            "Bosmina"= sum(Bosmina.longirostrus, na.rm = T),
            "Chydoridae"= sum(Chydorus.sp., na.rm = T),
            "Ascomorpha"= sum(Ascomorpha.sp., na.rm = T),
            "Asplanchna"= sum(Asplanchna.sp., na.rm = T),
            "Brachionus"= sum(Brachionus.sp., na.rm = T),
            "Collotheca"= sum(Collotheca.pelagica, na.rm = T),
            "Conochiloides"= sum(Conochiloides.Conochilus.spp., na.rm = T),
            "Euchlanis"= sum(Euchlanis.sp., na.rm = T),
            "Filinia"= sum(Filinia.terminalis, na.rm = T),
            "Gastropus"= sum(Gastropus.sp., na.rm = T),
            "Kellicottia"= sum(Kellicotia.longispinus, na.rm = T),
            "Keratella"= sum(Keratella.sp., na.rm = T),
            "Monostyla"= sum(Monostyla.sp...see.Lecane., na.rm = T),
            "NotholcaS"= sum(Notholca.squamula, na.rm = T),
            "NotholcaF"= sum(Notholca.folicea, na.rm = T),
            "Trichocerca"= sum(Trichocerca.sp., na.rm = T),
            "Tylotrocha"= sum(Tylotrocha.monopus, na.rm = T),
            "rot.unknown"= sum(rot.unknown, na.rm = T),
            "polyartha"= sum(polyartha.sp, na.rm = T),
            "Ostracoda"= sum(ostrocod, na.rm = T),
            "Chaoborus"= sum(Chaoborus..albicus..flavicans, na.rm = T), #
            "Chronomid"= sum(Midge, na.rm = T),
            "Ephemeroptera"= sum(Ephemeroptera, na.rm = T),
            "Odonata"= sum(Odonata, na.rm = T),
            "Plecoptera"= sum(Plecoptera, na.rm = T),
            "Mite"= sum(Mite, na.rm = T))

summary(d1)

# for richnness analysis convert all values <1 to 1
d2 <- as.data.frame(d1) %>%
  mutate_if(is.numeric, ~1 * (. != 0)) 
d3 <- dplyr::select(d2, 3:45)
names(d3)



# With this command, you`ll perform a NMDS and plot the results
d3 %>%
  metaMDS(trace = F) %>%
  ordiplot(type = "none") %>%
  text("sites")

# First step is to calculate a distance matrix. See PCOA for more information about the distance measures
# Here we use bray-curtis distance, which is recommended for abundance data
dist <- vegdist(d3,  method = "jaccard")

# In this part, we define a function NMDS.scree() that automatically 
# performs a NMDS for 1-10 dimensions and plots the nr of dimensions vs the stress
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# Use the function that we just defined to choose the optimal nr of dimensions
NMDS.scree(dist)


# Because the final result depends on the initial 
# random placement of the points 
# we`ll set a seed to make the results reproducible
set.seed(2)

# Here, we perform the final analysis and check the result
NMDS1 <- metaMDS(dist, k = 2, trymax = 100, trace = F)
# Do you know what the trymax = 100 and trace = F means?
# Let's check the results
NMDS1

# If you don`t provide a dissimilarity matrix, metaMDS automatically applies Bray-Curtis. So in our case, the results would have to be the same
NMDS2 <- metaMDS(d3, k = 2, trymax = 100, trace = F)
NMDS2

plot(NMDS1, type = "t")


NMDS3 <- metaMDS(d3, k = 2, trymax = 100, trace = F, autotransform = FALSE, distance="bray")
plot(NMDS3)
plot(NMDS3, display = "sites", type = "n")
points(NMDS3, display = "sites", col = "red", cex = 1.25)
text(NMDS3, display ="species")

# Alternatively, you can use the functions ordiplot and orditorp
ordiplot(NMDS3, type = "n")
orditorp(NMDS3, display = "species", col = "red", air = 0.01)
orditorp(NMDS3, display = "sites", cex = 1.1, air = 0.01)