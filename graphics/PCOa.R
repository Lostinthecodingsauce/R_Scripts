library(dplyr)
library(ggplot2)
library(readr)
library(MASS)
library(ape)
library(vegan)
library(factoextra)

#Read in data
Cr_Data<-read.csv("CRUSADE_RA_NT.dat")

#filter data for only the groups I want
microb_f_abun <-Cr_Data %>%
  mutate_if(is.numeric, sqrt)%>% #Normalizing data (arcsin(sqrt(x)))
  mutate_if(is.numeric, asin) %>%
  dplyr::filter(DNA.Source == "Bacterioplankton") %>%
  dplyr::filter(!Organism == "Water Control", !Organism == "State Water", !Organism == "Start Water")

#remove original dataframe
rm(Cr_Data)

#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- as.matrix(microb_f_abun %>%
                        dplyr::select(-c(Label:DNA.Source)),
                      dinames = list(paste("Label", 1:31, sep = ",")))

veg_bray <- vegdist(microb_f, "bray")

pc_scores<-pcoa(veg_bray)

pc_scores$Sample <- microb_f_abun$Label
pc_scores$Water_Treatment <- microb_f_abun$Water.Treatment
pc_scores$DNA_source <- microb_f_abun$DNA.Source
pc_scores$Organism <- microb_f_abun$Organism

barplot(pc_scores$values$Relative_eig[1:10])

# S3 method for pcoa
biplot.pcoa(pc_scores, pc_scores$Organism)

biplot(pc_scores, Y=NULL, col = microb_f_abun$Organism, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1)

pco_scores <- as.data.frame(pc_scores$vectors)
write.csv(pco_scores, "PCo_scores.dat")
write.csv(microb_f_abun, "PCo_reference.dat")
