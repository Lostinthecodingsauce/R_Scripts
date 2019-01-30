#Data table needs to have microbial OTUS/class/etc. as response variables and Organism or grouping variable in one column
#Written for CRUSADE--Zaq 10/5/2018
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(janitor)
library(ggplot2)

#read in Data
#cr_raw_pie <- read.csv("Crusade_pie.dat")
cr_raw_pie <- class_filter%>% #This stupid line comes directly form another script... Im dumb
  gather(Code, RA, 3:33)%>%
  separate(Code, paste("C", 1:4, sep = "."), sep = "_")%>%
  rename(Organism = `C.1`,
         Water = `C.2`,
         Timepoint = `C.3`,
         DNA.Source = `C.4`,
         Class = `Taxonomy.3`,
         Order = `Taxonomy.4`)%>%
  mutate(Organism = case_when(Organism == "CC" ~ "Calcium Carbonate Control",
                              Organism == "CH" ~ "Water Control",
                              Organism == "HR" ~ "Hydrolithon reinboldii",
                              Organism == "PO" ~ "Porolithon onkodes",
                              TRUE ~ "Start Water"))%>%
  mutate(Water = case_when(Water == "1F" ~ "Filtered",
                             Water == "2F" ~ "Filtered",
                             Water == "3F" ~ "Filtered",
                             TRUE ~ "Unfiltered"))%>%
  mutate(DNA.Source = case_when(DNA.Source == "CCA" ~ "CCA Microbiome",
                                  TRUE ~ "Bacterioplankton"))

#Filter by DNA source and group by organism
#Only from CCA Microbiome. No Water control or Start Water
cr_CCA_combine <-cr_raw_pie%>%
  dplyr::filter(!DNA.Source == "Bacterioplankton")%>%
  dplyr::filter(!Organism == "Start Water",
                !Organism == "Water Control")%>%
  group_by(Organism, Order, Class)%>%
    summarize_if(is.numeric, mean)%>%
  ungroup(cr_CCA_combine)

#Only from Bacterioplankton No Water control or Start Water
cr_bact_filt_combine <-cr_raw_pie%>%
  dplyr::filter(DNA.Source == "Bacterioplankton")%>%
  dplyr::filter(Water == "Filtered")%>%
  dplyr::filter(!Organism == "Start Water",
                !Organism == "Water Control")%>%
  group_by(Organism, Order, Class)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup(cr_bact_filt_combine)

cr_bact_unf_combine <-cr_raw_pie%>%
    dplyr::filter(DNA.Source == "Bacterioplankton")%>%
    dplyr::filter(!Water == "Filtered")%>%
    dplyr::filter(!Organism == "Start Water",
                  !Organism == "Water Control")%>%
  group_by(Organism, Order, Class)%>%
  summarize_if(is.numeric, mean)%>%
    ungroup(cr_bact_unf_combine)

#Filtering to only X Rabundances, tried to limit to ~10 different classes
# This can be used to find the 
# bact_filt_04 <- cr_bact_filt_combine%>%
#   gather(Class, RA, 2:65)%>%
#   group_by(Organism, Class)%>%
#   dplyr::filter(RA, sum(RA[RA>0.01]))%>%
#   ungroup(bact_filt_04)
# 
# bact_unf_04 <- cr_bact_unf_combine%>%
#   gather(Class, RA, 2:65)%>%
#   group_by(Organism, Class)%>%
#   dplyr::filter(RA, sum(RA[RA>0.01]))%>%
#   ungroup(bact_unf_04)
# 
# CCA_04 <- cr_CCA_combine%>%
#   gather(Class, RA, 2:65)%>%
#   group_by(Organism, Class)%>%
#   dplyr::filter(RA, sum(RA[RA>0.02]))%>%
#   ungroup(CCA_04)
# 

#Calculating 1-sum(RA) of bacts that did not make cutoff
# other_unf <- bact_unf_04%>%
#   group_by(Organism)%>%
#   mutate(RA, 1-sum(RA))
# 
# other_filt <- bact_filt_04%>%
#   group_by(Organism)%>%
#   mutate(RA, 1-sum(RA))
# 
# other_CCA <- CCA_04%>%
#   group_by(Organism)%>%
#   mutate(RA, 1-sum(RA))

# View(other_unf)
# View(other_filt)
# View(other_CCA)

#Write data table based off of pie chart setup
#To do this in jmp export. But I think R makes better Pie Charts
# write.csv(other_filt, "Crusade_pie_FBact.dat")
# write.csv(other_unf, "Crusade_pie_UBact.dat")
# write.csv(other_CCA, "Crusade_pie_CCA.dat")

#Filter each Treatment by organism
PO_filt <- cr_bact_filt_combine%>%
  dplyr::filter(Organism == "Porolithon onkodes")

HR_filt <- cr_bact_filt_combine%>%
  dplyr::filter(Organism == "Hydrolithon reinboldii")

PO_unf <- cr_bact_unf_combine%>%
  dplyr::filter(Organism == "Porolithon onkodes")

HR_unf <- cr_bact_unf_combine%>%
  dplyr::filter(Organism == "Hydrolithon reinboldii")

PO_CCA <- cr_CCA_combine%>%
  dplyr::filter(Organism == "Porolithon onkodes")

HR_CCA <- cr_CCA_combine%>%
  dplyr::filter(Organism == "Hydrolithon reinboldii")


#Visualize Pie Charts
#Filtered treatments
PO_filt_ggp <-ggplot(PO_filt, aes(x= "", y=RA, fill = Class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "Porolithon onkodes filtered") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

HR_filt_ggp <-ggplot(HR_filt, aes(x= "", y=RA, fill = Class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "Hydrolithon reinboldii filtered") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

PO_filt_pie <- PO_filt_ggp + coord_polar("y", start = 0)
HR_filt_pie <- HR_filt_ggp + coord_polar("y", start = 0)

#Unfiltered Treatments
PO_unf_ggp <-ggplot(PO_unf, aes(x= "", y=RA, fill = Class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "Porolithon onkodes unf") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

HR_unf_ggp <-ggplot(HR_unf, aes(x= "", y=RA, fill = Class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "HR unf") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

PO_unf_pie <- PO_unf_ggp + coord_polar("y", start = 0)
HR_unf_pie <- HR_unf_ggp + coord_polar("y", start = 0)

#CCA Treatments
PO_CCA_ggp <-ggplot(PO_CCA, aes(x= "", y=RA, fill = Class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "Porolithon onkodes CCA") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

HR_CCA_ggp <-ggplot(HR_CCA, aes(x= "", y=RA, fill = Class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "HR CCA") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

PO_CCA_pie <- PO_CCA_ggp + coord_polar("y", start = 0)
HR_CCA_pie <- HR_CCA_ggp + coord_polar("y", start = 0)

#Saving Pie Charts as pdfs to wd
#Filtered
# ggsave(filename = "PO_filt.pdf", plot= PO_filt_pie)
# ggsave(filename = "HR_filt.pdf", plot= HR_filt_pie)
# 
# #Unfiltered
# ggsave(filename = "PO_unfilt.pdf", plot= PO_unf_pie)
# ggsave(filename = "HR_unfilt.pdf", plot= HR_unf_pie)
# 
# #CCA
# ggsave(filename = "PO_CCAt.pdf", plot= PO_CCA_pie)
# ggsave(filename = "HR_CCA.pdf", plot= HR_CCA_pie)
# 
# 
