## Written by Zach Quinlan 12/17/2018
## This should be able to take the raw outputs from mothur and give a working df

# Libraries for cleaning
library(DescTools)
library(tidyverse)
library(data.table)

# Libraries for NMDS, PERMANOVA, PCoA
library(vegan)
library(MASS)
library(factoextra)
library(ape)


# Read in data frames -----------------------------------------------------
##  You need both taxon and RA outputs from mothur. They should be in csv format. Just need to change the ending to .dat
## This part is untested. SHould work...

taxon_file <- read_csv("your_taxon_file.dat")
ra_file <- read_csv("your_relative_abundance_table.dat")

# Join two tables and seperate taxonomy -----------------------------------
left_join(taxon_file, ra_file, by = NULL)
separate(Taxonomy, paste(taxonomy, 1:7, sep = "."), sep = ";")


# Cleaning all naming columns ---------------------------------------------
## At this point I can begin using custom functions
cleaned <- crusade_df%>%
  dplyr::mutate(Order = case_when(Order %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Order)))%>%
  dplyr::mutate(Family = case_when(Order == "unclassified" ~ "",
                                   TRUE ~ as.character(Family)))%>%
  dplyr::mutate(Genus = case_when(Order == "unclassified" ~ "",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Order == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(Family = case_when(Family %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                   TRUE ~ as.character(Family)))%>%
  dplyr::mutate(Genus = case_when(Family == "unclassified" ~ "",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Family == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(Genus = case_when(Genus %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Genus == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(OTU = case_when(OTU %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "sp",
                                TRUE ~ as.character(OTU)))%>%
  group_by(Order, Family, Genus, OTU) %>%
  summarize_if(is.numeric, sum)%>%
  ungroup(data_sum_OTU)%>%
  unite(OFGO, c(Order, Family, Genus, OTU), sep = ";", remove = TRUE)

transformed <- cleaned%>%  
  mutate_if(is.numeric, sqrt)%>%#Normalizing data (arcsin(sqrt(x)))
  mutate_if(is.numeric, asin)

transposed <-transformed%>%
  gather(sample_code, RA, 2:32)%>%
  spread(OFGO, RA)


# Designing working df for stats and such ---------------------------------
working_df <- transposed%>%
  separate(sample_code, paste("C", 1:4, sep = "."), sep = "_", remove = FALSE)%>%
  rename(Organism = `C.1`,
         Water = `C.2`,
         Timepoint = `C.3`,
         DNA.Source = `C.4`)%>%
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
                                TRUE ~ "Bacterioplankton"))%>%
  unite(Name, c(Organism, Water, Timepoint, DNA.Source), sep = "      ", remove = FALSE)


# Subsetted dataframes -------------------------------------------------------

# Bacterioplankton whole df
bacterioplankton_df <- working_df %>%
  dplyr::filter(Timepoint == "T3") %>%
  dplyr::filter(DNA.Source == "Bacterioplankton") %>%
  dplyr::filter(!Organism == "Start Water") %>%
  dplyr::select(-c(Name, Timepoint, DNA.Source))


# Bacterioplanknton unfiltered
bacterio_unfilt <- bacterioplankton_df %>%
  dplyr::filter(Water == "Unfiltered")

bacterio_un_abun <- bacterio_unfilt %>%
  dplyr::select(-c(Water, Organism, sample_code))


# Bacterioplankton Filtered
bacterio_filt <- bacterioplankton_df %>%
  dplyr::filter(!Water == "Unfiltered")

bacterio_f_abun <- bacterio_filt %>%
  dplyr::select(-c(Water,Organism, sample_code))


# CCA Microbiome
cca_microb <- working_df %>%
  dplyr::filter(Timepoint == "T3") %>%
  dplyr::filter(!DNA.Source == "Bacterioplankton") %>%
  dplyr::filter(!Organism == "Start Water") %>%
  dplyr::select(-c(Name, Timepoint, DNA.Source, Water))

cca_microb_abun <- cca_microb %>%
  dplyr::select(-c(Organism, sample_code))



# Running PERMANOVAs ------------------------------------------------------
# Bacterioplankton unfiltered
adonis(bacterio_un_abun ~ Organism, bacterio_unfilt, perm=1000, method="bray", set.seed(100))

# Bacterioplankton filtered
adonis(bacterio_f_abun ~ Organism, bacterio_filt, perm=1000, method="bray", set.seed(100))

# CCA Microbiome
adonis(cca_microb_abun ~ Organism, cca_microb, perm=1000, method="bray", set.seed(100))

# Making PCoA's -----------------------------------------------------------

#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- as.matrix(cca_microb %>%
                        dplyr::select(-c(Organism, sample_code)),
                      dinames = list(paste("", 1:31, sep = ",")))

veg_bray <- vegdist(microb_f, "bray") #Bray-curtis distances

pc_scores<-pcoa(veg_bray) #Calculating scores

#This plots Eigenvalues
#They will allow you to choose the best axes to show how your data varies
ylimit = c(0, 1.1*max(pc_scores$values$Relative_eig))

Eigan <- barplot(pc_scores$values$Relative_eig[1:10], ylim= ylimit)
# Add values to the bars
text(x = Eigan, y = pc_scores$values$Relative_eig[1:10], label = pc_scores$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

# S3 method for pcoa
biplot(pc_scores, Y=NULL, col = microb_f_abun$Organism, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1)

pco_scores <- as.data.frame(pc_scores$vectors)
pco_scores$Sample.Code <- cca_microb$sample_code     # This will add reference labels to the PCoA scores

write.csv(pco_scores, "PCo_scores.dat")           # This will print the Scores to remake the figure in JMP because JMP is bae