# Shannon and Pielous
# These indecies work best on raw count data not RA data and make sure not to transform the data
# When it is standardized you insert bias into the dataset
library(vegan)
library(dplyr)

# read in data
microb_df <- read.csv("CR_for_Richness:Eveness.dat")


# Calculating indecies ----------------------------------------------------

# edit column names to contain no spaces
microb_abun <- microb_df %>% 
  dplyr::filter(!Label.1 == "CH",
                !Label.1 == "ST")
  
Rich <- microb_abun %>%
  dplyr::select(-c(1:2))

##  Main datasheet
indecies <- microb_abun%>%
  dplyr::select(c(1))

# Richness index (Shannan-Weaver Index)
Crusade_Shannon <-diversity(Rich)

# Eveness index (Pielous Eveness Index)
Crusade_Pielous <-Crusade_Shannon/log(specnumber(Rich))

indecies$Shannon <- Crusade_Shannon
indecies$Pielous <- Crusade_Pielous

# Designing report with all indecies and correct column names
indecies_report <-
  tidyr::separate(indecies, Label,sep = "_", into=(c("Organism", "Water", "Timepoint","DNA.Source")))%>%
  mutate(Organism = case_when(Organism == "CC" ~ "Calcium Carbonate Control",
                              Organism == "HR" ~ "Hydrolithon reinboldii",
                              TRUE ~ "Porolithon onkodes"))%>%
  mutate(Water = case_when(Water == "1F" ~ "Filtered",
                           Water == "2F" ~ "Filtered",
                           Water == "3F" ~ "Filtered",
                           TRUE ~ "Unfiltered"))%>%
  mutate(DNA.Source = case_when(DNA.Source == "CCA" ~ "CCA Microbiome",
                                TRUE ~ "Bacterioplankton"))

##  Auxillary spreadsheets are commented out. Main one with all indecies and names is still included
# write.csv(Crusade_Shannon,"Crusade_Shannon.dat")
# write.csv(Crusade_Pielous, "Crusade_Pielous.dat")
# write.csv(microb_abun, "Crusade_abun.dat")
write.csv(indecies_report, "Crusade_indecies.dat")


# ANOVAs ------------------------------------------------------------------
indecies_filt <- indecies_report%>%
  dplyr::filter(Water == "Filtered")

indecies_filt_HR <- indecies_filt%>%
  dplyr::filter(Organism == "Hydrolithon reinboldii")

indecies_filt_PO <- indecies_filt%>%
  dplyr::filter(!Organism == "Hydrolithon reinboldii")

summary(aov(Shannon ~ DNA.Source, data = indecies_filt))

t.test(Shannon ~ DNA.Source, data = indecies_filt)

t.test(Pielous ~ DNA.Source, data = indecies_filt)
