# Attempt to run a PERMANOVA
library(vegan)
library(dplyr)

#read in data
crusade_df<-read.csv("CRUSADE_RA_T.dat")

#clean data?
#Need to have data completely cleaned of character data. Only have numerics.

crusade_filt <- crusade_df %>%
  mutate()%>%
  dplyr::filter(Timepoint == "T3") %>%
  dplyr::filter(DNA.Source == "Bacterioplankton") %>%
  # dplyr::filter(Water == "Unfiltered") %>%
  dplyr::filter(!Organism == "Start Water") %>%
  dplyr::select(-c(Label,Timepoint ,DNA.Source))

crusade_abun <- crusade_filt %>%
  dplyr::select(-c(Water,Organism))

adonis(crusade_abun ~ Organism, crusade_filt, perm=1000, method="bray", set.seed(100))


