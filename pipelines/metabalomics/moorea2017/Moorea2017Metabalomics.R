rm(list = ls())


# Libraries needed --------------------------------------------------------

#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)

#PCoA, PERMANOVA
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)

# Loading in Data and Joining Tables ---------------------------------------------------------


analog_hits <- read_tsv("Analog_hits.tsv")%>%
  rename(Scan_Number = 'Scan_analog')
true_hits <-read_tsv("library_hits_true.tsv")%>%
  rename(Scan_Number = 'Scan')
old_hits <- read_csv("Moorea_SPEDOM_MSMS_Polished_April2018.csv")%>%
  rename(Scan_Number = 'row ID')
Node_info <- read_tsv("Node_Info.clustersummary.tsv")

## Joining tables

new_run <- full_join(full_join(Node_info, true_hits, "Scan_Number"),analog_hits, "Scan_Number")

both_runs <- full_join(new_run, old_hits, "Scan_Number")

both_runs_by_old <- right_join(new_run, old_hits, "Scan_Number")

scan_IDs <- as.data.frame(old_hits$Scan_Number)

# Transposing and selecting only sample data
samples <- both_runs_by_old%>%
  select(c('Feature_Name', '% of Total(D_WA_1_T0D)':'% of Total(R_TR_3_TFN)'))%>%
  gather(sample_ID, RA, 2:278)%>%
  spread(Feature_Name, RA)

# Some Samples had - instead of _
samples$sample_ID <- gsub("-", "_", samples$sample_ID)
samples$sample_ID <- gsub("% of Total\\(", "", samples$sample_ID)
samples$sample_ID <- gsub(")", "", samples$sample_ID)

# renaming columns and row metadata
moorea_wdf <-samples%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint", "DOM_source"), sep = "_")%>%
  filter(!Experiment %like% "Blank")%>%
  filter(!Experiment == "SE",
         !Experiment == "SR")%>%
  dplyr::mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       Experiment == "M" ~ "mordor",
                                       Experiment == "R" ~ "RR3",
                                       TRUE ~ as.character(Experiment)))%>%
  dplyr::mutate(Organism = case_when(Organism == "AH" ~ "Amansia 13c",
                                     Organism == "AL" ~ "Algae_labile",
                                     Organism == "AQ" ~ "Aquil ASW",
                                     Organism == "AR" ~ "Algae_recalcitrant",
                                     Organism == "AW" ~ "Algae_control",
                                     Organism == "CC" ~ "CCA",
                                     Organism == "CL" ~ "Coral_labile",
                                     Organism == "CR" ~ "Coral_recalcitrant",
                                     Organism == "CW" ~ "Coral_control",
                                     Organism == "DH" ~ "Dictyota 13c",
                                     Organism == "DL" ~ "Dictyota 12c",
                                     Organism == "DT" ~ "Dictyota",
                                     Organism == "OW" ~ "Outside control",
                                     Organism == "PL" ~ "Porites lobata",
                                     Organism == "PV" ~ "Pocillopora verrucosa",
                                     Organism == "TR" ~ "Turf",
                                     Organism == "TW" ~ "Tent control",
                                     Organism == "WA" ~ "Water control",
                                     TRUE ~ as.character(Organism)))

moorea_wdf <- apply(moorea_wdf[6:10031], 2, function(x) gsub("%", "", x))
# moorea_wdf$Organismal_clade <- moorea_wdf$Organism
# 
# moorea_wdf <- moorea_wdf%>%
#   dplyr::mutate(Organismal_clade = case_when(Organismal_clade == "Outside control" ~ "water",
#                                              Organismal_clade == "Tent control" ~ "water",
#                                              Organismal_clade == "Water control" ~ "water",
#                                              Organismal_clade == "Porites lobata" ~ "coral",
#                                              Organismal_clade == "Pocillopora verrucosa" ~ "coral",
#                                              Organismal_clade %like any% c(1:20) ~ "environmental",
#                                              Organismal_clade == "Dictyota" ~ "algae",
#                                              Organismal_clade == "Turf" ~ "algae",
#                                              Organismal_clade == "Coral_control" ~ "coral",
#                                              Organismal_clade == "Coral_labile" ~ "coral",
#                                              Organismal_clade == "Coral_recalcitrant" ~ "coral",
#                                              Organismal_clade == "Algae_control" ~ "algae + coral",
#                                              Organismal_clade == "Algae_labile" ~ "algae + coral",
#                                              Organismal_clade == "Algae_recalcitrant" ~ "algae + coral",
#                                              Organismal_clade == "Algae" ~ "algae",
#                                              TRUE ~ as.character(Organismal_clade)))%>%
#   dplyr::filter(!Organismal_clade == "Blank")
#   
# moorea_wdf["Timepoint"][is.na(moorea_wdf["Timepoint"])] <- 0

dorcierr <- moorea_wdf%>%
  filter(Experiment == "dorcierr")%>%
  select(-"DOM_source")

mordor <- moorea_wdf%>%
  filter(Experiment == "mordor")

big3 <- moorea_wdf%>%
  filter(!Experiment == "LoRDI",
         !Experiment == "SPIFFy")

# dorcierr further cleaning -----------------------------------------------
dorcierr_wdf <- dorcierr%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = -1)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))

dorcierr_final <- dorcierr_wdf%>%
  filter(Timepoint == "TF")

dorcierr_transformations <- dorcierr_wdf%>%
  gather()

# PCoA --------------------------------------------------------------------

#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- dorcierr_final%>%
  select(-c(1:5))

veg_bray <- vegdist(microb_f, "bray") #Bray-curtis distances

pc_scores<-pcoa(veg_bray) #Calculating scores

#This plots Eigenvalues
#They will allow you to choose the best axes to show how your data varies
ylimit = c(0, 1.1*max(pc_scores$values$Relative_eig))

Eigan <- barplot(pc_scores$values$Relative_eig[1:10], ylim= ylimit)
# Add values to the bars
text(x = Eigan, y = pc_scores$values$Relative_eig[1:10], label = pc_scores$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

# Plot using ggplot2
pco_scores <- as.data.frame(pc_scores$vectors)

ggplot(pco_scores, mapping = aes(Axis.1, Axis.2, col = dorcierr_final$Organism, 
                                 shape = factor(dorcierr_final$DayNight))) +
  geom_point(stat = "identity")+
  scale_shape_manual(values = c(1,19)) +
  scale_color_brewer(type = "qual", palette = 2)

# PERMANOVA ---------------------------------------------------------------

adonis(microb_f ~ Organismal_clade, moorea_wdf_fix, perm=1000, method="bray", set.seed(100))

pairwiseAdonis::pairwise.adonis(microb_f, dorcierr_wdf$Organism, p.adjust.m = "BH")
# Writing CSVs ------------------------------------------------------------

write_csv(moorea_wdf, "moorea_metab_FeatureID_samples.csv")

write_csv(both_runs, "Moorea_2017_new_and_old_runs.csv")
write_delim(scan_IDs, "Feature_ID.txt", delim = " ", col_names = FALSE)

write_csv(both_runs_by_old, "new_by_old_AllMeta.csv")

write_csv(dorcierr, "Dorcierr_metab_wdf.csv")

