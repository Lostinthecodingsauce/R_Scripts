rm(list = ls())


# Libraries needed --------------------------------------------------------

#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)

#PCoA, PERMANOVA
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)

# Loading in Data and Joining Tables ---------------------------------------------------------

analog_hits <- readxl::read_excel("Analog_hits.xlsx")%>%
  rename(Scan_Number = 'Scan_analog')
true_hits <- readxl::read_excel("library_hits_true.xlsx")%>%
  rename(Scan_Number = 'Scan')
old_hits <- readxl::read_excel("Moorea_SPEDOM_MSMS_Polished_April2018 (Autosaved).xlsx")%>%
  rename(Scan_Number = 'row ID')
Node_info <- readxl::read_excel("Node_Info.clustersummary.xlsx")

## Joining tables

new_run <- full_join(full_join(Node_info, true_hits, "Scan_Number"),analog_hits, "Scan_Number")

both_runs <- full_join(new_run, old_hits, "Scan_Number")

both_runs_by_old <- right_join(new_run, old_hits, "Scan_Number")

scan_IDs <- as.data.frame(old_hits$Scan_Number)

# Select raw data for Daniels analysis
raw_samples <- both_runs_by_old%>%
  select(c('Feature_Name','% of Total(D_WA_1_T0D)':'% of Total(D_Blank_DI)'))

# Transposing and selecting only sample data
samples <- both_runs_by_old%>%
  select(c('Feature_Name', '% of Total(D_WA_1_T0D)':'% of Total(R_TR_3_TFN)'))%>%
  gather(sample_ID, RA, 2:278)%>%
  spread(Feature_Name, RA)

# Some Samples had - instead of _
samples$sample_ID <- gsub("-", "_", samples$sample_ID)
samples$sample_ID <- gsub("% of Total\\(", "", samples$sample_ID)
samples$sample_ID <- gsub(")", "", samples$sample_ID)

##Converts all RA columns to numeric and percentages
moorea_percentages <- as.data.frame(sapply(samples[2:10027], function(x) as.numeric(gsub("%", "", x))/100))

#Transforms RA columns to asin(sqrt)
moorea_asin_sqrt <- as.data.frame(sapply(moorea_percentages, function(x) asin(sqrt(x))))

#Data matrix with both RA and asin(sqrt(asin)) data
moorea_asin_samples <- moorea_asin_sqrt%>%
  add_column("sample_ID" = samples$sample_ID, .before = 1)
moorea_percent_samples <- moorea_percentages%>%
  add_column("sample_ID" = samples$sample_ID, .before = 1)
moorea_percent_transformed <- full_join(moorea_percent_samples, moorea_asin_samples, by = "sample_ID",  suffix = c(".RA", ".asin(sqrt)"))

## renaming columns and row metadata
## from here on columns 6:10031 are RA and 10032:20053 are transformed
moorea_wdf <-moorea_percent_transformed%>%
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
# If low on memory. run the lines below: 
# write_csv(dorcierr_wdf, "dorcier_wdf_lowmemory.csv")
# rm(list = ls())
# dorcierr_wdf <- read_csv("dorcier_wdf_lowmemory.csv")

dorcierr_final <- dorcierr_wdf%>%
  filter(Timepoint == "TF")

## Need to take the mean of T0 to subtract from TF because of the differences in replicate #
dorcierr_transformations <- dorcierr_wdf%>%
  select(c(1:10031))%>% 
  gather(feature_name, RA, 6:10031)%>%
  spread(Timepoint, RA)

dorcierr_t0 <- dorcierr_transformations%>%
  select(-c(TF))%>%
  filter(!Replicate == 3,
         !Replicate == 4)%>%
  group_by(Organism, DayNight, feature_name)%>%
  summarize_if(is.numeric, mean)

dorcierr_newt0 <- left_join(dorcierr_transformations, dorcierr_t0, by = c("Organism", "DayNight", "feature_name"))%>%
  select(-c(T0.x, Replicate.y))%>%
  rename(Replicate = 'Replicate.x',
         T0 = 'T0.y')

dorcierr_newt0$change <- (dorcierr_newt0$TF - dorcierr_newt0$T0)
  
dorcierr_transformed <- dorcierr_newt0%>%
  select(-c(T0,TF))%>%
  spread(feature_name, change)


# PCoA --------------------------------------------------------------------

#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- dorcierr_transformed%>%
  select(c(5:10030))

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

pairwiseAdonis::pairwise.adonis(microb_f, dorcierr_final$Organism, p.adjust.m = "BH")


# ANOVA and Post-Hoc -------------------------------------------------------------------

anova_variables <- c("Organism", "DayNight", "Organism*DayNight")

two_way_model <- sapply(dorcierr_transformed[5:10030], function(x) aov(x ~ dorcierr_transformed[["Organism"]]*dorcierr_transformed[["DayNight"]]))

p_values <- as.data.frame(sapply(dorcierr_transformed[5:10030], function(x) summary(
  aov(x ~ dorcierr_transformed[["Organism"]]*dorcierr_transformed[["DayNight"]]))[[1]][1:3,'Pr(>F)']))

# Add test names as the first column
anova_p_values <- p_values%>%
  add_column(anova_variables, .before = 1)

# Now the data is made Tidy and we filter to only significant values
anova_dr_tidy <- anova_p_values%>%
  gather(Clade, F_value, 2:10027)%>%
  filter(F_value, F_value <0.05)

organism_anova <- anova_dr_tidy%>%
  filter(anova_variables == "Organism")

organism_sig_clades <- as.vector(organism_anova$Clade)

Tukey <- sapply(dorcierr_transformed[5:10030], function(x) TukeyHSD(
  aov(x ~ dorcierr_transformed[["Organism"]]*dorcierr_transformed[["DayNight"]])))

tukey_sig <- as.data.frame(Tukey)

## Tukey Organism
tukey_organism <- tukey_sig%>%
  rownames_to_column(var = "Anova_test")%>%
  filter(Anova_test == 'dorcierr_transformed[["Organism"]]')%>%
  select(-1)%>%
  gather(feature_name, p_value, 1:ncol(tukey_sig))%>%
  filter(feature_name %in% c(organism_sig_clades))

tukey_organism$p_value <-  tukey_organism$p_value%>%
  gsub("c\\(", "", .)%>%
  gsub(")", "", .)

tukey_organism_ps <- tukey_organism%>%
  separate(p_value, paste("c", 1:60, sep = "."), sep = ",")


# Finding Largest "Player"  ---------------------------------------------------

dorcierr_template_day <- dorcierr_transformed%>%
  gather(feature_name, RA, 5:ncol(dorcierr_transformed))%>%
  filter(feature_name %in% c(organism_sig_clades))%>%
  filter(DayNight == "Day")%>%
  select(-c(Experiment, DayNight, Replicate))%>%
  mutate(RA = as.numeric(RA))

dorcierr_means <- dorcierr_template_day%>%
  group_by(Organism, feature_name)%>%
  summarize_if(is.numeric, mean)%>%
  spread(Organism, RA)%>%
  group_by(feature_name)

dorcierr_means$max <- apply(dorcierr_means[2:7], 1, max)



# Writing CSVs ------------------------------------------------------------

write_csv(moorea_wdf, "moorea_metab_FeatureID_samples.csv")

write_csv(both_runs, "Moorea_2017_new_and_old_runs.csv")
write_delim(scan_IDs, "Feature_ID.txt", delim = " ", col_names = FALSE)

write_csv(both_runs_by_old, "new_by_old_AllMeta.csv")

write_csv(dorcierr, "Dorcierr_metab_wdf.csv")

