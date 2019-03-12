# rm(list = ls())


# Libraries needed --------------------------------------------------------

#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)

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
  rename(Scan_Number = 'Scan_old')
Node_info <- readxl::read_excel("Node_Info.clustersummary.xlsx")

Node_check <- left_join(old_hits, Node_info, by = "Scan_Number")%>%
  select(c('Scan_Number', 'componentindex', 'componentindex_old'))

## Joining tables

new_run <- full_join(full_join(Node_info, true_hits, "Scan_Number"),analog_hits, "Scan_Number")

both_runs <- full_join(new_run, old_hits, "Scan_Number")

both_runs_by_old <- right_join(new_run, old_hits, "Scan_Number")

scan_IDs <- as.data.frame(old_hits$Scan_Number)


# Designing tables for future analyses ------------------------------------

# Select Feature name and Network it came from
network <- both_runs_by_old%>%
  select(c(Feature_Name,`componentindex_old`, `Compound_Name_analog`, `LibraryID`))%>%
  rename(feature_name = `Feature_Name`,
         cluster = `componentindex_old`,
         compound = `Compound_Name_analog`)

# Transposing and selecting only sample data
sample_raw <- both_runs_by_old%>%
  select(c('Feature_Name', '% of Total(D_WA_1_T0D)':'% of Total(R_TR_3_TFN)'))

sample_percentages <- as.data.frame(sapply(sample_raw[2:ncol(sample_raw)], function(x) as.numeric(gsub("%", "", x))/100))%>%
  add_column(Feature_Name = sample_raw$Feature_Name, .before = 1)

samples <- sample_percentages%>%
  add_column(max = apply(sample_raw[2:ncol(sample_raw)],1,max))%>%
  filter(max, max > 0.0001)%>%
  select(-"max")%>%
  gather(sample_ID, RA, 2:278)%>%
  spread(Feature_Name, RA)

# Some Samples had - instead of _
samples$sample_ID <- samples$sample_ID%>%
  gsub("-", "_", .)%>%
  gsub("% of Total\\(", "", .)%>%
  gsub(")", "", .)

#Transforms RA columns to asin(sqrt)
moorea_asin_sqrt <- as.data.frame(sapply(samples[2:ncol(samples)], function(x) asin(sqrt(x))))

#Data matrix with both RA and asin(sqrt(asin)) data
moorea_asin_samples <- moorea_asin_sqrt%>%
  add_column("sample_ID" = samples$sample_ID, .before = 1)
moorea_percent_transformed <- full_join(samples, moorea_asin_samples, by = "sample_ID",  suffix = c(".RA", ".asin(sqrt)"))

## renaming columns and row metadata
## from here on columns 6:10031 are RA and 10032:20053 are transformed
moorea_wdf <-moorea_percent_transformed%>%
  mutate(sample_ID = case_when(sample_ID == "SPIFFy_7_B 2" ~ "SPIFFy_6_B",
                               TRUE ~ as.character(sample_ID)))%>%
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

 spiffy <- moorea_wdf%>%
   filter(Experiment == "SPIFFy")%>%
   dplyr::select(-c("Timepoint", "DOM_source"))
   


# Spiffy workshopping ---------------------------------------------------------
spiffy_wdf <- spiffy%>%
   add_column(reef_area = spiffy$Organism, .before = 3)%>%
   rename(`Site` =`Organism`)%>%
   mutate(reef_area = case_when(reef_area == 1 ~ "bay",
                                reef_area == 2 ~ "bay",
                                reef_area == 3 ~ "bay",
                                reef_area == 4 ~ "bay",
                                reef_area == 5 ~ "bay",
                                reef_area == 6 ~ "backreef",
                                reef_area == 7 ~ "backreef",
                                reef_area == 8 ~ "backreef",
                                reef_area == 9 ~ "backreef",
                                reef_area == 10 ~ "backreef",
                                reef_area == 11 ~ "backreef",
                                reef_area == 12 ~ "backreef",
                                TRUE ~ "Forereef"))


spiffy_no_bay <-spiffy_wdf%>%
  filter(!reef_area == "bay")

aov_spiffy <- as.data.frame(sapply(spiffy_no_bay[8859:ncol(spiffy_no_bay)], function(x) summary(aov(x ~ spiffy_no_bay[["reef_area"]]))[[1]][1,'Pr(>F)']))

anova_spiffy <- aov_spiffy%>%
  rownames_to_column(var = "feature_name")%>%
  rename(`f_value`= `sapply(spiffy_no_bay[8859:ncol(spiffy_no_bay)], function(x) summary(aov(x ~ spiffy_no_bay[["reef_area"]]))[[1]][1, "Pr(>F)"])`)

anova_spiffy$FDR_f <- p.adjust(anova_spiffy$f_value, method = "BH")

anova_sigs_spiffy <- anova_spiffy%>%
  filter(FDR_f, FDR_f < 0.05)

spiffy_sig_names <-as.vector(anova_sigs_spiffy$feature_name)

spiffy_only_sigs <- spiffy_no_bay%>%
  select(c(1:4, spiffy_sig_names))%>%
  group_by(reef_area)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(feature_name, RA, 2:332)%>%
  spread(reef_area, RA)

spiffy_means_only <-spiffy_only_sigs%>%
  select(-1)

spiffy_only_sigs$MaxNames <- colnames(spiffy_means_only)[max.col(spiffy_means_only, ties.method="first")]

spiffy_only_sigs$max <- apply(spiffy_only_sigs[2:3], 1, max)

spiffy_two_locations_pvalues_max <- spiffy_only_sigs%>%
  add_column(f_value = anova_sigs_spiffy$f_value)%>%
  add_column(FDR_f = anova_sigs_spiffy$FDR_f)%>%
  filter(FDR_f, FDR_f < 0.05)

spiffy_two_locations_pvalues_max$feature_name <-spiffy_two_locations_pvalues_max$feature_name%>%
  gsub(".asin\\(sqrt)", "", .)



spiffy_twolocals_network <- left_join(
  spiffy_two_locations_pvalues_max, network, by = "feature_name"
)


spiffy_rare <- spiffy_no_bay%>%
  select(Site, reef_area, 5:8858)%>%
  unite(site_name, c(Site, reef_area), sep = "_", remove = TRUE)%>%
  group_by(site_name)%>%
  summarize_if(is.numeric, mean)%>%
  gather(feature_name, RA, 2:8855)%>%
  spread(site_name, RA)%>%
  add_column(max = apply(.[2:ncol(.)], 1, max))%>%
  filter(max, max < 0.0001)
  
spiffy_rare$feature_name <- gsub(".RA", "", spiffy_rare$feature_name)  
  
spiffy_rare_byreeflocal <- spiffy_rare%>%
  ungroup()%>%
  select(-max)%>%
  gather(site_name, RA, 2:16)%>%
  separate(site_name, c("site", "reef_area"), sep = "_")%>%
  select(-site)%>%
  group_by(feature_name, reef_area)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(reef_area, RA)
  

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
  select(c(1:8858))%>% 
  gather(feature_name, RA, 6:8858)%>%
  spread(Timepoint, RA)

dorcierr_t0 <- dorcierr_transformations%>%
  select(-c(TF))%>%
  filter(!Replicate == 3,
         !Replicate == 4)%>%
  group_by(Organism, DayNight, feature_name)%>%
  summarize_if(is.numeric, mean)

dorcierr_newt0 <- left_join(dorcierr_transformations, dorcierr_t0, by = c("Organism", "DayNight", "feature_name"))%>%
  select(-c(T0.x))%>%
  rename(T0 = 'T0.y')

dorcierr_newt0$change <- (dorcierr_newt0$TF - dorcierr_newt0$T0)
  
dorcierr_transformed <- dorcierr_newt0%>%
  select(-c(T0,TF))%>%
  spread(feature_name, change)

dorcierr_day_transformed <- dorcierr_transformed%>%
  filter(DayNight == "Day")


dorcierr_exudates_day <- dorcierr_wdf%>%
  filter(DayNight == "Day",
         Timepoint == "T0",
         !Replicate == 3,
         !Replicate == 4)

# PCoA --------------------------------------------------------------------

#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- dorcierr_final%>%
  select(c(10032:ncol(dorcierr_final)))

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

ggplot(pco_scores, mapping = aes(Axis.1, Axis.2, col = dorcierr_wdf$Organism, 
                                 shape = factor(dorcierr_wdf$DayNight))) +
  geom_point(stat = "identity")+
  scale_shape_manual(values = c(1,19)) +
  scale_color_brewer(type = "qual", palette = 2)

# NMDS --------------------------------------------------------------------

set.seed(1234)
microb_NMDS <- metaMDS(microb_f)

# produce Shepard plot
# large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. 
stressplot(microb_NMDS)

# plot the NMDS
# open circles are communities, red crosses are species
plot(microb_NMDS)

ordiplot(microb_NMDS, type = "n")
orditorp(microb_NMDS, display = "sites", air = 0.25)


# create a plot with convex hulls connecting vertices of the points made by a treatment

treatment <- c(microb_f_abun$Organism)
samples <- c(microb_f_abun$Label)
##Different Plotting options. Comment in to plot

# ordiplot(microb_NMDS, type = "n",
#          main = "NMDS of Incubation Vessels")
# ordihull(microb_NMDS, groups = treatment, draw = "polygon", label = F)
# #orditorp(microb_NMDS, display = "sites", labels = samples, air = 0.05)
# orditorp(microb_NMDS, display = "sites", air = 0.2)
# 
# 
# # spider plot
# ordiplot(microb_NMDS, type = "n",
#          main = "NMDS of Incubation Vessels")
# orditorp(microb_NMDS, display = "sites", air = 0.2)
# ordispider(microb_NMDS, groups = treatment)
# 
# 
# # ellipse plot
# ordiplot(microb_NMDS, type = "n",
#          main = "NMDS of Incubation Vessels")
# orditorp(microb_NMDS, display = "sites", air = 0.2)
# ordiellipse(microb_NMDS, groups = treatment)
# 
# 
# 
# # MST plot
# 
# dist_microb <- vegdist(microb)
# clust_microb <- hclust(dist_microb, method = "complete")
# 
# 
# ordiplot(microb_NMDS, type = "n",
#          main = "NMDS of Incubation Vessels")
# orditorp(microb_NMDS, display = "sites", air = 0.2)
# ordicluster(microb_NMDS, cluster = clust_microb)


# plot NMDS output in ggplot

df.scores <- as.data.frame(scores(microb_NMDS)) #Using the scores function from vegan to extract the site scores and convert to a data.frame

df.scores$site <- rownames(df.scores) # create a column of site names, from the rownames of data.scores

# add grouping variables
df.scores$group <- dorcierr_final$Organism
df.scores$DNA_source <- dorcierr_final$DayNight
# df.scores$Organism <- dorcierr_final$Timepoint


species.scores <- as.data.frame(scores(microb_NMDS, "species")) #Using the scores function from vegan to extract the species scores and convert to a data.frame

species.scores$species <- rownames(species.scores)



# plot NMDS

ggplot(data = df.scores, mapping = aes(NMDS1, NMDS2, col = dorcierr_final$Organism, 
                     shape = factor(dorcierr_final$DayNight))) +
  geom_point() +
  coord_equal() +
  scale_shape_manual(values = c(0, 19)) +
  theme_bw()

# PERMANOVA ---------------------------------------------------------------

adonis(microb_f ~ Organismal_clade, moorea_wdf_fix, perm=1000, method="bray", set.seed(100))

pairwiseAdonis::pairwise.adonis(microb_f, dorcierr_final$Organism, p.adjust.m = "BH")

# One-Way Anova -----------------------------------------------------------
p_values_oneway <- as.data.frame(sapply(dorcierr_day_transformed[5:ncol(dorcierr_day_transformed)], function(x) summary(
  aov(x ~ dorcierr_day_transformed[["Organism"]]))[[1]][1,'Pr(>F)']))

# Now the data is made Tidy and we filter to only significant values
anova_dr_tidy <- p_values_oneway%>%
  rownames_to_column(var = "feature_name")%>%
  rename(f_values = 2)%>%
  filter(f_values, f_values < 0.05)

oneway_anova_sig <- as.vector(anova_dr_tidy$feature_name)

sigs_only_day <- dorcierr_day_transformed%>%
  select(c(1:4, oneway_anova_sig))%>%
  add_column("grouping" = dorcierr_day_transformed$Organism, .before = 2)%>%
  mutate(grouping = case_when(grouping == "Pocillopora verrucosa" ~ "Coral",
                              grouping == "Porites lobata" ~ "Coral",
                              grouping == "Dictyota" ~ "Fleshy Algae",
                              grouping == "Turf" ~ "Fleshy Algae",
                              TRUE ~ as.character(grouping)))

write_csv(sigs_only_day, "dorcierr_day_forDunnetts.csv")

## Post-Hoc tests run in JMP cause R is stupid and I hate it
dunnetts_pvalue <- read_csv("~/Documents/SDSU/Moorea_2017/dorcierr/Stats/dorcierr_day_Dunnetts.csv")
  
tukey_pvalue <-read_csv("~/Documents/SDSU/Moorea_2017/dorcierr/Stats/dorcierr_day_Tukey.csv")

dunnetts_sig <- dunnetts_pvalue%>%
  rename(`FDR_p` =`FDR Adj PValue`)%>%
  filter(FDR_p, FDR_p < 0.05)

# Different Bins ----------------------------------------------------------
## making a primary producer pool data frame and determining what increases across all five primary producers
primary_producer_pool_almost <- dunnetts_sig%>%
  select(-c(1,4,5))%>%
  mutate(FDR_p = case_when(!FDR_p == "NA" ~ 1))%>%
  spread(Level, FDR_p)

primary_producer_pool <- primary_producer_pool_almost%>%
  add_column(sum = apply(primary_producer_pool_almost[2:6], 1, sum))%>%
  filter(sum, sum == 5)

primary_producer_compounds <- as.vector(primary_producer_pool$Y)  

pp_increase <- sigs_only_day%>%
  select(1, c(primary_producer_compounds))%>%
  group_by(Experiment)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-1)%>%
  gather(feature_name, change, 1:100)%>%
  filter(change, change > 0.00)


pp_decrease <- sigs_only_day%>%
  select(1, c(primary_producer_compounds))%>%
  group_by(Experiment)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  select(-1)%>%
  gather(feature_name, change, 1:100)%>%
  filter(change, change < 0.00)

pp_decrease$feature_name <- pp_decrease$feature_name%>%
  gsub(".RA", "", .)

pp_increase$feature_name <- pp_increase$feature_name%>%
  gsub(".RA", "", .)

## Trying to match pp_decrease to spiffy
spiffy_pp_decrease <- inner_join(pp_decrease, spiffy_two_locations_pvalues_max, by = "feature_name")
spiffy_pp_increase <- inner_join(pp_increase, spiffy_two_locations_pvalues_max, by = "feature_name")

overlap_primary_producers_increase <- spiffy_pp_increase%>%
  add_column("increased_dorcierr" = "yes", .after = 1)

overlap_primary_producers_decrease <-spiffy_pp_decrease%>%
  add_column("increased_dorcierr" = "no", .after = 1)

overlap_primary_producers <- left_join(
  bind_rows(overlap_primary_producers_increase, overlap_primary_producers_decrease), network, by = "feature_name")%>%
  select(-c(11,12))%>%
  separate(feature_name, "component_ID", sep = "_", remove = FALSE)
  
## Possible candidates of super labile compounds. 
## Decrease in Dorcierr and found in low abundance in Spiffy (RA< 0.0001) in back/forereef
spiffy_pvalues_max$feature_name <-spiffy_pvalues_max$feature_name%>%
  gsub(".asin\\(sqrt)", "", .)

noverlap_primary_producers_decrease <- inner_join(pp_decrease, spiffy_rare_byreeflocal, by = "feature_name")

super_labile_primary_producer_candidates <- left_join(noverlap_primary_producers_decrease, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  separate(feature_name, "component_ID", sep = "_", remove = FALSE)

## Possible candidates of super recalcitrant compounds. Increase in Dorcierr and found in Spiffy
# super_recalcitrant_primary_producers


# BULK REEF EXUDATES [NEEDS WORK. NOT FILTERED BY SIG SPECIES>CONTROL (Level)]------------------------------------------------------
## bulk_reef_exudates is a list of everything which is siginificant by dunnetts
bulk_reef_exudates <- as.vector(dunnetts_sig%>%
  group_by(Y)%>%
  summarize_if(is.numeric, mean)%>%
    select(Y,2)%>%
    spread(Y, 2)%>%
    select(-c(primary_producer_compounds))%>%
    gather(Y, x))$Y
  
## sigs_only_day is all feature RA in the day which were significant by an anova
## there are 1380 unique features which were significant in any one species from water control BY GROUPING
bulk_reef_exudates_day_increase <- sigs_only_day%>%
  select(grouping, c(bulk_reef_exudates))%>%
  group_by(grouping)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  gather(feature_name, change, 2:1381)%>%
  filter(change, change > 0.00)

bulk_reef_exudates_day_decrease <- sigs_only_day%>%
  select(grouping, c(bulk_reef_exudates))%>%
  group_by(grouping)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  gather(feature_name, change, 2:1381)%>%
  filter(change, change < 0.00)

bulk_reef_exudates_day_decrease$feature_name <- bulk_reef_exudates_day_decrease$feature_name%>%
  gsub(".RA", "", .)

bulk_reef_exudates_day_increase$feature_name <- bulk_reef_exudates_day_increase$feature_name%>%
  gsub(".RA", "", .)

## sigs_only_day is all feature RA in the day which were significant by an anova
## there are 1380 unique features which were significant in any one species from water control BY SPECIES
species_bre_day_increase <- sigs_only_day%>%
  select(Organism, c(bulk_reef_exudates))%>%
  filter(!Organism == "Water control")%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%  
  gather(feature_name, change, 2:1281)%>%
  filter(change, change > 0.00)

species_bre_day_decrease <- sigs_only_day%>%
  select(Organism, c(bulk_reef_exudates))%>%
  filter(!Organism == "Water control")%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(feature_name, change, 2:1281)%>%
  filter(change, change < 0.00)
  
species_bre_day_increase$feature_name <- gsub(".RA", "", species_bre_day_increase$feature_name)

species_bre_day_decrease$feature_name <- gsub(".RA", "", species_bre_day_decrease$feature_name)

## By species decrease by network
sbre_day_decrease_network <- left_join(species_bre_day_decrease, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  add_column("counting" = as.numeric("1"))%>%
  spread(Organism, counting)
  
sbre_day_decrease_network[is.na(sbre_day_decrease_network)] <- 0

sbreddecrease_features <-sbre_day_decrease_network%>%
  group_by(feature_name, cluster)%>%
  summarise_at(vars(4:8), sum)%>%
  ungroup()

# Making all_networks for notes -------------------------------------------
## organismal subsets of reef_biomarkers
## coral
bre_day_coral_decrease <- left_join(bulk_reef_exudates_day_decrease, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  filter(grouping == "Coral")%>%
  add_column("coral" = as.numeric("1"))%>%
  group_by(cluster)%>%
  summarise_at(vars(coral), sum)%>%
  ungroup()%>%
  arrange(desc(coral))

bre_day_coral_increase <- left_join(bulk_reef_exudates_day_increase, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  filter(grouping == "Coral")%>%
  add_column("coral" = as.numeric("1"))%>%
  group_by(cluster)%>%
  summarise_at(vars(coral), sum)%>%
  ungroup()%>%
  arrange(desc(coral))

coral_networks <- full_join(bre_day_coral_decrease, bre_day_coral_increase, by = "cluster", suffix = c("_decrease", "_increase"))

## CCA
bre_day_CCA_decrease <- left_join(bulk_reef_exudates_day_decrease, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  filter(grouping == "CCA")%>%
  add_column("cca" = as.numeric("1"))%>%
  group_by(cluster)%>%
  summarise_at(vars(cca), sum)%>%
  ungroup()%>%
  arrange(desc(cca))

bre_day_CCA_increase <- left_join(bulk_reef_exudates_day_increase, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  filter(grouping == "CCA")%>%
  add_column("cca" = as.numeric("1"))%>%
  group_by(cluster)%>%
  summarise_at(vars(cca), sum)%>%
  ungroup()%>%
  arrange(desc(cca))

cca_networks <- full_join(bre_day_CCA_decrease, bre_day_CCA_increase, network, by = "cluster", suffix = c("_decrease", "_increase"))

## All Algae
bre_day_algae_decrease <- left_join(bulk_reef_exudates_day_decrease, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  filter(!grouping == "Coral",
         !grouping == "Water control")%>%
  add_column("algae" = as.numeric("1"))%>%
  group_by(cluster)%>%
  summarise_at(vars(algae), sum)%>%
  ungroup()%>%
  arrange(desc(algae))

bre_day_algae_increase <- left_join(bulk_reef_exudates_day_increase, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  filter(!grouping == "Coral",
         !grouping == "Water control")%>%
  add_column("algae" = as.numeric("1"))%>%
  group_by(cluster)%>%
  summarise_at(vars(algae), sum)%>%
  ungroup()%>%
  arrange(desc(algae))

algae_networks <- full_join(bre_day_algae_decrease, bre_day_algae_increase, by = "cluster", suffix = c("_decrease", "_increase"))

## Fleshy Algae
bre_day_fleshy_decrease <- left_join(bulk_reef_exudates_day_decrease, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  filter(grouping == "Fleshy Algae")%>%
  add_column("fleshy" = as.numeric("1"))%>%
  group_by(cluster)%>%
  summarise_at(vars(fleshy), sum)%>%
  ungroup()%>%
  arrange(desc(fleshy))

bre_day_fleshy_increase <- left_join(bulk_reef_exudates_day_increase, network, by = "feature_name")%>%
  select(-c('compound', 'LibraryID'))%>%
  filter(grouping == "Fleshy Algae")%>%
  add_column("fleshy" = as.numeric("1"))%>%
  group_by(cluster)%>%
  summarise_at(vars(fleshy), sum)%>%
  ungroup()%>%
  arrange(desc(fleshy))

fleshy_network <- full_join(bre_day_fleshy_decrease, bre_day_fleshy_increase, by = "cluster", suffix = c("_decrease", "_increase"))


all_networks <- full_join(
    full_join(
      full_join(coral_networks, fleshy_network, by = "cluster"), 
      cca_networks, by = "cluster"), 
    algae_networks, by = "cluster")

all_networks[is.na(all_networks)] <- 0

all_networks$max <- apply(all_networks[2:ncol(all_networks)], 1, sum)




# Dorcierr Exudation ------------------------------------------------------

exudates_day_dunnetts <- read_csv("~/Documents/SDSU/Moorea_2017/dorcierr/Stats/Dorcierr_exudates_day_Dunnetts.csv")%>%
  rename(`FDR_p` = `FDR Adj PValue`)

sig_exudates <- exudates_day_dunnetts%>%
  filter(FDR_p, FDR_p < 0.05)%>%
  select(-c(1,4,5))%>%
  add_column("binary" = as.numeric("1"))%>%
  rename(`feature_name` = `Y`)

sig_exudates$feature_name <- gsub(".asin\\(sqrt)", "", sig_exudates$feature_name)

sig_exudates_spread <- sig_exudates%>%
  select(-FDR_p)%>%
  spread(Level, binary)%>%
  add_column("coral" = .$`Pocillopora verrucosa`+ .$`Porites lobata`)%>%
  add_column("calcifiers" = .$coral + .$CCA)%>%
  add_column("fleshy" = .$Turf + .$Dictyota)%>%
  add_column("algae"= .$fleshy + .$CCA)%>%
  add_column("primary_producers" = .$coral +.$algae)%>%
  add_column("pv_cca" = .$`Pocillopora verrucosa`+ .$CCA)%>%
  add_column("pv_turf" = .$`Pocillopora verrucosa`+ .$Turf)%>%
  add_column("pv_dictyota" = .$`Pocillopora verrucosa`+ .$Dictyota)%>%
  add_column("pl_cca" = .$`Porites lobata`+ .$CCA)%>%
  add_column("pl_turf" = .$`Porites lobata`+ .$Turf)%>%
  add_column("pl_dictyota" = .$`Porites lobata`+ .$Dictyota)%>%
  add_column("cca_turf" = .$CCA+ .$Turf)%>%
  add_column("cca_dictyota" = .$CCA + .$Dictyota)%>%
  add_column("coral_turf" = .$coral + .$Turf)%>%
  add_column("coral_dictyota" = .$coral + .$Dictyota)%>%
  add_column("pv_fleshy" = .$fleshy + .$`Pocillopora verrucosa`)%>%
  add_column("pl_fleshy" = .$fleshy + .$`Porites lobata`)%>%
  add_column("not_PV" = .$algae + .$`Porites lobata`)%>%
  add_column("not_PL" = .$algae+ .$`Pocillopora verrucosa`)%>%
  add_column("not_CCA" = .$coral + .$fleshy)%>%
  add_column("not_turf" = .$calcifiers + .$Dictyota)%>%
  add_column("not_Dicyota" = .$calcifiers + .$Turf)

sig_exudates_spread[is.na(sig_exudates_spread)] <- 0

sig_exudates_numeric <- sig_exudates_spread%>%
  select(-1)
  
sig_exudates_numeric$MaxNames <- colnames(sig_exudates_numeric)[max.col(sig_exudates_numeric, ties.method="first")]  

exudates_ranked <- sig_exudates_numeric%>%
  add_column(sig_exudates_spread$Y, .before = 1)%>%
  select(1, 29)%>%
  rename(`feature_name` = `sig_exudates_spread$Y`)

exudate_networks_ranked <- left_join(exudates_ranked, network, by = "feature_name")%>%
  select(-c(4,5))%>%
  separate(feature_name, "component_ID", sep = "_", remove = FALSE)

## networks for features
exudates_network_features <- left_join(sig_exudates, network, by = "feature_name")%>%
  select(-c(3,6,7))%>%
  group_by(feature_name, cluster)%>%
  summarize_at(vars(3), sum)

# Two Way-ANOVA and Post-Hoc -------------------------------------------------------------------

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
## Daytime
dorcierr_template_day <- dorcierr_transformed%>%
  gather(feature_name, RA, 5:ncol(dorcierr_transformed))%>%
  filter(feature_name %in% c(organism_sig_clades))%>%
  filter(DayNight == "Day")%>%
  select(-c(Experiment, DayNight, Replicate))%>%
  mutate(RA = as.numeric(RA))

dorcierr_template_day$feature_name <- gsub(".RA", "", dorcierr_template_day$feature_name)

dorcierr_means_day <- dorcierr_template_day%>%
  group_by(Organism, feature_name)%>%
  summarize_if(is.numeric, mean)%>%
  spread(Organism, RA)%>%
  group_by(feature_name)

dorcierr_means_day$max <- apply(dorcierr_means_day[2:7], 1, max)

means_only_day <- dorcierr_means_day%>%
  ungroup()%>%
  select(-c(1, max))

dorcierr_means_day$MaxNames <- colnames(means_only_day)[max.col(means_only_day, ties.method="first")]

dorcierr_meanmax_organism_day <- left_join(dorcierr_means_day, network, by = "feature_name")%>%
  arrange(MaxNames, cluster)%>%
  filter(!compound == "NA")
  
organism_day_neg <- dorcierr_meanmax_organism_day%>%
  filter(max, max < 0.00)

organism_day_pos <- dorcierr_meanmax_organism_day%>%
  filter(max, max > 0.00)

## Nighttime
dorcierr_template_night <- dorcierr_transformed%>%
  gather(feature_name, RA, 5:ncol(dorcierr_transformed))%>%
  filter(feature_name %in% c(organism_sig_clades))%>%
  filter(DayNight == "Night")%>%
  select(-c(Experiment, DayNight, Replicate))%>%
  mutate(RA = as.numeric(RA))

dorcierr_template_night$feature_name <- gsub(".RA", "", dorcierr_template_night$feature_name)

dorcierr_means_night <- dorcierr_template_night%>%
  group_by(Organism, feature_name)%>%
  summarize_if(is.numeric, mean)%>%
  spread(Organism, RA)%>%
  group_by(feature_name)

dorcierr_means_night$max <- apply(dorcierr_means_night[2:7], 1, max)

means_only_night <- dorcierr_means_night%>%
  ungroup()%>%
  select(-c(1, max))

dorcierr_means_night$MaxNames <- colnames(means_only_night)[max.col(means_only_night, ties.method="first")]

dorcierr_meanmax_organism_night <- left_join(dorcierr_means_night, network, by = "feature_name")%>%
  arrange(MaxNames, cluster)%>%
  filter(!compound == "NA")

organism_night_neg <- dorcierr_meanmax_organism_night%>%
  filter(max, max < 0.00)

organism_night_pos <- dorcierr_meanmax_organism_night%>%
  filter(max, max > 0.00)

test <- inner_join(organism_night_pos, organism_day_neg, by = "feature_name", suffix = c(".night", ".day"))%>%
  select(-c(2:8, 12:18))

# Super Labile Primary Producer Compounds ---------------------------------
labile_pp_RA <- super_labile_primary_producer_candidates

labile_pp_RA$feature_name <- paste0(labile_pp_RA$feature_name, ".RA")

labile_pp_RA_vector <- as.vector(labile_pp_RA$feature_name)

super_labile_pp_42 <- as.vector(labile_pp_RA%>%
                                  filter(cluster == "42"))$feature_name

super_labile_pp_90 <-as.vector(labile_pp_RA%>%
                                 filter(cluster == "90"))$feature_name

super_labile_pp_623 <- as.vector(labile_pp_RA%>%
                                   filter(cluster == "623"))$feature_name

dorcierr_labile_pp_42 <- dorcierr_wdf%>%
  select(1:5, c(super_labile_pp_42))

dorcierr_labile_pp_90 <- dorcierr_wdf%>%
  select(1:5, c(super_labile_pp_90))

dorcierr_labile_pp_623 <- dorcierr_wdf%>%
  select(1:5, c(super_labile_pp_623))


dorcierr_labile_pp_general <- dorcierr_wdf%>%
  select(1:5, c(labile_pp_RA_vector))

write_csv(dorcierr_labile_pp_general, "dorcierr_RA_pp_all.csv")
# Significantly eaten -----------------------------------------------------
sbreddecrease_features_dotRA <- sbreddecrease_features

sbreddecrease_features_dotRA$feature_name <- paste0(sbreddecrease_features_dotRA$feature_name, ".RA")

network_596_eaten <- as.vector(sbreddecrease_features_dotRA%>%
                                 filter(cluster == "596"))$feature_name

dorcierr_labile_calcifiers_596 <- dorcierr_wdf%>%
  select(1:5, c(network_596_eaten))

day_transformed_transposed <- dorcierr_day_transformed%>%
  gather(feature_name, RA, 5:8857)%>%
  select(-c(1,3,4))

slab_org <- left_join(labile_pp_RA, day_transformed_transposed, by = "feature_name")%>%
  select(-c(1:5))%>%
  group_by(cluster, Organism)%>%
  summarise_at(vars(RA), mean)%>%
  ungroup()%>%
  spread(cluster, RA)

# Sifnificantly produced --------------------------------------------------
sig_exudate_names <- exudates_network_features

sig_exudate_names$feature_name <- paste0(sig_exudate_names$feature_name, ".RA")

#Networks 42, 90, 623, 596
network_42_exuded <- as.vector(sig_exudate_names%>%
                                  filter(cluster == "42"))$feature_name
network_90_exuded <- as.vector(sig_exudate_names%>%
                                  filter(cluster == "90"))$feature_name
network_623_exuded <- as.vector(sig_exudate_names%>%
                                  filter(cluster == "623"))$feature_name

network_596_exuded <- as.vector(sig_exudate_names%>%
                                 filter(cluster == "596"))$feature_name
network_7_exuded <- as.vector(sig_exudate_names%>%
                                  filter(cluster == "7"))$feature_name
network_253_exuded <- as.vector(sig_exudate_names%>%
                                filter(cluster == "253"))$feature_name

## PP production
dorcierr_production_pp_42 <- dorcierr_exudates_day%>%
  select(1:5, c(network_42_exuded))

dorcierr_production_pp_90 <- dorcierr_exudates_day%>%
  select(1:5, c(network_90_exuded))

dorcierr_production_pp_623 <- dorcierr_exudates_day%>%
  select(1:5, c(network_623_exuded))

## Calcifiers
dorcierr_production_calfiers_596 <- dorcierr_exudates_day%>%
  select(1:5, c(network_596_exuded))

dorcierr_production_calfiers_253 <- dorcierr_exudates_day%>%
  select(1:5, c(network_253_exuded))

# Writing CSVs ------------------------------------------------------------

## working data frames
write_csv(moorea_wdf, "moorea_metab_FeatureID_samples.csv")
write_csv(both_runs, "Moorea_2017_new_and_old_runs.csv")
write_delim(scan_IDs, "Feature_ID.txt", delim = " ", col_names = FALSE)
write_csv(both_runs_by_old, "new_by_old_AllMeta.csv")

write_csv(dorcierr_wdf, "Dorcierr_metab_wdf.csv")
write_csv(spiffy_wdf, "spiffy_wdf.csv") 
write_csv(dorcierr_exudates_day, "Dorcierr_exudates_day.csv")

# Change in network IDs
write_csv(Node_check, "Moorea_network_changes.csv")

# Significant compounds
write_csv(spiffy_two_locations_pvalues_max, "spiffy_significance_byreeflocation.csv")

write_csv(overlap_primary_producers, "overlap_between_dorcierr_primaryprods_spiffy.csv")
write_csv(super_labile_primary_producer_candidates, "super_labile_primary_producer_candidates_dorcierr.csv")

## Super Labile compounds for bar charts TF-T0
# write_csv(dorcierr_labile_pp_42, "dorcierr_labile_pp_42.csv")
# write_csv(dorcierr_labile_pp_90, "dorcierr_labile_pp_90.csv")
# write_csv(dorcierr_labile_pp_623, "dorcierr_labile_pp_623.csv")
# write_csv(dorcierr_labile_calcifiers_596, "dorcierr_labile_calcifiers_596.csv")

write_csv(slab_org, "super_labile_for_graphing.csv")

## Super labile compound production (T0)
write_csv(dorcierr_production_pp_42, "dorcierr_production_pp_42.csv")
write_csv(dorcierr_production_pp_90, "dorcierr_production_pp_90.csv")
write_csv(dorcierr_production_pp_623, "dorcierr_production_pp_623.csv")

write_csv(dorcierr_production_calfiers_596, "dorcierr_production_calcifiers_596.csv")
write_csv(dorcierr_production_calfiers_7, "dorcierr_production_calcifiers_7.csv")
write_csv(dorcierr_production_calfiers_253, "dorcierr_production_calfiers_253.csv")


