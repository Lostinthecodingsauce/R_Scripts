## Moorea2017Metabalomics removing blanks
## Written March 11 2019 for analyzing Mo'orea metabalomic data

# Loading libraries -------------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)

#PCoA, PERMANOVA
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)


# Reading in CSVs ---------------------------------------------------------
# True hits and analog hits are exported CSVs from GNPS
# True hits are more strictly matched to the library
true_hits <- read_tsv("Library-hits.tsv")%>%
  rename("feature_number" = '#Scan#')

analog_hits <- read_tsv("Analog-hits.tsv")%>%
  rename("feature_number" = '#Scan#')

# Node info includes networking information about each feature
node_info <- read_tsv("Node_info.tsv")%>%
  rename('feature_number' = 'cluster index',
         'network' = 'componentindex')

# Canopus tries to classify each feature
canopus_anotations <- read_csv("SIRIUS_etc/converted/Canopus_classes.csv")

# Sirius and Zodiac both try to assign molecular formulas to all the features
sirius_zodiac_anotations <- read_csv("SIRIUS_etc/converted/SIRIUS_Zodiac_converted.csv")%>%
  rename(feature_number = 1)%>%
  dplyr::select(-c(14:ncol(.)))

# Feature table has all features found within the experiments and blanks
# The columns need to be changed to the actual experiment sample codes
# Feature_table_raw is the raw export from MZMine
feature_table_raw <- read_csv("Morrea_Feayures-Table_all_Gap-Filled5.csv")%>%
  rename('feature_number' = 'row ID')

ms_sample_codes <- read_csv("Mo'orea 2017 Mass spec sample codes - Sheet1.csv")%>%
  rename('run_code' = 'Sample ID',
         'sample_code' = 'Sample Name')

## highest probability canopus annotation
# canopus gives many different classifications and the percent chance that the feature falls into that category
# This section pulls the feature whcih has the highest possiblity to be it.
canopus_best_match <- as.data.frame(canopus_anotations$name)%>%
  rename('feature_number' = 1)
canopus_best_match$canopus_probability <- apply(canopus_anotations[2:ncol(canopus_anotations)], 1, max)
canopus_best_match$canopus_annotation <- colnames(canopus_anotations)[max.col(canopus_anotations[2:1288],
                                                                              ties.method="first")]
# Combines canopus, sirus, and zodiac
super_computer_annotations <- full_join(canopus_best_match, sirius_zodiac_anotations, by = "feature_number")

## join library hits, analog hits and super computer predictions
metadata <- full_join(node_info, 
                      full_join(super_computer_annotations,
                                full_join(feature_table_raw[1:4],
                                          full_join(true_hits, analog_hits, by = "feature_number",
                                                    suffix = c("Library_", "Analog_")),
                                          by = "feature_number"),
                                by = "feature_number"),
                      by = "feature_number")

metadata$feature_number <- as.character(metadata$feature_number)

## making feature table so we can remove blanks
# have to change the MS codes for sample codes
# dplyr::selecting for Dorcierr, Mordor, RR3, Spiffy
feature_table_temp <- feature_table_raw%>%
  dplyr::select(-X326)%>%
  dplyr::select(-c(2:4))%>%
  gather(run_code, ion_charge, 2:ncol(.))%>%
  spread(feature_number, ion_charge)

feature_table_temp$run_code <- feature_table_temp$run_code%>%
  gsub(".mzXML Peak area", "", .)%>%
  gsub("_MSMS", "", .)

feature_table_dirty <- left_join(ms_sample_codes, feature_table_temp, by = "run_code")%>%
  dplyr::select(-run_code)%>%
  filter(!sample_code %like any% c("%_XAD", "%_C18", "SR%", "LoRDI%", "SE%"))%>%
  gather(feature_number, ion_charge, 2:20743)%>%
  spread(sample_code, ion_charge)
  
## defining what columns the samples are in
ions_samples <- 10:259

## defining different blanks
ions_blanks <- c(2:8, 260)


## flagging background features and subtraction from samples ------------------------------------------------
# The idea here is to flag and remove features where max(blanks)*.5 >= mean(samples)
ions <- feature_table_dirty

ions$max_blanks <- apply(ions[ions_blanks], 1, max)

ions$mean_samples <- apply(ions[ions_samples], 1, mean)


# This section finds the features where mean area under the curve across all samples is larger than .5 * max of the blanks
# The non background features are saved into a vector so that they can be filtered from the master database
# mean(samples)*0.5 > max_blanks
no_background <-ions%>%
  filter(mean_samples, mean_samples*0.5 > max_blanks)

feature_table_no_background <- left_join(no_background[1], feature_table_dirty, by = "feature_number")

## Removing transient features ---------------------------------------------
# Transient features are defined as features who's area under the curve is not more than 2E5 in at least 3 samples
# This was determined by comparing gap filled and non-gap filled data and dplyr::selecting for the lowest peak area in non-gap filled data
# This gives the assumption that anything lower than 2E5 is noise. See supplemental figure
feature_table_no_background <- cbind(feature_table_no_background,
                                  trans_feature_finder = 
                                    rowSums(feature_table_no_background[ions_samples] > 2E5))

feature_table_no_back_trans <- feature_table_no_background%>%
  filter(trans_feature_finder, trans_feature_finder >= 3)%>%
  dplyr::select(-trans_feature_finder)

## Calculate Total Ion Charge (TIC) and relativize samples -----------------
## FIND WHERE SCAN NUMBER IS GOING SO THAT THESE CAN ALL BE COMIBINED AGAIN
feature_table_RA_temp <- feature_table_no_back_trans%>%
  gather(sample_name, peak_area, 2:260)%>%
  spread(feature_number, peak_area)%>%
  add_column(TIC = apply(.[2:13539], 1, sum))

# IonCharge/TIC = RA values for every feature peak area /sum(sample peak areas)
feature_relative_abundance <- 
  as.data.frame(
    sapply(
      feature_table_RA_temp[2:13539],
      function(x) x/feature_table_RA_temp$TIC))%>%
  add_column(sample_name = feature_table_RA_temp$sample_name, .before = 1)%>%
  gather(feature_number, RA, 2:13539)%>%
  spread(sample_name, RA)

# angular transformation of RA
feature_asin_sqrt <-
  as.data.frame(
    sapply(
      feature_relative_abundance[2:260],
      function(x) asin(sqrt(x))))%>%
  add_column(feature_number = feature_table_no_back_trans$feature_number)

# All three joined together (Peak area, RA, asin(sqrt))
# Node and network info = [1:25], CANOPUS = [25:27], SIRIUS/ZODIAC = [28:40], Library Hits = [41:67], analog hits = []
# Blanks.RA = [112:119, 370], Area under the curve = [120:369], 
# RA (blanks included) = [371:630], asin(sqrt) (blanks included) = [631:888]
feature_table_combined <- 
  right_join(metadata,
               left_join(
                 left_join(
                   feature_table_no_background, feature_relative_abundance,
                   by = "feature_number", suffix = c("" , ".RA")),
                 feature_asin_sqrt, by = "feature_number", suffix = c("", ".asin(sqrt)")),
               by = "feature_number")

feature_table_wdf <- feature_table_combined%>%
  dplyr::select(feature_number, everything())

write_csv(feature_table_wdf, "moorea_feature_table_master_post_filtered.csv")
# Making Mo’orea working data frame for stats -----------------------------
moorea_transposed <- feature_table_wdf%>%
  dplyr::select(feature_number, `Blank_Lot_6350565_01.asin(sqrt)`:ncol(.))%>%
  gather(sample_ID, angular, 2:260)%>%
  spread(feature_number, angular)


moorea_transposed$sample_ID <- moorea_transposed$sample_ID%>%
  gsub(".asin\\(sqrt)", "", .)%>%
  gsub("-", "_", .)

moorea_wdf <- moorea_transposed%>%
  mutate(sample_ID = case_when(sample_ID == "SPIFFy_7_B 2" ~ "SPIFFy_6_B",
                               TRUE ~ as.character(sample_ID)))%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "Blank")%>%
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

# moorea_wdf["Timepoint"][is.na(moorea_wdf["Timepoint"])] <- 0

dorcierr <- moorea_wdf%>%
  filter(Experiment == "dorcierr")

mordor <- moorea_wdf%>%
  filter(Experiment == "mordor")

spiffy <- moorea_wdf%>%
  filter(Experiment == "SPIFFy")%>%
  dplyr::dplyr::select("Timepoint")


# Spiffy Workshopping -----------------------------------------------------

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
  dplyr::select(c(1:4, spiffy_sig_names))%>%
  group_by(reef_area)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(feature_name, RA, 2:332)%>%
  spread(reef_area, RA)

spiffy_means_only <-spiffy_only_sigs%>%
  dplyr::select(-1)

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
  dplyr::select(Site, reef_area, 5:8858)%>%
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
  dplyr::select(-max)%>%
  gather(site_name, RA, 2:16)%>%
  separate(site_name, c("site", "reef_area"), sep = "_")%>%
  dplyr::select(-site)%>%
  group_by(feature_name, reef_area)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(reef_area, RA)


# dorcierr cleanup and subsetting -----------------------------------------------
dorcierr_wdf <- dorcierr%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = -1)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))
# If low on memory. run the lines below: 
# write_csv(dorcierr_wdf, "dorcier_wdf_lowmemory.csv")
# rm(list = ls())
# dorcierr_wdf <- read_csv("dorcier_wdf_lowmemory.csv")

### Dorcierr data frame subsets:
dorcierr_final <- dorcierr_wdf%>%
  filter(Timepoint == "TF")

dorcierr_exudates_day <- dorcierr_wdf%>%
  filter(DayNight == "Day",
         Timepoint == "T0",
         !Replicate == 3,
         !Replicate == 4)

## Need to take the mean of T0 to subtract from TF because of the differences in replicate
## This section needs to be redone because it is outdated
dorcierr_transformations <- dorcierr_wdf%>%
  dplyr::select(c(1:16956))%>% 
  gather(feature_name, RA, 5:16956)%>%
  spread(Timepoint, RA)

dorcierr_t0 <- dorcierr_transformations%>%
  dplyr::select(-c(TF))%>%
  filter(!Replicate == 3,
         !Replicate == 4)%>%
  group_by(Organism, DayNight, feature_name)%>%
  summarize_if(is.numeric, mean)

dorcierr_newt0 <- left_join(dorcierr_transformations, dorcierr_t0, by = c("Organism", "DayNight", "feature_name"))%>%
  dplyr::select(-c(T0.x))%>%
  rename(T0 = 'T0.y')

dorcierr_newt0$change <- (dorcierr_newt0$TF - dorcierr_newt0$T0)

dorcierr_transformed <- dorcierr_newt0%>%
  dplyr::select(-c(T0,TF))%>%
  spread(feature_name, change)

dorcierr_day_transformed <- dorcierr_transformed%>%
  filter(DayNight == "Day")


# Dorcierr Day Exudates (T0) ----------------------------------------------
# average organism exudate dataframe
average_RA <- dorcierr_exudates_day%>%
  gather(feature_number, 6:ncol(.))

# Oneway anova significant features
p_values_oneway_day_exudates_dorc <- 
  as.data.frame(
    sapply(
      dorcierr_exudates_day[6:ncol(dorcierr_exudates_day)],
      function(x) summary(aov(x ~ dorcierr_exudates_day[["Organism"]]))[[1]][1,'Pr(>F)']))

sig_one_way_day_exudates_dorc <- p_values_oneway_day_exudates_dorc%>%
  rownames_to_column(var = "feature_number")%>%
  rename('p_value' = 2)%>%
  filter(p_value, p_value < 0.05)

sig_features_day_exudates_dorc <- as.factor(sig_one_way_day_exudates_dorc$feature_number)

# Dunnetts test testing ---------------------------------------------------
## glht will run but gives confusing output. Trying to mull through it
# Comes from Multcomp which gives problems with select function.
set.seed(20140123)
dorcierr_exudates_day_dunnetts <- dorcierr_exudates_day%>%
  unite(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_name, asin, 5:13542)%>%
  spread(org_rep, asin)%>%
  add_column(sum = apply(.[5:16], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(org_rep, asin, 5:16)%>%
  spread(feature_name, asin)%>%
  separate(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)

dunnett_model_exday <- lapply(dorcierr_exudates_day_dunnetts[6:11758], function(y, f)
  summary(glht(aov(y ~ f, dorcierr_exudates_day_dunnetts), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(dorcierr_exudates_day$Organism), "Water control")))

p_values_dunnett_exudates_day <- as.data.frame(
      sapply(
        dunnett_model_exday,
        function(x) x$test$pvalues))
  
  
p_values_dunnett_exudates_day <- cbind("Organism" = c("CCA", "Dictyota", "Pocillopora verrucosa", "Porites Lobata", "Turf"), 
                                       p_values_dunnett_exudates_day)

dunnett_FDR_exudates_day <- p_values_dunnett_exudates_day%>%
  dplyr::select(c(sig_features_day_exudates_dorc))%>%
  gather(feature_number, p_value, 2:ncol(.))

dunnett_FDR_exudates_day$FDR_f <-p.adjust(dunnett_FDR_exudates_day$p_value, method = "BH")

dunnett_sig_exudates_day <- dunnett_FDR_exudates_day%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)


  
  