## Moorea2017Metabalomics removing blanks, running stats, etc.
## Written March 11th 2019 - March 27th 2019 for analyzing Mo'orea metabalomic data
## Editted specifically for RR3 September 23rd 2019 for re-anlysis and publication

# LOADING -- libraries -------------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)
library(CHNOSZ)

#PCoA, PERMANOVA
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)


# LOADING -- dataframes ---------------------------------------------------------
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

chemont_anotations <- read_csv("categories.canopus.strings.nelsonMarch2019.CSV")%>%
  rename('canopus_annotation' = 'name')

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

# PRE-CLEANING -- Canopus---------------------------

canopus_annotation_names <- canopus_anotations%>%
  gather(canopus_annotation, canopus_probability, 2:ncol(.))

canopus_chemonnt_tidy <- left_join(canopus_annotation_names, chemont_anotations, by = "canopus_annotation")

# SET -- CANOPUS filters --------------------------------------------------
# Canopus annotations which are between levels 4-6 AND above 80% probability

canopus_filtered_tidy <- canopus_chemonnt_tidy%>%
  rename('feature_number' = 'name')%>%
  group_by(feature_number)%>%
  do(filter(., canopus_probability >= 0.80))%>%
  do(filter(., level > 3))%>%
  do(filter(., level <7))%>%
  # do(filter(., level == max(level)))%>%
  do(filter(., canopus_probability == max(canopus_probability)))%>%
  do(filter(., level == max(level)))%>%
  do(filter(., nchar(CLASS_STRING) == max(nchar(CLASS_STRING))))%>%
  do(filter(., nchar(canopus_annotation) == max(nchar(canopus_annotation))))%>%
  ungroup()

write_csv(canopus_filtered_tidy, "~/Documents/SDSU/Moorea_2017/190312_new_fusion/canopus_filtered_tidy.csv")

# SET -- Zodiac Score -----------------------------------------------------
# Combines canopus, sirus, and zodiac
super_computer_annotations <- full_join(canopus_filtered_tidy, sirius_zodiac_anotations, by = "feature_number")%>%
  filter(ZodiacScore >= 0.98,
         quality == "Good")

# PRE-CLEANING -- metadata and raw faeture table --------------------------
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
# dplyr::selecting for RR3 
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
  filter(sample_code %like any% c("%Blank%","R_%"))%>%
  gather(feature_number, ion_charge, 2:ncol(.))%>%
  spread(sample_code, ion_charge)

# DEFINING -- Samples and blanks ------------------------------------------
## defining what columns the samples are in
ions_samples <- 14:56

## defining different blanks
ions_blanks <- c(2:8, 10)


# BACKGROUND FEATURES -- flagging and removing ------------------------------------------------
# Background features are defined as features where max(blanks) >= 0.5*mean(samples)
ions <- feature_table_dirty

ions$max_blanks <- apply(ions[ions_blanks], 1, max)

ions$mean_samples <- apply(ions[ions_samples], 1, mean)


# This section finds the features where mean area under the curve across all samples is larger than 2 * max of the blanks
# The non background features are saved into a vector so that they can be filtered from the master database
# mean(samples)*0.5 > max_blanks
no_background <-ions%>%
  filter(mean_samples, mean_samples*0.5 > max_blanks)

feature_table_no_background <- left_join(no_background[1], feature_table_dirty, by = "feature_number")

# TRANSIENT FEATURES -- flagging and removing ---------------------------------------------
# Transient features are defined as features who's area under the curve is not more than 1E5 in at least 3 samples
# This was determined by comparing gap filled and non-gap filled data and selecting for the lowest peak area in non-gap filled data
# This gives the assumption that anything lower than 1E5 is noise. See supplemental figure
feature_table_no_background <- cbind(feature_table_no_background,
                                     trans_feature_finder = 
                                       rowSums(feature_table_no_background[ions_samples] > 1E5))

feature_table_no_back_trans <- feature_table_no_background%>%
  filter(trans_feature_finder, trans_feature_finder >= 3)%>%
  dplyr::select(-trans_feature_finder)


# CLASSIFYING -- Ambient or Exudate features ------------------------------
# Calculating exudates as mean(TF)/Mean(Water T0) > 1
# Replacing all 0's with 1000.
ambient_exudate_features_temp <- feature_table_no_back_trans%>%
  gather(sample_name, peak_area, 2:ncol(.))%>%
  spread(feature_number, peak_area)%>%
  filter(!sample_name %like% "%Blank%")%>%
  separate(sample_name, c("Experiment", "Organism", "Replicate", "Time"), sep = "_")%>%
  separate(Time, c("Timepoint", "DayNight"), sep = 2)%>%
  unite(sample, c("Organism", "Timepoint", "DayNight"), sep = "_")%>%
  group_by(sample)%>%
  summarize_if(is.numeric, mean)%>%
  gather(feature_number, peak_area, 2:ncol(.))%>%
  separate(sample, c("Organism", "Timepoint", "DayNight"), sep = "_")%>%
  unite(sample, c("Organism", "Timepoint"), sep = "_")%>%
  spread(sample, peak_area)

ambient_exudate_features_temp[ambient_exudate_features_temp == 0] <- 1000

ambient_exudate_features <- as.data.frame(
  sapply(
    ambient_exudate_features_temp[3:ncol(ambient_exudate_features_temp)],
    function(x) log2(x/ambient_exudate_features_temp$WA_T0)))%>%
  add_column(feature_number = ambient_exudate_features_temp$feature_number, .before = 1)%>%
  add_column(DayNight = ambient_exudate_features_temp$DayNight, .before = 1)%>%
  add_column(max = apply(.[3:9],1, max))

exudates_day <- as.vector(ambient_exudate_features%>%
                            filter(max > 1,
                                   DayNight == "D"))$feature_number

exudates_night <- as.vector(ambient_exudate_features%>%
                              filter(max > 1,
                                     DayNight == "N"))$feature_number

# RELATIVIZING DATA -- Total Ion Charge (TIC) and relativize samples -----------------
feature_table_RA_temp <- feature_table_no_back_trans%>%
  gather(sample_name, peak_area, 2:ncol(.))%>%
  spread(feature_number, peak_area)%>%
  add_column(TIC = apply(.[2:ncol(.)], 1, sum))

feature_table_log10_temp <- feature_table_no_back_trans%>%
  gather(sample_name, peak_area, 2:ncol(.))%>%
  spread(feature_number, peak_area)

feature_log10 <- as.data.frame(
  sapply(
    feature_table_log10_temp[2:ncol(feature_table_log10_temp)],
    function(x) log10(x)))%>%
  add_column(sample_name = feature_table_RA_temp$sample_name, .before = 1)%>%
  gather(feature_number, RA, 2:ncol(.))%>%
  spread(sample_name, RA)

# IonCharge/TIC = RA values for every feature peak area /sum(sample peak areas)
feature_relative_abundance <- 
  as.data.frame(
    sapply(
      feature_table_RA_temp[2:ncol(feature_table_RA_temp)],
      function(x) x/feature_table_RA_temp$TIC))%>%
  add_column(sample_name = feature_table_RA_temp$sample_name, .before = 1)%>%
  gather(feature_number, RA, 2:ncol(feature_table_RA_temp))%>%
  spread(sample_name, RA)

# angular transformation of RA
feature_asin_sqrt <-
  as.data.frame(
    sapply(
      feature_relative_abundance[2:ncol(feature_relative_abundance)],
      function(x) asin(sqrt(x))))%>%
  add_column(feature_number = feature_relative_abundance$feature_number)

# Build feature table working data frame ----------------------------------
# All three joined together (Peak area, RA, asin(sqrt))
feature_table_combined <- left_join(
  left_join(
    right_join(metadata, feature_table_no_back_trans, by = 'feature_number'),
    feature_relative_abundance, by = 'feature_number', suffix = c("", ".RA")),
  feature_asin_sqrt, by = "feature_number", suffix = c("", ".asin(sqrt)"))

feature_table_wdf <- feature_table_combined%>%
  dplyr::select(feature_number, everything())

write_csv(feature_table_wdf, "RR3_feature_table_master_post_filtered.csv")
# Making Mo’orea working data frame for stats -----------------------------
moorea_transposed <- feature_table_wdf%>%
  dplyr::select(feature_number, 226:ncol(.))%>%
  gather(sample_ID, angular, 2:ncol(.))%>%
  spread(feature_number, angular)

moorea_transposed$sample_ID <- moorea_transposed$sample_ID%>%
  gsub(".asin\\(sqrt)", "", .)%>%
  gsub("-", "_", .)

blanks_wdf <- moorea_transposed%>%
  filter(sample_ID %like% "%Blank%",
         !sample_ID == 'D_Blank_DI')%>%
  mutate(sample_ID = case_when(sample_ID == "Blank_Lot_6350565_01"~ "Blank_Blank_635056501",
                               sample_ID == "Blank_SD_01_A" ~ "Blank_Blank_SD01A",
                               sample_ID == "Blank_SD_01_B" ~ "Blank_Blank_SD01B",
                               sample_ID == "Blank_SD_LoRDI" ~ "Blank_Blank_SDLoRDI",
                               sample_ID == "Blank_SD_PPL" ~ "Blank_Blank_SDPPL",
                               sample_ID == "Blank? look up name on PPL" ~ "Blank_Blank_unknown",
                               sample_ID == "D_Blank" ~ "Blank_Blank_D",
                               TRUE ~ "Blank_Blank_Spiffy"))%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")

moorea_wdf <- moorea_transposed%>%
  mutate(sample_ID = case_when(sample_ID == "SPIFFy_7_B 2" ~ "SPIFFy_6_B",
                               TRUE ~ as.character(sample_ID)))%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank")%>%
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