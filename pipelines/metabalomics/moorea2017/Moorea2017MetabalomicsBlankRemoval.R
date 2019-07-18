## Moorea2017Metabalomics removing blanks, running stats, etc.
## Written March 11th 2019 - March 27th 2019 for analyzing Mo'orea metabalomic data

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


# Reading in Dataframes ---------------------------------------------------------
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

# Combining metadata and creating feature_table ---------------------------
## highest probability canopus annotation
# canopus gives many different classifications and the percent chance that the feature falls into that category
# This section pulls the feature whcih has the highest possiblity to be it.
# We want only canopus annotations which are level 5 AND above 70% probability
# canopus_annotation_names <- canopus_anotations%>%
#   gather(canopus_annotation, canopus_probability, 2:1288)
# 
# canopus_chemonnt_tidy <- left_join(canopus_annotation_names, chemont_anotations, by = "canopus_annotation")
# 
# canopus_filtered_tidy <- canopus_chemonnt_tidy%>%
#   rename('feature_number' = 'name')%>%
#   group_by(feature_number)%>%
#   do(filter(., canopus_probability >= 0.95))%>%
#   do(filter(., level == max(level)))%>%
#   do(filter(., canopus_probability == max(canopus_probability)))%>%
#   do(filter(., nchar(CLASS_STRING) == max(nchar(CLASS_STRING))))%>%
#   do(filter(., nchar(canopus_annotation) == max(nchar(canopus_annotation))))%>%
#   ungroup()
## Have to filter out the features which are not in the main spreadsheet
canopus_anotations_only <- canopus_anotations%>%
  gather(canopus_annotation, prob, 2:ncol(.))%>%
  spread(name, prob)

## Have to order the level in Descending Order
chemont_ordered <- chemont_anotations%>%
  filter(canopus_annotation %like any% canopus_annotations_only)
  arrange(-level)

ordered_annotations <- as.vector(chemont_ordered$canopus_annotation)

sorted <- left_join(canopus_anotations_only, chemont_ordered[1], by = "canopus_annotation")%>%
  gather(name, probability, 2:ncol(.))%>%
  spread(canopus_annotation, probability)


canopus_max <- sorted%>%
  dplyr::select(-1)

canopus_high_percentage <- canopus_anotations%>%
  add_column(canopus_annotation = colnames(canopus_max)[max.col(canopus_max, ties.method = "first")], .before=2)%>%
  add_column(probability = apply(canopus_max, 1, max), .after = 2)%>%
  dplyr::select(1:3)


canopus_chemont_high_percentage <- left_join(canopus_high_percentage, chemont_anotations, by = "canopus_annotation")%>%
  rename(feature_number = 1)

# Combines canopus, sirus, and zodiac
super_computer_annotations <- full_join(canopus_chemont_high_percentage, sirius_zodiac_anotations, by = "feature_number")

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

# Build feature table working data frame ----------------------------------
# All three joined together (Peak area, RA, asin(sqrt))
# Node and network info = [1:25], CANOPUS = [25:27], SIRIUS/ZODIAC = [28:40], Library Hits = [41:67], analog hits = []
# Blanks.RA = [112:119, 370], Area under the curve = [120:369], 
# RA (blanks included) = [371:630], asin(sqrt) (blanks included) = [631:888]
feature_table_combined <- left_join(
  left_join(
    right_join(metadata, feature_table_no_back_trans, by = 'feature_number'),
    feature_relative_abundance, by = 'feature_number', suffix = c("", ".RA")),
  feature_asin_sqrt, by = "feature_number", suffix = c("", ".asin(sqrt)"))

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


dorcierr <- moorea_wdf%>%
  filter(Experiment == "dorcierr")

mordor <- moorea_wdf%>%
  filter(Experiment == "mordor")

spiffy <- moorea_wdf%>%
  filter(Experiment == "SPIFFy")%>%
  dplyr::select(-"Timepoint")

rr3 <- moorea_wdf%>%
  filter(Experiment == "RR3")


# Making RA dataframes ----------------------------------------------------
moorea_transposed_RA <- feature_table_wdf%>%
  dplyr::select(feature_number, `Blank_Lot_6350565_01.RA`:`SPIFFy_Blank.RA`)%>%
  gather(sample_ID, angular, 2:260)%>%
  spread(feature_number, angular)

moorea_transposed_RA$sample_ID <- moorea_transposed_RA$sample_ID%>%
  gsub(".RA", "", .)%>%
  gsub("-", "_", .)

moorea_wdf_RA <- moorea_transposed_RA%>%
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

dorcierr_RA <- moorea_wdf_RA%>%
  filter(Experiment == "dorcierr")

dorcierr_wdf_RA <- dorcierr_RA%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = -1)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))

dorcierr_exudates_day_RA <- dorcierr_wdf_RA%>%
  filter(DayNight == "Day",
         Timepoint == "T0",
         !Replicate == 3,
         !Replicate == 4)

dorcierr_exudates_night_RA <- dorcierr_wdf_RA%>%
  filter(DayNight == "Night",
         Timepoint == "T0",
         !Replicate == 3,
         !Replicate == 4)

rr3_day_RA <- moorea_wdf_RA%>%
  filter(Experiment == "dorcierr")%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = -1)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  filter(Timepoint == 'TF',
         DayNight == "Day")

rr3_night_RA <- moorea_wdf_RA%>%
  filter(Experiment == "dorcierr")%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = -1)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  filter(Timepoint == 'TF',
         DayNight == "Night")


# Spiffy Cleaning and Subsetting-----------------------------------------------------
spiffy_almost_wdf <- spiffy%>%
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

spiffy_rep_noise <- spiffy_almost_wdf%>%
  gather(feature_number, asin, 5:ncol(.))%>%
  add_column("binary" = 1)%>%
  group_by(Site, feature_number)%>%
  filter(., asin != 0.00)%>%
  summarize_if(is.numeric, sum)%>%
  filter(binary == 2)%>%
  ungroup()

spiffY_real_features <- as.vector(spiffy_rep_noise$feature_number)

spiffy_wdf <- spiffy_almost_wdf%>%
  dplyr::select(c(1:4, spiffY_real_features))


spiffy_no_bay <-spiffy_wdf%>%
  filter(!reef_area == "bay")

spiffy_no_bay_RA <- moorea_wdf_RA%>%
  filter(Experiment == "SPIFFy")%>%
  dplyr::select(-"Timepoint")%>%
  add_column(reef_area = .$Organism, .before = 3)%>%
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
                               TRUE ~ "Forereef"))%>%
  filter(!reef_area == "bay")



# Spiffy Oneway Anova -----------------------------------------------------
## One way anova looking at only backreef and fore reef
aov_spiffy <- as.data.frame(sapply(spiffy_no_bay[5:ncol(spiffy_no_bay)], function(x) summary(aov(x ~ spiffy_no_bay[["reef_area"]]))[[1]][1,'Pr(>F)']))

anova_spiffy <- aov_spiffy%>%
  rownames_to_column(var = "feature_name")%>%
  rename(`f_value`= `sapply(spiffy_no_bay[5:ncol(spiffy_no_bay)], function(x) summary(aov(x ~ spiffy_no_bay[["reef_area"]]))[[1]][1, "Pr(>F)"])`)

anova_spiffy$FDR_p <- p.adjust(anova_spiffy$f_value, method = "BH")

anova_sigs_spiffy <- anova_spiffy%>%
  filter(FDR_p, FDR_p < 0.05)

spiffy_sig_names <-as.vector(anova_sigs_spiffy$feature_name)

spiffy_only_sigs <- spiffy_no_bay%>%
  dplyr::select(c(1:4, spiffy_sig_names))%>%
  group_by(reef_area)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(feature_name, RA, 2:1480)%>%
  spread(reef_area, RA)

# Spiffy largest player table ---------------------------------------------
## Making spiffy max table
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

# Spiffy significant columns ----------------------------------------------
## Making spiffy columns to be added to feature_table_master
spiffy_sig_columns <-anova_sigs_spiffy%>%
  dplyr::select(-f_value)%>%
  add_column('diff_back-fore_spiffy' = spiffy_only_sigs$backreef-spiffy_only_sigs$Forereef)%>%
  add_column('%diff_spiffy' = (spiffy_only_sigs$backreef-spiffy_only_sigs$Forereef)/(spiffy_only_sigs$backreef+spiffy_only_sigs$Forereef))%>%
  rename('feature_number' = 'feature_name')

## Finding rare features in Spiffy
## Needs to be modified to find RA
spiffy_rare <- spiffy_no_bay_RA%>%
  unite(site_name, c(Site, reef_area), sep = "_", remove = TRUE)%>%
  group_by(site_name)%>%
  summarize_if(is.numeric, mean)%>%
  gather(feature_name, RA, 2:13539)%>%
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
  spread(reef_area, RA)%>%
  add_column('rare_spiffy' = "1", .before = 2)%>%
  add_column('rare_back-fore_diff' = .$backreef-.$Forereef, .after = 2)%>%
  add_column('rare_back-fore_%diff' = (.$backreef-.$Forereef)/(.$backreef+.$Forereef))%>%
  dplyr::select(-c('backreef', 'Forereef'))%>%
  rename('feature_number' = 'feature_name')



# Dorcierr cleanup and subsetting -----------------------------------------------
dorcierr_wdf <- dorcierr%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = -1)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))
# If low on memory. run the lines below: 
# write_csv(dorcierr_wdf, "dorcier_wdf_lowmemory.csv")
# rm(list = ls())
# dorcierr_wdf <- read_csv("dorcier_wdf_lowmemory.csv")

### Dorcierr data frame subsets:
dorcierr_exudates_day <- dorcierr_wdf%>%
  filter(DayNight == "Day",
         Timepoint == "T0",
         !Replicate == 3,
         !Replicate == 4)

dorcierr_remins_day <- dorcierr_wdf%>%
  filter(DayNight == "Day")

dorcierr_exudates_night <- dorcierr_wdf%>%
  filter(DayNight == "Night",
         Timepoint == "T0")

dorcierr_remins_day <- dorcierr_wdf%>%
  filter(!DayNight == "Night")
## Need to take the mean of T0 to subtract from TF because of the differences in replicate
## This section needs to be redone because it is outdated
dorcierr_transformations <- dorcierr_wdf_RA%>%
  gather(feature_number, RA, 6:13543)%>%
  spread(Timepoint, RA)

dorcierr_t0 <- dorcierr_transformations%>%
  dplyr::select(-c(TF))%>%
  filter(!Replicate == 3,
         !Replicate == 4)%>%
  group_by(Organism, DayNight, feature_number)%>%
  summarize_if(is.numeric, mean)

dorcierr_newt0 <- left_join(dorcierr_transformations, dorcierr_t0, by = c("Organism", "DayNight", "feature_number"))%>%
  dplyr::select(-c(T0.x))%>%
  rename(T0 = 'T0.y')%>%
  add_column(change = .$TF - .$T0)%>%
  add_column(percent_change = (.$TF - .$T0)/(.$TF + .$T0))

dorcierr_change <- dorcierr_newt0%>%
  dplyr::select(-c(percent_change, T0, TF, Replicate, Experiment))%>%
  group_by(feature_number, DayNight, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(Organism, change)%>%
  rename('CCA_diff_Tf-T0_dorcierr' = 'CCA',
         'Dictyota_diff_Tf-T0_dorcierr' = 'Dictyota',
         'Pocillpora verrucosa_diff_Tf-T0_dorcierr' = 'Pocillopora verrucosa',
         'Porites_lobata_diff_Tf-T0_dorcierr' = 'Porites lobata',
         'Turf_diff_Tf-T0_dorcierr' = 'Turf',
         'Water control_diff_Tf-T0_dorcierr' = 'Water control')
  
dorcierr_percent_change <- dorcierr_newt0%>%
  dplyr::select(-c(change, T0, TF, Replicate, Experiment))%>%
  group_by(feature_number, DayNight, Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(Organism, percent_change)%>%
  rename('CCA_%diff_Tf-T0_dorcierr' = 'CCA',
         'Dictyota_%diff_Tf-T0_dorcierr' = 'Dictyota',
         'Pocillpora verrucosa_%diff_Tf-T0_dorcierr' = 'Pocillopora verrucosa',
         'Porites_lobata_%diff_Tf-T0_dorcierr' = 'Porites lobata',
         'Turf_%diff_Tf-T0_dorcierr' = 'Turf',
         'Water control_%diff_Tf-T0_dorcierr' = 'Water control')


dorcierr_transformed <- left_join(dorcierr_change, dorcierr_percent_change, by = c('feature_number', 'DayNight'))

## Now we have filtered to just the Day samples and transposed so that we will have delta added into the master
dorcierr_day_transformed_RA <- dorcierr_transformed%>%
  filter(DayNight == "Day")


# Dorcierr Day Exudates (T0) ANOVA----------------------------------------------
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

# Dorcierr Day transformations (T0 to TF) -------------------------------------
## dorcierr remineralizations T0 and TF to test each species difference between T0 and TF
# This filters out anything that is 0's across the board
dorcierr_remineralizations_day_anova <- dorcierr_remins_day%>%
  unite(org_rep, c("Organism", "Replicate", "Timepoint"), sep = "_", remove = TRUE)%>%
  gather(feature_number, asin, 4:13541)%>%
  spread(org_rep, asin)%>%
  add_column(sum = apply(.[4:39], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(org_rep, asin, 4:39)%>%
  separate(org_rep, c("Organism", "Replicate", "Timepoint"), sep = "_", remove = TRUE)%>%
  unite(feat_org, c("feature_number", "Organism"), sep = "_", remove = TRUE)

#Oneway model for TF - T0
oneway_day_remins_dorc <- dorcierr_remineralizations_day_anova%>%
  group_by(feat_org)%>%
  do(p_value =  summary(aov(asin ~ Timepoint, data = .))[[1]][1,'Pr(>F)'])

#P_values
p_values_oneway_day_remins_dorc <- oneway_day_remins_dorc%>%
  separate(feat_org, c("feature_number", "Organism"), sep = "_")%>%
  filter(!p_value == 'NaN')%>%
  add_column(FDR_f = p.adjust(.$p_value, method = "BH"))%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)%>%
  rename('CCA_P_value_D_t0_tf_day_ttest' = 'CCA',
         'Dictyota_P_value_D_t0_tf_day_ttest' = 'Dictyota',
         'Pocillopora verrucosa_P_value_D_t0_tf_day_ttest' = 'Pocillopora verrucosa',
         'Porites lobata_P_value_D_t0_tf_day_ttest' = 'Porites lobata',
         'Turf_P_value_D_t0_tf_day_ttest' = 'Turf',
         'Water control_P_value_D_t0_tf_day_ttest' = 'Water control')

sig_features_timepoint_anova <-as.vector(p_values_oneway_day_remins_dorc$feature_number)

dorcierr_day_transformed_RA_sig <- dorcierr_day_transformed_RA%>%
  filter(feature_number %like any% c(sig_features_timepoint_anova))

timepoint_anova_sigs <- p_values_oneway_day_remins_dorc%>%
  add_column('CCA_diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`CCA_diff_Tf-T0_dorcierr`, .after = 2)%>%
  add_column('CCA_%diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`CCA_%diff_Tf-T0_dorcierr`, .after = 3)%>%
  add_column('Dictyota_diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`Dictyota_diff_Tf-T0_dorcierr`, .after = 5)%>%
  add_column('Dictyota_%diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`Dictyota_%diff_Tf-T0_dorcierr`, .after = 6)%>%
  add_column('Pocillpora verrucosa_diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`Pocillpora verrucosa_diff_Tf-T0_dorcierr`, .after = 8)%>%
  add_column('Pocillpora verrucosa_%diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`Pocillpora verrucosa_%diff_Tf-T0_dorcierr`, .after = 9)%>%
  add_column('Porites_lobata_diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`Porites_lobata_diff_Tf-T0_dorcierr`, .after = 11)%>%
  add_column('Porites_lobata_%diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`Porites_lobata_%diff_Tf-T0_dorcierr`, .after = 12)%>%
  add_column('Turf_diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`Turf_diff_Tf-T0_dorcierr`, .after = 14)%>%
  add_column('Turf_%diff_Tf-T0_D_day' = dorcierr_day_transformed_RA_sig$`Turf_%diff_Tf-T0_dorcierr`, .after = 15)%>%
  add_column('Water control_diff_TF-T0_D_day' = dorcierr_day_transformed_RA_sig$`Water control_diff_Tf-T0_dorcierr`, .after = 17)%>%
  add_column('Water control_%diff_TF-T0_D_day' = dorcierr_day_transformed_RA_sig$`Water control_%diff_Tf-T0_dorcierr`, .after = 18)


# Dorcierr Dunnetts Day ---------------------------------------------------
# Comes from Multcomp which gives problems with select function.
set.seed(2005)
## dorcierr exudates only looks at T0
## There are some NA's or 0s which interfere with the dunnett model
dorcierr_exudates_day_dunnetts <- dorcierr_exudates_day%>%
  dplyr::select(c(1:5, sig_features_day_exudates_dorc))%>%
  unite(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_name, asin, 5:3311)%>%
  spread(org_rep, asin)%>%
  add_column(sum = apply(.[5:16], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(org_rep, asin, 5:16)%>%
  spread(feature_name, asin)%>%
  separate(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)

dunnett_model_exday <- lapply(dorcierr_exudates_day_dunnetts[6:2851], function(y, f)
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

rm(dunnett_model_exday)  
  
p_values_dunnett_exudates_day <- cbind("Organism" = c("CCA_WA_P_value_D_exudates_day_dunnett",
                                                      "Dictyota_WA_P_value_D_exudates_day_dunnett", 
                                                      "Pocillopora verrucosa_WA_P_value_D_exudates_day_dunnett",
                                                      "Porites Lobata_WA_P_value_D_exudates_day_dunnett",
                                                      "Turf_WA_P_value_D_exudates_day_dunnett"),
                                       p_values_dunnett_exudates_day)

dunnett_FDR_exudates_day <- p_values_dunnett_exudates_day%>%
  gather(feature_number, p_value, 2:ncol(.))

dunnett_FDR_exudates_day$FDR_f <-p.adjust(dunnett_FDR_exudates_day$p_value, method = "BH")

dunnett_sig_exudates_day <- dunnett_FDR_exudates_day%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)

dunnett_sig_features <- as.vector(dunnett_sig_exudates_day$feature_number)

# average organism exudate dataframe
exudate_average_RA <- dorcierr_exudates_day_RA%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  gather(feature_number, RA, 2:ncol(.))%>%
  spread(Organism, RA)

exudate_difference_RA <- as.data.frame(
  sapply(
    exudate_average_RA[2:6],
    function(x) x-exudate_average_RA$`Water control`))%>%
  add_column(feature_number = exudate_average_RA$feature_number, .before = 1)

exudate_percent_difference_RA <- as.data.frame(
  sapply(exudate_average_RA[2:6],
         function(x) (x-exudate_average_RA$`Water control`)/(x+exudate_average_RA$`Water control`)))%>%
  add_column(feature_number = exudate_average_RA$feature_number, .before = 1)

exudate_diff_columns_RA <- left_join(exudate_difference_RA, exudate_percent_difference_RA, by = 'feature_number', suffix = c("_diff", "_%diff"))%>%
  filter(feature_number %like any% c(dunnett_sig_features))

## adding in difference and percent difference to sig pvalues
dunnett_sig_columns <- dunnett_sig_exudates_day%>%
  add_column('CCA_diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$CCA_diff, .after = 2)%>%
  add_column('CCA_%diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$`CCA_%diff`, .after = 3)%>%
  add_column('Dictyota_diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$Dictyota_diff, .after = 5)%>%
  add_column('Dictyota_%diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$`Dictyota_%diff`, .after = 6)%>%
  add_column('Pocillopora verrucosa_diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$`Pocillopora verrucosa_diff`, .after = 8)%>%
  add_column('Pocillopora verrucosa_%diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$`Pocillopora verrucosa_%diff`, .after = 9)%>%
  add_column('Porites lobata_diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$`Porites lobata_diff`, .after = 11)%>%
  add_column('Porites lobata_%diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$`Porites lobata_%diff`, .after = 12)%>%
  add_column('Turf_diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$Turf_diff, .after = 14)%>%
  add_column('Turf_%diff_D_exudates_day_dunnett' = exudate_diff_columns_RA$`Turf_%diff`, .after = 15)

# Dorcierr Night Exudates (T0) Night --------------------------------------
# Oneway anova significant features
p_values_oneway_night_exudates_dorc <- 
  as.data.frame(
    sapply(
      dorcierr_exudates_night[6:ncol(dorcierr_exudates_night)],
      function(x) summary(aov(x ~ dorcierr_exudates_night[["Organism"]]))[[1]][1,'Pr(>F)']))

sig_one_way_night_exudates_dorc <- p_values_oneway_night_exudates_dorc%>%
  rownames_to_column(var = "feature_number")%>%
  rename('p_value' = 2)%>%
  filter(p_value, p_value < 0.05)

sig_features_night_exudates_dorc <- as.factor(sig_one_way_night_exudates_dorc$feature_number)

# Dorcierr Night Dunnetts -------------------------------------------------
## dorcierr exudates only looks at T0
## There are some NA's or 0s which interfere with the dunnett model
dorcierr_exudates_night_dunnetts <- dorcierr_exudates_night%>%
  dplyr::select(c(1:5, sig_features_night_exudates_dorc))%>%
  unite(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_name, asin, 5:ncol(.))%>%
  spread(org_rep, asin)%>%
  add_column(sum = apply(.[5:16], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(org_rep, asin, 5:16)%>%
  spread(feature_name, asin)%>%
  separate(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)

dunnett_model_exnight <- lapply(dorcierr_exudates_night_dunnetts[6:ncol(dorcierr_exudates_night_dunnetts)], function(y, f)
  summary(glht(aov(y ~ f, dorcierr_exudates_night_dunnetts), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(dorcierr_exudates_night$Organism), "Water control")))

p_values_dunnett_exudates_night <- as.data.frame(
  sapply(
    dunnett_model_exnight,
    function(x) x$test$pvalues))

rm(dunnett_model_exnight)  

p_values_dunnett_exudates_night <- cbind("Organism" = c("CCA_WA_P_value_D_exudates_night_dunnett",
                                                      "Dictyota_WA_P_value_D_exudates_night_dunnett", 
                                                      "Pocillopora verrucosa_WA_P_value_D_exudates_night_dunnett",
                                                      "Porites Lobata_WA_P_value_D_exudates_night_dunnett",
                                                      "Turf_WA_P_value_D_exudates_night_dunnett"),
                                       p_values_dunnett_exudates_night)

dunnett_FDR_exudates_night <- p_values_dunnett_exudates_night%>%
  gather(feature_number, p_value, 2:ncol(.))

dunnett_FDR_exudates_night$FDR_f <-p.adjust(dunnett_FDR_exudates_night$p_value, method = "BH")

dunnett_sig_exudates_night <- dunnett_FDR_exudates_night%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)

dunnett_night_sig_features <- as.vector(dunnett_sig_exudates_night$feature_number)

# average organism exudate dataframe
exudate_night_average_RA <- dorcierr_exudates_night_RA%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  gather(feature_number, RA, 2:ncol(.))%>%
  spread(Organism, RA)

exudate_night_difference_RA <- as.data.frame(
  sapply(
    exudate_night_average_RA[2:6],
    function(x) x-exudate_night_average_RA$`Water control`))%>%
  add_column(feature_number = exudate_night_average_RA$feature_number, .before = 1)

exudate_night_percent_difference_RA <- as.data.frame(
  sapply(exudate_night_average_RA[2:6],
         function(x) (x-exudate_night_average_RA$`Water control`)/(x+exudate_night_average_RA$`Water control`)))%>%
  add_column(feature_number = exudate_night_average_RA$feature_number, .before = 1)

exudate_night_diff_columns_RA <- left_join(exudate_night_difference_RA, exudate_night_percent_difference_RA, by = 'feature_number', suffix = c("_diff", "_%diff"))%>%
  filter(feature_number %like any% c(dunnett_night_sig_features))

## adding in difference and percent difference to sig pvalues
dunnett_night_sig_columns <- dunnett_sig_exudates_night%>%
  add_column('CCA_diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$CCA_diff, .after = 2)%>%
  add_column('CCA_%diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$`CCA_%diff`, .after = 3)%>%
  add_column('Dictyota_diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$Dictyota_diff, .after = 5)%>%
  add_column('Dictyota_%diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$`Dictyota_%diff`, .after = 6)%>%
  add_column('Pocillopora verrucosa_diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$`Pocillopora verrucosa_diff`, .after = 8)%>%
  add_column('Pocillopora verrucosa_%diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$`Pocillopora verrucosa_%diff`, .after = 9)%>%
  add_column('Porites lobata_diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$`Porites lobata_diff`, .after = 11)%>%
  add_column('Porites lobata_%diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$`Porites lobata_%diff`, .after = 12)%>%
  add_column('Turf_diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$Turf_diff, .after = 14)%>%
  add_column('Turf_%diff_D_exudates_night_dunnett' = exudate_night_diff_columns_RA$`Turf_%diff`, .after = 15)

# RR3 cleaning ------------------------------------------------------------
rr3_wdf <- rr3%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = -1)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  filter(Timepoint == 'TF')


rr3_day <- rr3_wdf%>%
  filter(DayNight == "Day")

rr3_night <- rr3_wdf%>%
  filter(DayNight == "Night")
# RR3 Day Anova Organism ------------------------------------------------------
p_values_oneway_day_rr3 <- 
  as.data.frame(
    sapply(
      rr3_day[6:ncol(rr3_day)],
      function(x) summary(aov(x ~ rr3_day[["Organism"]]))[[1]][1,'Pr(>F)']))

sig_one_way_day_rr3 <- p_values_oneway_day_rr3%>%
  rownames_to_column(var = "feature_number")%>%
  rename('p_value' = 2)%>%
  filter(p_value, p_value < 0.05)

sig_features_day_rr3 <- as.factor(sig_one_way_day_rr3$feature_number)


# RR3 Day Dunnetts ------------------------------------------------------------
rr3_day_dunnetts <- rr3_day%>%
  dplyr::select(c(1:5, sig_features_day_rr3))%>%
  unite(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_name, asin, 5:ncol(.))%>%
  spread(org_rep, asin)%>%
  add_column(sum = apply(.[5:22], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(org_rep, asin, 5:22)%>%
  spread(feature_name, asin)%>%
  separate(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)

dunnett_model_rrday <- lapply(rr3_day_dunnetts[6:ncol(rr3_day_dunnetts)], function(y, f)
  summary(glht(aov(y ~ f, rr3_day_dunnetts), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(rr3_day_dunnetts$Organism), "Water control")))

p_values_dunnett_rr3_day <- as.data.frame(
  sapply(
    dunnett_model_rrday,
    function(x) x$test$pvalues))

p_values_dunnett_rr3_day <- cbind("Organism" = c("CCA_WA_P_value_rr_day_dunnett",
                                                 "Dictyota_WA_P_value_rr_day_dunnett", 
                                                 "Pocillopora verrucosa_WA_P_value_rr_day_dunnett",
                                                 "Porites Lobata_WA_P_value_rr_day_dunnett",
                                                 "Turf_WA_P_value_rr_day_dunnett"),
                                  p_values_dunnett_rr3_day)

dunnett_FDR_rr3_day <- p_values_dunnett_rr3_day%>%
  gather(feature_number, p_value, 2:ncol(.))

dunnett_FDR_rr3_day$FDR_f <-p.adjust(dunnett_FDR_rr3_day$p_value, method = "BH")

dunnett_sig_rr3_day <- dunnett_FDR_rr3_day%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)

dunnett_rr3_sig_features <- as.vector(dunnett_sig_rr3_day$feature_number)


# RR3_RA Day columns ---------------------------------------------------------
rr3_average_RA <- rr3_day_RA%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  gather(feature_number, RA, 2:ncol(.))%>%
  spread(Organism, RA)

rr3_difference_RA <- as.data.frame(
  sapply(
    rr3_average_RA[2:6],
    function(x) x-rr3_average_RA$`Water control`))%>%
  add_column(feature_number = rr3_average_RA$feature_number, .before = 1)

rr3_percent_difference_RA <- as.data.frame(
  sapply(rr3_average_RA[2:6],
         function(x) (x-rr3_average_RA$`Water control`)/(x+rr3_average_RA$`Water control`)))%>%
  add_column(feature_number = exudate_average_RA$feature_number, .before = 1)

rr3_diff_columns_RA <- left_join(rr3_difference_RA, rr3_percent_difference_RA, by = 'feature_number', suffix = c("_diff", "_%diff"))%>%
  filter(feature_number %like any% c(dunnett_rr3_sig_features))


## adding in difference and percent difference to sig pvalues
rr3_dunnett_sig_columns <- dunnett_sig_rr3_day%>%
  add_column('CCA_diff_rr3_day_dunnett' = rr3_diff_columns_RA$CCA_diff, .after = 2)%>%
  add_column('CCA_%diff_rr3_day_dunnett' = rr3_diff_columns_RA$`CCA_%diff`, .after = 3)%>%
  add_column('Dictyota_diff_rr3_day_dunnett' = rr3_diff_columns_RA$Dictyota_diff, .after = 5)%>%
  add_column('Dictyota_%diff_rr3_day_dunnett' = rr3_diff_columns_RA$`Dictyota_%diff`, .after = 6)%>%
  add_column('Pocillopora verrucosa_diff_rr3_day_dunnett' = rr3_diff_columns_RA$`Pocillopora verrucosa_diff`, .after = 8)%>%
  add_column('Pocillopora verrucosa_%diff_rr3_day_dunnett' = rr3_diff_columns_RA$`Pocillopora verrucosa_%diff`, .after = 9)%>%
  add_column('Porites lobata_diff_rr3_day_dunnett' = rr3_diff_columns_RA$`Porites lobata_diff`, .after = 11)%>%
  add_column('Porites lobata_%diff_rr3_day_dunnett' = rr3_diff_columns_RA$`Porites lobata_%diff`, .after = 12)%>%
  add_column('Turf_diff_rr3_day_dunnett' = rr3_diff_columns_RA$Turf_diff, .after = 14)%>%
  add_column('Turf_%diff_rr3_day_dunnett' = rr3_diff_columns_RA$`Turf_%diff`, .after = 15)


# RR3 Day Tukeys -------------------------------------------------------------
tukey_model_rrday <- sapply(rr3_day_dunnetts[6:ncol(rr3_day_dunnetts)], function(x)
  TukeyHSD(aov(x ~ rr3_day_dunnetts$Organism, data = rr3_day_dunnetts), p.adjust.methods = "BH"))

p_values_tukey_rr3_day <- as.data.frame(tukey_model_rrday)%>%
  rownames_to_column(var = "variable")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% '%p.adj%')%>%
  filter(value, value < 0.05)

p_values_tukey_rr3_day$feature_info <- p_values_tukey_rr3_day$feature_info%>%
  gsub("X", "", .)%>%
  gsub("rr3_day_dunnetts.Organism.p.adj", "adj_p_rr3_tukey_day", .)

sig_tukey_rr3_day <- p_values_tukey_rr3_day%>%
  separate(feature_info, c("feature_number", "test_info"), sep = "\\.")%>%
  unite(column_names, c(variable, test_info), sep = "_")

tukey_sig_figs <- as.vector(sig_tukey_rr3_day%>% spread(column_names, value))$feature_number

## This part is cheesy. I am letting Tukey calculate diff for me for RA values. P_values not pulled form this test
rr3_day_filtered_RA <- rr3_day_RA%>%
  dplyr::select(1:5, c(tukey_sig_figs))

diff_model_rrra <- sapply(rr3_day_filtered_RA[6:ncol(rr3_day_filtered_RA)], function(x)
  TukeyHSD(aov(x ~ rr3_day_filtered_RA$Organism, data = rr3_day_filtered_RA)))

diff_tukey_rr3_day <- as.data.frame(diff_model_rrra)%>%
  rownames_to_column(var = "variable")%>%
  gather(feature_info, value, 2:54153)%>%
  filter(feature_info %like% '%diff%')

diff_tukey_rr3_day$feature_info <- diff_tukey_rr3_day$feature_info%>%
  gsub("X", "", .)%>%
  gsub("rr3_day_filtered_RA.Organism.diff", "diff_rr3_tukey_day", .)

sig_diff_tukey_rr3_day <- diff_tukey_rr3_day%>%
  separate(feature_info, c("feature_number", "test_info"), sep = "\\.")%>%
  unite(column_names, c(variable, test_info), sep = "_")%>%
  filter(feature_number %like any% c(tukey_sig_figs))

## combining diff and p_value columns
rr3_tukey_sig_columns <- bind_rows(sig_tukey_rr3_day, sig_diff_tukey_rr3_day)%>%
  arrange(column_names)%>%
  spread(column_names, value)

# RR3 night Anova Organism ------------------------------------------------------
p_values_oneway_night_rr3 <- 
  as.data.frame(
    sapply(
      rr3_night[6:ncol(rr3_night)],
      function(x) summary(aov(x ~ rr3_night[["Organism"]]))[[1]][1,'Pr(>F)']))

sig_one_way_night_rr3 <- p_values_oneway_night_rr3%>%
  rownames_to_column(var = "feature_number")%>%
  rename('p_value' = 2)%>%
  filter(p_value, p_value < 0.05)

sig_features_night_rr3 <- as.factor(sig_one_way_night_rr3$feature_number)


# RR3 night Dunnetts ------------------------------------------------------------
rr3_night_dunnetts <- rr3_night%>%
  dplyr::select(c(1:5, sig_features_night_rr3))%>%
  unite(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_name, asin, 5:ncol(.))%>%
  spread(org_rep, asin)%>%
  add_column(sum = apply(.[5:22], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(org_rep, asin, 5:22)%>%
  spread(feature_name, asin)%>%
  separate(org_rep, c("Organism", "Replicate"), sep = "_", remove = TRUE)

dunnett_model_rrnight <- lapply(rr3_night_dunnetts[6:ncol(rr3_night_dunnetts)], function(y, f)
  summary(glht(aov(y ~ f, rr3_night_dunnetts), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(rr3_night_dunnetts$Organism), "Water control")))

p_values_dunnett_rr3_night <- as.data.frame(
  sapply(
    dunnett_model_rrnight,
    function(x) x$test$pvalues))

p_values_dunnett_rr3_night <- cbind("Organism" = c("CCA_WA_P_value_rr_night_dunnett",
                                                 "Dictyota_WA_P_value_rr_night_dunnett", 
                                                 "Pocillopora verrucosa_WA_P_value_rr_night_dunnett",
                                                 "Porites Lobata_WA_P_value_rr_night_dunnett",
                                                 "Turf_WA_P_value_rr_night_dunnett"),
                                  p_values_dunnett_rr3_night)

dunnett_FDR_rr3_night <- p_values_dunnett_rr3_night%>%
  gather(feature_number, p_value, 2:ncol(.))

dunnett_FDR_rr3_night$FDR_f <-p.adjust(dunnett_FDR_rr3_night$p_value, method = "BH")

dunnett_sig_rr3_night <- dunnett_FDR_rr3_night%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)

dunnett_rr3_night_sig_features <- as.vector(dunnett_sig_rr3_night$feature_number)


# RR3_RA night columns ---------------------------------------------------------
rr3_average_night_RA <- rr3_night_RA%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  gather(feature_number, RA, 2:ncol(.))%>%
  spread(Organism, RA)

rr3_difference_night_RA <- as.data.frame(
  sapply(
    rr3_average_night_RA[2:6],
    function(x) x-rr3_average_night_RA$`Water control`))%>%
  add_column(feature_number = rr3_average_night_RA$feature_number, .before = 1)

rr3_percent_difference_night_RA <- as.data.frame(
  sapply(rr3_average_night_RA[2:6],
         function(x) (x-rr3_average_night_RA$`Water control`)/(x+rr3_average_night_RA$`Water control`)))%>%
  add_column(feature_number = rr3_average_night_RA$feature_number, .before = 1)

rr3_night_diff_columns_RA <- left_join(rr3_difference_night_RA, rr3_percent_difference_night_RA, by = 'feature_number', suffix = c("_diff", "_%diff"))%>%
  filter(feature_number %like any% c(dunnett_rr3_night_sig_features))


## adding in difference and percent difference to sig pvalues
rr3_night_dunnett_sig_columns <- dunnett_sig_rr3_night%>%
  add_column('CCA_diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$CCA_diff, .after = 2)%>%
  add_column('CCA_%diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$`CCA_%diff`, .after = 3)%>%
  add_column('Dictyota_diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$Dictyota_diff, .after = 5)%>%
  add_column('Dictyota_%diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$`Dictyota_%diff`, .after = 6)%>%
  add_column('Pocillopora verrucosa_diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$`Pocillopora verrucosa_diff`, .after = 8)%>%
  add_column('Pocillopora verrucosa_%diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$`Pocillopora verrucosa_%diff`, .after = 9)%>%
  add_column('Porites lobata_diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$`Porites lobata_diff`, .after = 11)%>%
  add_column('Porites lobata_%diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$`Porites lobata_%diff`, .after = 12)%>%
  add_column('Turf_diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$Turf_diff, .after = 14)%>%
  add_column('Turf_%diff_rr3_night_dunnett' = rr3_night_diff_columns_RA$`Turf_%diff`, .after = 15)


# RR3 Night Tukeys -------------------------------------------------------------
tukey_model_rrnight <- sapply(rr3_night_dunnetts[6:ncol(rr3_night_dunnetts)], function(x)
  TukeyHSD(aov(x ~ rr3_night_dunnetts$Organism, data = rr3_night_dunnetts), p.adjust.methods = "BH"))

p_values_tukey_rr3_night <- as.data.frame(tukey_model_rrnight)%>%
  rownames_to_column(var = "variable")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% '%p.adj%')%>%
  filter(value, value < 0.05)

p_values_tukey_rr3_night$feature_info <- p_values_tukey_rr3_night$feature_info%>%
  gsub("X", "", .)%>%
  gsub("rr3_night_dunnetts.Organism.p.adj", "adj_p_rr3_tukey_night", .)

sig_tukey_rr3_night <- p_values_tukey_rr3_night%>%
  separate(feature_info, c("feature_number", "test_info"), sep = "\\.")%>%
  unite(column_names, c(variable, test_info), sep = "_")

tukey_night_sig_figs <- as.vector(sig_tukey_rr3_night%>% spread(column_names, value))$feature_number

## This part is cheesy. I am letting Tukey calculate diff for me for RA values. P_values not pulled form this test
rr3_night_filtered_RA <- rr3_night_RA%>%
  dplyr::select(1:5, c(tukey_night_sig_figs))

diff_model_rrra_night <- sapply(rr3_night_filtered_RA[6:ncol(rr3_night_filtered_RA)], function(x)
  TukeyHSD(aov(x ~ rr3_night_filtered_RA$Organism, data = rr3_night_filtered_RA)))

diff_tukey_rr3_night <- as.data.frame(diff_model_rrra_night)%>%
  rownames_to_column(var = "variable")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% '%diff%')

diff_tukey_rr3_night$feature_info <- diff_tukey_rr3_night$feature_info%>%
  gsub("X", "", .)%>%
  gsub("rr3_night_filtered_RA.Organism.diff", "diff_rr3_tukey_night", .)

sig_diff_tukey_rr3_night <- diff_tukey_rr3_night%>%
  separate(feature_info, c("feature_number", "test_info"), sep = "\\.")%>%
  unite(column_names, c(variable, test_info), sep = "_")%>%
  filter(feature_number %like any% c(tukey_night_sig_figs))

## combining diff and p_value columns
rr3_night_tukey_sig_columns <- bind_rows(sig_tukey_rr3_night, sig_diff_tukey_rr3_night)%>%
  arrange(column_names)%>%
  spread(column_names, value)

# Adding all stat columns into the feature_table_wdf ----------------------
feature_table_wdf_stats <- full_join(
  full_join(
    full_join(
      full_join(
        full_join(
          full_join(
            full_join(
              full_join(
                full_join(feature_table_wdf, dunnett_sig_columns, by = "feature_number"),
                dunnett_night_sig_columns, by = "feature_number"),
              timepoint_anova_sigs, by = "feature_number"),
            spiffy_sig_columns, by = "feature_number"),
          spiffy_rare_byreeflocal, by = "feature_number"),
        rr3_dunnett_sig_columns, by = "feature_number"),
      rr3_tukey_sig_columns, by = "feature_number"),
    rr3_night_dunnett_sig_columns, by = "feature_number"),
  rr3_night_tukey_sig_columns, by = "feature_number")


write_csv(feature_table_wdf_stats, "feature_table_master_post_stats.csv")
