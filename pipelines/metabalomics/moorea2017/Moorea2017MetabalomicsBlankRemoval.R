## Moorea2017Metabalomics removing blanks
## Written March 11 2019 for analyzing Mo'orea metabalomic data

# Loading libraries -------------------------------------------------------
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
  select(-c(14:ncol(.)))

# Feature table has all features found within the experiments and blanks
# The columns need to be changed to the actual experiment sample codes
feature_table_raw <- read_csv("Morrea_Feayures-Table_all_Gap-Filled5.csv")
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
                                full_join(true_hits, analog_hits, by = "feature_number", suffix = c("Library_", "Analog_")), 
                                by = "feature_number"), 
                      by = "feature_number")

metadata$feature_number <- as.character(metadata$feature_number)

## making feature table so we can remove blanks
# have to change the MS codes for sample codes
# selecting for Dorcierr, Mordor, RR3, Spiffy
feature_table_temp <- feature_table_raw%>%
  select(-X326)%>%
  select(-c(2:4))%>%
  rename('feature_number' = 'row ID')%>%
  gather(run_code, ion_charge, 2:ncol(.))%>%
  spread(feature_number, ion_charge)

feature_table_temp$run_code <- feature_table_temp$run_code%>%
  gsub(".mzXML Peak area", "", .)%>%
  gsub("_MSMS", "", .)

feature_table_dirty <- left_join(ms_sample_codes, feature_table_temp, by = "run_code")%>%
  select(-run_code)%>%
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


#This section finds the features where mean area under the curve across all samples is larger than .5 * max of the blanks
#The non background features are saved into a vector so that they can be filtered from the master database
no_background <-ions%>%
  filter(mean_samples, mean_samples > 2*max_blanks)

feature_table_no_background <- left_join(no_background[1], feature_table_dirty, by = "feature_number")

## Removing transient features ---------------------------------------------
# Transient features are defined as features who's area under the curve is not more than 2E5 in at least 3 samples
# This was determined by comparing gap filled and non-gap filled data and selecting for the lowest peak area in non-gap filled data
# This gives the assumption that anything lower than 2E5 is noise. See supplemental figure
feature_table_no_background <- cbind(feature_table_no_background,
                                  trans_feature_finder = 
                                    rowSums(feature_table_no_background[ions_samples] > 2E5))%>%
  filter(trans_feature_finder, trans_feature_finder >= 3)%>%
  select(-trans_feature_finder)

## Calculate Total Ion Charge (TIC) and relativize samples -----------------
## FIND WHERE SCAN NUMBER IS GOING SO THAT THESE CAN ALL BE COMIBINED AGAIN
feature_table_no_background$TIC <- apply(feature_table_no_background[ions_samples], 1, sum)

# IonCharge/TIC = RA values for every feature/sample
feature_relative_abundance <- 
  as.data.frame(
    sapply(
      feature_table_no_background[ions_samples],
      function(x) x/feature_table_no_background$TIC))%>%
  add_column(feature_number = feature_table_no_background$feature_number)

# angular transformation of RA
feature_asin_sqrt <-
  as.data.frame(
    sapply(
      feature_relative_abundance[1:250],
      function(x) asin(sqrt(x))))%>%
  add_column(feature_number = feature_table_no_background$feature_number)

# All three joined together (Peak area, RA, asin(sqrt))
# Node and network info = [1:22], CANOPUS = [23:25], SIRIUS/ZODIAC = [26:38], Library Hits = [39:65], analog hits = []
# Blanks = [109:116], Area under the curve = [117:367], RA = [369:618], asin(sqrt) = [619:868]
feature_table_combined <- 
  right_join(metadata,
            left_join(
              left_join(
                feature_table_no_background, feature_relative_abundance,
                by = "feature_number", suffix = c("" , ".RA")),
              feature_asin_sqrt, by = "feature_number", suffix = c("", ".asin(sqrt)")), by = "feature_number")

feature_table_wdf <- feature_table_combined%>%
  select(feature_number, everything())


# Making Mo’orea working data frame for stats -----------------------------

moorea_transposed <- feature_table_wdf%>%
  select(1, 619:868)%>%
  gather(sample_ID, angular, 2:251)%>%
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
  dplyr::select("Timepoint")


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

## Need to take the mean of T0 to subtract from TF because of the differences in replicate #
dorcierr_transformations <- dorcierr_wdf%>%
  select(c(1:16956))%>% 
  gather(feature_name, RA, 5:16956)%>%
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


# Dorcierr Day Exudates (T0) ----------------------------------------------
p_values_oneway_day_exudates_dorc <- 
  as.data.frame(
    sapply(
      dorcierr_exudates_day[6:ncol(dorcierr_exudates_day)],
      function(x) summary(aov(x ~ dorcierr_exudates_day[["Organism"]]))[[1]][1,'Pr(>F)']))

sig_one_way_day_exudates_dorc <- p_values_oneway_day_exudates_dorc%>%
  rownames_to_column(var = "feature_number")%>%
  rename('p-value' = 2)

dorcierr_exudates_day$Organism <- as.factor(dorcierr_exudates_day$Organism)

oneway_model_day_exudates_dorc <-
  sapply(
    dorcierr_exudates_day[6:ncol(dorcierr_exudates_day)],
    function(x)  aov(x ~ dorcierr_exudates_day[["Organism"]]))

Dunnetts_day_exudates_dorc <-
  sapply(
    c(dorcierr_exudates_day[6:ncol(dorcierr_exudates_day)], dorcierr_exudates_day$Organism),
    function(x, y)  
      summary(
        glht(
          aov(x ~ dorcierr_exudates_day[["Organism"]]),
          linfct = mcp(y = "Dunnett")))) 

dunnett_test <- dorcierr_exudates_day%>%
  gather(feature_number, asin, 6:16957)%>%
  group_by(feature_number)%>% 
    do(d_fit = glht(aov(asin ~ Organism, data = .),
            linfct = mcp(Organism = "Dunnett")))

dunnett_test <- dorcierr_exudates_day%>%
  gather(feature_number, asin, 6:13443)%>%
  group_by(feature_number)%>% 
  do(d_fit = DunnettTest(.$asin, .$Organism, control = "Water control"))


t2_fits <- grow_t2 %>% 
  group_by(cSubjectID) %>% 
  do(fitWeight_BR = lm(Weight_average_lin_interp  = Age + log(Age) + I(1/Age),
                       data = grow_t2))

