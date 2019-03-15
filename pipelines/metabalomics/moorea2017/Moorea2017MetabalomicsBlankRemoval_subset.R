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
true_hits <- read_tsv("true_hits.tsv")%>%
  rename(scan_number = 'Scan')

analog_hits <- read_tsv("analog_hits.tsv")%>%
  rename(scan_number = 'Scan')

# Node info includes networking information about each feature
node_info <- read_csv("node_info.csv")

# Canopus tries to classify each feature
canopus_anotations <- read_csv("SIRIUS_etc/converted/Canopus_classes.csv")

# Sirius and Zodiac both try to assign molecular formulas to all the features
sirius_zodiac_anotations <- read_csv("SIRIUS_etc/converted/SIRIUS_Zodiac_converted.csv")%>%
  rename(scan_number = 'id')%>%
  select(-c('ZodiacMF2':'treeSize101'))

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
  rename('scan_number' = 1)
canopus_best_match$canopus_probability <- apply(canopus_anotations[2:ncol(canopus_anotations)], 1, max)
canopus_best_match$canopus_annotation <- colnames(canopus_anotations)[max.col(canopus_anotations[2:1288],
                                                                              ties.method="first")]
# Combines canopus, sirus, and zodiac
super_computer_annotations <- full_join(canopus_best_match, sirius_zodiac_anotations, by = "scan_number")

## join library hits, analog hits and super computer predictions
metadata <- full_join(node_info, full_join(super_computer_annotations, 
                                           full_join(true_hits, analog_hits, by = scan_number, suffix = c("true_hits", "analog_hits")), 
                                           by = scan_number), by = scan_number)

## making feature table so we can remove blanks
# have to change the MS codes for sample codes
# selecting for Dorcierr, Mordor, RR3, Spiffy
feature_table_temp <- feature_table_raw%>%
  select(-X326)%>%
  select(-c(2:4))%>%
  rename('scan_number' = 'row ID')%>%
  gather(run_code, ion_charge, 2:ncol(.))%>%
  spread(scan_number, ion_charge)

feature_table_temp$run_code <- feature_table_temp$run_code%>%
  gsub(".mzXML Peak area", "", .)%>%
  gsub("_MSMS", "", .)

feature_table_dirty <- left_join(ms_sample_codes, feature_table_temp, by = "run_code")%>%
  select(-run_code)%>%
  filter(!sample_code %like any% c("%_XAD", "%_C18", "SR%", "LoRDI%", "%LoRDI", "SE%"))%>%
  gather(scan_number, ion_charge, 2:20743)%>%
  spread(sample_code, ion_charge)

## defining what columns the samples are in
ions_samples <- 9:258

## defining different blanks
ions_blanks <- c(2:8, 259)


## flagging background features and subtraction from samples ------------------------------------------------
# The idea here is to flag and remove features where max(log10(blanks))*.5 >= mean(log10(samples))

# Log10 transformation produced -inf, had to be turned into 0's
ions_log10 <- feature_table_dirty%>%
  mutate_if(is.numeric, log10)

is.na(ions_log10) <- sapply(ions_log10, is.infinite)

ions_log10[is.na(ions_log10)] <- 0

ions_log10$max_blanks <- apply(ions_log10[2:8], 1, max)

ions_log10$mean_samples <- apply(ions_log10[10:259], 1, mean)

#This section finds the features where mean area under the curve across all samples is larger than .5 * max of the blanks
#The non background features are saved into a vector so that they can be filtered from the master database
no_background <-ions_log10%>%
  filter(mean_samples, mean_samples > 0.5*max_blanks)

feature_table_no_background <- left_join(no_background[1], feature_table_dirty, by = "scan_number")

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
  add_column(scan_number = feature_table_no_background$scan_number)

# angular transformation of RA
feature_asin_sqrt <-
  as.data.frame(
    sapply(
      feature_relative_abundance[1:250],
      function(x) asin(sqrt(x))))%>%
  add_column(scan_number = feature_table_no_background$scan_number)

# All three joined together (Ion Charge, RA, asin(sqrt))
# Ion Charge = 10:260, RA = 261:510, asin(sqrt) = 511:760
feature_table_wdf <- 
  left_join(
    left_join(
      feature_table_no_background, feature_relative_abundance, by = "scan_number", suffix = c("" , ".RA")),
    feature_asin_sqrt, by = "scan_number", suffix = c("", ".asin(sqrt)"))

