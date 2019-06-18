## Base metabalomics cleaning pipeline
## Written by Zach Quinlan March 2019
## Editted June 2019

## Before you can run this you have to change any of the read.tsv or read.csv 
##    to the file tpye you are wanting to read.
## If you are not using CANOPUS, SIRIUS or Zodiac make sure you comment out those lines

# Loading libraries -------------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)


## Define your dataframes ---------------------------------------------------

# True hits library AS TSV
library_hits <- "Library-hits.tsv"

# Analog hits AS TSV
analog_hit_df <- "Analog-hits.tsv"

# Node info AS TSV
node_info_df <- "Node_info.tsv"

## CANOPUS will have two dataframes
# Canopus classes which is the type of chemical-ish that the feature most likely is. AS CSV
canopus_classes <- "Canopus_classes.csv"

#Canopus chemont annotations is the dataframe with strings for each chemical fucntional groups. AS CSV
canopus_functional_groups <- "categories.canopus.strings.CSV"

#Sirius and zodiac dataframe. Someitmes will come pre-combined as it is in this pipeline. AS CSV
sirius_zodiac_raw <- "SIRIUS_Zodiac_converted.csv"

## Feature table exported from MzMine with areas under the peak. AS CSV
gap_filled_feature_table <- "Feature_table_GapFilled.csv"

# MS sample codes. AS CSV
sample_codes <- "MS_SampleCodes.csv"

## Now in the MS sample codes data frame define where the samples and blanks are
## However, add 1 to the row number.
## So if you want to select the first through the 10th row it should actually be 2:11

# define row of all samples
ions_samples <- 10:259

## define row of all blanks
ions_blanks <- c(2:8, 260)


# Reading in Dataframes ---------------------------------------------------------
# True hits and analog hits are exported TSVs from GNPS
# True hits are more strictly matched to the library
# True hits, analog hits and node info data sheets should all be included in every single run
true_hits <- read_tsv(library_hits)%>%
  rename("feature_number" = '#Scan#')

analog_hits <- read_tsv(analog_hit_dfs)%>%
  rename("feature_number" = '#Scan#')

# Node info includes networking information about each feature
node_info <- read_tsv(node_info_df)%>%
  rename('feature_number' = 'cluster index',
         'network' = 'componentindex')

# Canopus tries to classify each feature
# Canopus, Sirius and Zodiac are all computer learning programs which assign probable chemical formulas
# There are optional. If you did not run them then comment out the below lines
canopus_anotations <- read_csv(canopus_classes)

chemont_anotations <- read_csv(canopus_functional_groups)%>%
  rename('canopus_annotation' = 'name')

# Sirius and Zodiac both try to assign molecular formulas to all the features
# The first 13 columns are the top most likely annotations. If you want to also look at the less likely hits
# comment out the dplyr::select line. It will have a lot of columns 413-ish
sirius_zodiac_anotations <- read_csv(sirius_zodiac_raw)%>%
  rename(feature_number = 1)%>%
  dplyr::select(-c(14:ncol(.)))

# Feature table has all features found within the experiments and blanks
# The columns need to be changed to the actual experiment sample codes
# Feature_table_raw is the raw export from MZMine
feature_table_raw <- read_csv(gap_filled_feature_table)%>%    ##Change this to your .csv
  rename('feature_number' = 'row ID')

ms_sample_codes <- read_csv(sample_codes)%>%      ## Change this to your sample code .csv
  rename('run_code' = 'Sample ID',
         'sample_code' = 'Sample Name')

# Combining metadata and creating feature_table ---------------------------
## highest probability canopus annotation
# canopus gives many different classifications and the percent chance that the feature falls into that category
# This section pulls the feature whcih has the highest possiblity to be it.
# We want only canopus annotations which are level 5 AND above 70% probability
# Column one of the canopus anotations column 1 is the feature number so has to be removed and re-added
canopus_annotation_names <- canopus_anotations%>%
  gather(canopus_annotation, canopus_probability, 2:ncol(.))

canopus_chemonnt_tidy <- left_join(canopus_annotation_names, chemont_anotations, by = "canopus_annotation")

canopus_filtered_tidy <- canopus_chemonnt_tidy%>%
  rename('feature_number' = 'name')%>%
  group_by(feature_number)%>%
  do(filter(., canopus_probability >= 0.70))%>%
  do(filter(., level == max(level)))%>%
  do(filter(., canopus_probability == max(canopus_probability)))%>%
  do(filter(., nchar(CLASS_STRING) == max(nchar(CLASS_STRING))))%>%
  do(filter(., nchar(canopus_annotation) == max(nchar(canopus_annotation))))%>%
  ungroup()

# Combines canopus, sirus, and zodiac
super_computer_annotations <- full_join(canopus_filtered_tidy, sirius_zodiac_anotations, by = "feature_number")

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
feature_table_temp <- feature_table_raw%>%
  gather(run_code, ion_charge, 2:ncol(.))%>%
  spread(feature_number, ion_charge)

# run codes can get messed up some times and have to be re-editted so they have the same code
feature_table_temp$run_code <- feature_table_temp$run_code%>%
  gsub(".mzXML Peak area", "", .)%>%
  gsub("_MSMS", "", .)

feature_table_dirty <- left_join(ms_sample_codes, feature_table_temp, by = "run_code")%>%
  dplyr::select(-run_code)%>%
  #filter(!sample_code %like any% c("%_XAD", "%_C18", "SR%", "LoRDI%", "SE%"))%>%
  #The above line is a good place to remove any samples you dont care to look at
  #This would be the situation where you run someones samples in your run but dont want their blanks or dumb features
  gather(feature_number, ion_charge, 2:ncol(.))%>%
  spread(sample_code, ion_charge)


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
  gather(sample_name, peak_area, 2:ncol(.))%>%
  spread(feature_number, peak_area)%>%
  add_column(TIC = apply(.[2:ncol(.)], 1, sum))

# IonCharge/TIC = RA values for every feature peak area /sum(sample peak areas)
feature_relative_abundance <- 
  as.data.frame(
    sapply(
      feature_table_RA_temp[2:ncol(.)],
      function(x) x/feature_table_RA_temp$TIC))%>%
  add_column(sample_name = feature_table_RA_temp$sample_name, .before = 1)%>%
  gather(feature_number, RA, 2:ncol(.))%>%
  spread(sample_name, RA)

# angular transformation of RA
feature_asin_sqrt <-
  as.data.frame(
    sapply(
      feature_relative_abundance[2:ncol(.)],
      function(x) asin(sqrt(x))))%>%
  add_column(feature_number = feature_table_no_back_trans$feature_number)

# Build feature table working data frame ----------------------------------
# All three joined together (Peak area, RA, asin(sqrt))

feature_table_combined <- left_join(
  left_join(
    right_join(metadata, feature_table_no_back_trans, by = 'feature_number'),
    feature_relative_abundance, by = 'feature_number', suffix = c("", ".RA")),
  feature_asin_sqrt, by = "feature_number", suffix = c("", ".asin(sqrt)"))

feature_table_wdf <- feature_table_combined%>%
  dplyr::select(feature_number, everything())

write_csv(feature_table_wdf, "moorea_feature_table_master_post_filtered.csv")