## Script Written by Zach Quinlan 06/19/19
# Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019
# Only working on daytime exudation and remineralization
# Rewritten with changes to RR3 starting pipeline 10/11/2019
# 16s rRNA amplicon seque3nce data added in a upaded 10/18/2019

# LOADING -- packages -------------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)
library(CHNOSZ)
library(furrr)
library(future)

#PCoA, PERMANOVA
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)


# CORES -- setting processors available -----------------------------------
##Only used if future_mapping
num_cores <- availableCores() -1
# don't murder your compututer and save your self a core
# this is the parellel planning step (changes global env so this is plan for all parellel 
# work unless specificed otherwise)
plan(multiprocess, workers = num_cores) #defaults to sequential process, multiprocess is one option for parellel

# LOADING -- dataframes  ------------------------------------------------------
## FCM and fDOM data
dorc_fcm_fdom <- read_xlsx("~/Documents/SDSU/DORCIERR/Datasets/DORCIERR_fDOM_FCM.xlsx")%>%
  rename(sample_name =`Sample Name of DORCIERR_FCM_Final`)%>%
  rename('DayNight' = 'Experiment')%>%
  rename(Organism = 'Organsim')%>%
  filter(DayNight == "Day")%>%
  dplyr::select(-DayNight)%>%
  mutate(Organism = case_when(Organism == "Water" ~ "Water control",
                              TRUE ~ as.character(Organism)))

#DOC data
moorea_doc <- read_xlsx("~/Documents/SDSU/DORCIERR/Datasets/MO17_ExpSummary_DOC.2018.04.04.xlsx")%>%
  dplyr::select(1:2)%>%
  rename(sample_name = 1)%>%
  rename(DOC = 2)

# True hits and analog hits are exported CSVs from GNPS
# True hits are more strictly matched to the library
true_hits <- read_tsv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/Library-hits.tsv")%>%
  rename("feature_number" = '#Scan#')

analog_hits <- read_tsv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/Analog-hits.tsv")%>%
  rename("feature_number" = '#Scan#')

# Node info includes networking information about each feature
node_info <- read_tsv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/Node_info.tsv")%>%
  rename('feature_number' = 'cluster index',
         'network' = 'componentindex')

# Canopus tries to classify each feature
canopus_anotations <- read_csv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/SIRIUS_etc/converted/Canopus_classes.csv")

chemont_anotations <- read_csv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/categories.canopus.strings.nelsonMarch2019.CSV")%>%
  rename('canopus_annotation' = 'name')

# Sirius and Zodiac both try to assign molecular formulas to all the features
sirius_zodiac_anotations <- read_csv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/SIRIUS_etc/converted/SIRIUS_Zodiac_converted.csv")%>%
  rename(feature_number = 1)%>%
  dplyr::select(-c(14:ncol(.)))


# Feature table has all features found within the experiments and blanks
# The columns need to be changed to the actual experiment sample codes
# Feature_table_raw is the raw export from MZMine
feature_table_raw <- read_csv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/Morrea_Feayures-Table_all_Gap-Filled5.csv")%>%
  rename('feature_number' = 'row ID')

ms_sample_codes <- read_csv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/Mo'orea 2017 Mass spec sample codes - Sheet1.csv")%>%
  rename('run_code' = 'Sample ID',
         'sample_code' = 'Sample Name')
#Networking information for analyzing stats
network_id <- feature_metadata%>%
  dplyr::select('feature_number', 'network', 'LibraryID', 'Compound_NameAnalog_', 'ZodiacMF', 'ZodiacScore',
                'canopus_annotation', 'level', 'canopus_probability', 'CLASS_STRING')
network_id$feature_number <- as.character(network_id$feature_number)

#16s rRNA sequences
microbe_abundance_raw <- read_tsv("~/Documents/SDSU/DORCIERR/Datasets/microbes/MCR2017.16S.Nelson.Pipeline.October2019/abundance_table_100.shared.tsv")
microbe_taxonomy <- read_tsv("~/Documents/SDSU/DORCIERR/Datasets/microbes/MCR2017.16S.Nelson.Pipeline.October2019/annotations_100.taxonomy.tsv")


# CLEANING -- SIRIUS_Zodiac elemental composition of molecular formulas -------------------------------------------
networking_elements <- sirius_zodiac_anotations%>%
  filter(!ZodiacMF == "not_explainable")%>%
  group_by(feature_number)%>% 
  do(., rownames_to_column(as.data.frame(makeup(.$ZodiacMF, multiplier = 1), var = "element")))%>%
  spread(rowname, 3)

networking_elements[is.na(networking_elements)] <- 0

networking_energy <- networking_elements%>%
  dplyr::select(c(1, 'C', 'H', 'N', 'O', 'P', 'S'))%>%
  add_column(NOSC = (-((4*.$C + .$H - 3*.$N - 2*.$O + 5*.$P - 2*.$S)/.$C)+4))%>%
  add_column(dG = 60.3-28.5*.$NOSC)

# CLEANING -- Canopus---------------------------

canopus_annotation_names <- canopus_anotations%>%
  gather(canopus_annotation, canopus_probability, 2:ncol(.))

canopus_chemonnt_tidy <- left_join(canopus_annotation_names, chemont_anotations, by = "canopus_annotation")

# SET -- CANOPUS filters --------------------------------------------------
# Canopus annotations which are above level 3 AND 80% probability

canopus_filtered_tidy <- canopus_chemonnt_tidy%>%
  rename('feature_number' = 'name')%>%
  group_by(feature_number)%>%
  do(filter(., canopus_probability >= 0.80))%>%
  do(filter(., level > 3))%>%
  # do(filter(., level == max(level)))%>%
  do(filter(., canopus_probability == max(canopus_probability)))%>%
  do(filter(., level == max(level)))%>%
  do(filter(., nchar(CLASS_STRING) == max(nchar(CLASS_STRING))))%>%
  do(filter(., nchar(canopus_annotation) == max(nchar(canopus_annotation))))%>%
  ungroup()

write_csv(canopus_filtered_tidy, "~/Documents/SDSU/Moorea_2017/190312_new_fusion/canopus_filtered_tidy.csv")

# Combines canopus, sirus, and zodiac
super_computer_annotations <- full_join(full_join(canopus_filtered_tidy, 
                                                  sirius_zodiac_anotations, by = "feature_number"),
                                        networking_energy, by = "feature_number")%>%
  add_column(`characterization scores` = .$`quality`)%>%
  mutate(`characterization scores` = case_when(`characterization scores` != "Good" ~ "Bad quality",
                                               ZodiacScore < .98 ~ "Low probability",
                                               TRUE ~ as.character(`characterization scores`)))



# CLEANING -- metadata and filter out bad samples --------------------------
## join library hits, analog hits and super computer predictions
metadata <- full_join(node_info, 
                      full_join(super_computer_annotations,
                                full_join(feature_table_raw[1:4],
                                          full_join(true_hits, analog_hits, by = "feature_number",
                                                    suffix = c("Library_", "Analog_")),
                                          by = "feature_number"),
                                by = "feature_number"),
                      by = "feature_number")%>%
  add_column(binary_ID = .$LibraryID, .before = 1)%>%
  mutate(binary_ID = case_when(binary_ID != "N/A" ~ "1",
                               Compound_NameAnalog_ != "NA" ~ "2",
                               TRUE ~ as.character(binary_ID)))%>%
  add_column(combined_ID = .$LibraryID, .before = 1)%>%
  mutate(combined_ID = case_when(binary_ID == "1" ~ LibraryID,
                                 binary_ID == "2" ~ Compound_NameAnalog_,
                                 binary_ID == "N/A" ~ canopus_annotation,
                                 TRUE ~ as.character(binary_ID)))

metadata$feature_number <- as.character(metadata$feature_number)

networking <- metadata%>%
  dplyr::select(c(feature_number, network,combined_ID, binary_ID, 
                  canopus_annotation:CLASS_STRING, `characterization scores`, C:dG))%>%
  separate(CLASS_STRING, c("level 1", "level 2", "level 3",
                           "level 4", "level 5", "level 6", "level 7", "level 8"), sep = ";")%>%
  mutate(c_temp = case_when(C > 0 ~ "C",
                            TRUE ~ "_"),
         o_temp = case_when(O > 0 ~ "O",
                            TRUE ~ "_"),
         h_temp = case_when(H > 0 ~ "H",
                            TRUE~ "_"),
         n_temp = case_when(N > 0 ~ "N",
                            TRUE ~ "_"),
         p_temp = case_when(P > 0 ~ "P",
                            TRUE ~ "_"),
         s_temp = case_when(S > 0 ~ "S",
                            TRUE ~ "_"))%>%
  unite(simplified_makeup, c("c_temp", "h_temp", "o_temp", "n_temp", "p_temp", "s_temp"), sep = "")%>%
  mutate(simplified_makeup = gsub("_","", simplified_makeup),
         simplified_makeup = case_when(`characterization scores` != "Good" ~ "uncharacterized",
                                       simplified_makeup == "" ~ "uncharacterized",
                                       TRUE ~ as.character(simplified_makeup)))
         

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
  filter(sample_code %like any% c("%Blank%","R_%", "D_%", "M_%", "SPIFFy_%"))%>%
  filter(!sample_code %like any% c("%C18%", "%XAD%"))%>%
  gather(feature_number, ion_charge, 2:ncol(.))%>%
  spread(sample_code, ion_charge)

# DEFINING -- Samples and blanks ------------------------------------------
## defining what columns the samples are in
ions_samples <- 10:259

## defining different blanks
ions_blanks <- c(2:8, 260)


# FLAGGING -- BACKGROUND FEATURES flagging and removing ------------------------------------------------
# Background features are defined as features where max(blanks) >= 0.5*mean(samples)
ions <- feature_table_dirty

ions$max_blanks <- apply(ions[ions_blanks], 1, max)

ions$mean_samples <- apply(ions[ions_samples], 1, mean)


# This section finds the features where mean area under the curve across all samples is larger than 2 * max of the blanks
# The non background features are saved into a vector so that they can be filtered from the master database
# mean(samples)*0.5 > max_blanks
no_background <-ions%>%
  mutate(mean_samples = case_when(mean_samples*0.5 > max_blanks ~ "real",
                                  TRUE ~ "background"))%>%
  rename("background" = "mean_samples")

# FLAGGING -- TRANSIENT FEATURES flagging and removing ---------------------------------------------
# Transient features are defined as features who's area under the curve is not more than 2E5 in at least 3 samples
# This was determined by comparing gap filled and non-gap filled data and selecting for the lowest peak area in non-gap filled data
# This gives the assumption that anything lower than 2E5 is noise. See supplemental figure
feature_table_no_background_trans_finder <- feature_table_dirty%>%
  gather(sample, xic, 2:ncol(.))%>%
  separate(sample, "experiment", sep = "_", remove = FALSE, extra = "drop")%>%
  filter(!experiment %like% "%Blank%")%>%
  filter(!sample %like% "%Blank%")%>%
  mutate(experiment = case_when(experiment == "D" ~ "Dorcierr_transient",
                                experiment == "M" ~ "Mordor_transient",
                                experiment == "R" ~ "RR3_transient",
                                TRUE ~ "SPIFFy_transient"))%>%
  group_by(experiment)%>%
  nest()%>%
  mutate(data = map(data, ~ spread(.x, sample, xic)%>%
                      add_column(trans_feature_finder = rowSums(.[3:ncol(.)] > 2E5), .before = 2)%>%
                      mutate(transient = case_when(trans_feature_finder >= 3 ~ "real",
                                                   TRUE ~ "transient"))%>%
                      dplyr::select(c(feature_number, transient))))%>%
  unnest(data)%>%
  spread(experiment, transient)


feature_table_no_back_trans_filter <- full_join(feature_table_no_background_trans_finder, no_background, by = "feature_number")%>%
  dplyr::select(feature_number, background, everything())

# FILTERING -- out background and transient features ----------------------
dorcierr_real_features <- as.vector(feature_table_no_back_trans_filter%>%
                                      filter(background == "real")%>%
                                      filter(Dorcierr_transient == "real"))$feature_number

feature_table_no_back_trans <- feature_table_dirty%>%
  gather(sample, val, 2:ncol(.))%>%
  spread(feature_number, val)%>%
  filter(sample %like any% c("D_%", "%Blank%"))%>%
  dplyr::select(c(1, dorcierr_real_features))%>%
  gather(feature_number, val, 2:ncol(.))%>%
  spread(sample, val)

# RELATIVIZATION AND NORMALIZATION -- RA_asin -> Not grouped by ambient or exudate -----------------
feature_table_TIC <- feature_table_no_back_trans%>%
  # group_by(exudate_behavior)%>%
  # nest()%>%
  # mutate(data = map(data, ~ 
  gather(sample_name, peak_area, 2:ncol(.))%>%
  spread(feature_number, peak_area)%>%
  add_column(TIC = apply(.[2:ncol(.)], 1, sum), .before = 2)

feature_table_relnorm <- feature_table_TIC%>%
  # transformations = map(data, ~ 
  gather(.,feature_number, xic, 3:(ncol(.)))%>%
  mutate(RA = .$xic/.$TIC)%>%
  mutate(asin = asin(sqrt(RA)))%>%
  dplyr::select(-TIC)%>%
  gather(transformation, values, xic:asin)%>%
  arrange(transformation)%>%
  unite(sample_transformed, c("sample_name", "transformation"), sep = "_")%>%
  spread(sample_transformed, values)
  # unnest(transformations)


# DORCIERR feature_table --------------------------------------------------
feature_table_combined <- right_join(metadata, feature_table_relnorm, by = "feature_number")

dorcierr_table_wdf_temp <- feature_table_combined%>%
  dplyr::select(c(feature_number, everything()))

# NORMALIZATION -- NOSC and energy to C -----------------------------------
# Percent is Reltaive Abundance * C content of feature Final equation is (RA*C) / (sum(RA*C)) * NOSC
carbon_normalized_ra_NOSC <- dorcierr_table_wdf_temp%>%
  filter(`characterization scores` == "Good")%>%
  dplyr::select(c(feature_number, C, ends_with("_RA")))%>%
  gather(sample_name, ra, 3:ncol(.))%>%
  mutate(percent_total_C = ra*C)%>%
  group_by(sample_name)%>%
  mutate(sum_c = sum(percent_total_C))%>%
  mutate(carbon_norm_temp = percent_total_C/sum_c)%>%
  ungroup()%>%
  right_join(metadata%>%
               dplyr::select(c(feature_number, NOSC)),
             .,  by = "feature_number")%>%
  mutate(carbon_normalized_NOSC = carbon_norm_temp*NOSC)%>%
  dplyr::select(c(feature_number, sample_name, carbon_normalized_NOSC))%>%
  mutate(sample_name = gsub("_RA", "_RAC", sample_name))%>%
  spread(sample_name, carbon_normalized_NOSC)

# Percent is XIC * C content of feature 
carbon_normalized_xic_NOSC <- dorcierr_table_wdf_temp%>%
  filter(`characterization scores` == "Good")%>%
  dplyr::select(c(feature_number, C, ends_with("_xic")))%>%
  gather(sample_name, xic, 3:ncol(.))%>%
  mutate(percent_total_C = xic*C)%>%
  group_by(sample_name)%>%
  mutate(sum_c = sum(percent_total_C))%>%
  mutate(carbon_norm_temp = percent_total_C/sum_c)%>%
  ungroup()%>%
  right_join(metadata%>%
               dplyr::select(c(feature_number, NOSC)),
             .,  by = "feature_number")%>%
  mutate(carbon_normalized_NOSC = carbon_norm_temp*NOSC)%>%
  dplyr::select(c(feature_number, sample_name, carbon_normalized_NOSC))%>%
  mutate(sample_name = gsub("_xic", "_xicc", sample_name))%>%
  spread(sample_name, carbon_normalized_NOSC)

carbon_normalized_NOSC <- full_join(carbon_normalized_ra_NOSC, carbon_normalized_xic_NOSC, by = "feature_number")

# CLEANING-- adding carbon normalized values to wdf ------------------
dorcierr_features_wdf <- left_join(dorcierr_table_wdf_temp, carbon_normalized_NOSC, by = 'feature_number')

write_csv(dorcierr_features_wdf, "Dorcierr_feature_table_master_post_filtered.csv")

# PRE-CLEANING -- Making Dorcierr working data frame for stats -----------------------------
dorc_transposed <- dorcierr_features_wdf%>%
  dplyr::select(c(feature_number, ends_with("_asin")))%>%
  gather(sample_ID, angular, 2:ncol(.))%>%
  spread(feature_number, angular)

dorc_transposed$sample_ID <- dorc_transposed$sample_ID%>%
  gsub("_asin", "", .)%>%
  gsub("-", "_", .)

blanks_wdf <- dorc_transposed%>%
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

dorc_wdf <- dorc_transposed%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank")%>%
  dplyr::mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       Experiment == "M" ~ "mordor",
                                       Experiment == "R" ~ "RR3",
                                       TRUE ~ as.character(Experiment)))%>%
  dplyr::mutate(Organism = case_when(Organism == "CC" ~ "CCA",
                                     Organism == "DT" ~ "Dictyota",
                                     Organism == "PL" ~ "Porites lobata",
                                     Organism == "PV" ~ "Pocillopora verrucosa",
                                     Organism == "TR" ~ "Turf",
                                     Organism == "WA" ~ "Water control",
                                     TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))
# CLEANING -- Subsetting to FCM or fDOM -----------------------------------------------
fdom_wdf <- dorc_fcm_fdom%>%
  dplyr::select(c(4:'PARAFAC6'))%>%
  filter(!Timepoint == "T1",
         !Timepoint == "T2",
         !Timepoint == "T3",
         !Timepoint == "T4",
         !Timepoint == "T5",
         !Timepoint == "T6")

fdom_wdf$Replicate <- as.character(fdom_wdf$Replicate)

fcm_wdf <- dorc_fcm_fdom%>%
  dplyr::select(c(1:7,33:ncol(.)))

# CLEANING -- xic data ----------------------------------------
## Need to find difference between T0 and TF for each organism
feature_xic <- dorcierr_features_wdf%>%
  dplyr::select(c(feature_number, ends_with("_xic")))%>%
  gather(sample_ID, xic, 2:ncol(.))%>%
  mutate(sample_ID = gsub("_xic", "", sample_ID))%>%
  mutate(sample_ID = gsub("-", "_", sample_ID))%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank",
         !Timepoint == "T0")%>%
  dplyr::mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       Experiment == "M" ~ "mordor",
                                       Experiment == "R" ~ "RR3",
                                       TRUE ~ as.character(Experiment)))%>%
  dplyr::mutate(Organism = case_when(Organism == "CC" ~ "CCA",
                                     Organism == "DT" ~ "Dictyota",
                                     Organism == "PL" ~ "Porites lobata",
                                     Organism == "PV" ~ "Pocillopora verrucosa",
                                     Organism == "TR" ~ "Turf",
                                     Organism == "WA" ~ "Water control",
                                     TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))

dorc_time <- feature_xic%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~
                      dplyr::select(.x, -Replicate)%>%
                      group_by(Organism, feature_number, Timepoint)%>%
                      summarize_if(is.numeric, mean)%>%
                      ungroup()%>%
                      spread(Timepoint, xic)%>%
                      add_column(difference = .$TF-.$T0)%>%
                      dplyr::select(-c("TF", "T0"))),
         increase_over_time = map(data, ~
                                    filter(.x, difference > 0)%>%
                                    unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)),
         decrease_over_time = map(data, ~
                                    filter(.x, difference < 0)%>%
                                    unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)))

decrease_over_time <- dorc_time%>%
  dplyr::select(c(1, decrease_over_time))%>%
  unnest(decrease_over_time)%>%
  unite(combined, c(combined, DayNight), sep = "_")

increase_over_time <- dorc_time%>%
  dplyr::select(c(1, increase_over_time))%>%
  unnest(increase_over_time)%>%
  unite(combined, c(combined, DayNight), sep = "_")

##Finding Exudates and accumulites
organism_exudates <- feature_xic%>%
  group_by(DayNight, Timepoint)%>%
  nest()%>%
  mutate(data = map(data, ~
                      dplyr::select(.x, -c(Replicate))%>%
                      group_by(Organism, feature_number)%>%
                      summarize_if(is.numeric, mean)%>%
                      ungroup()%>%
                      spread(Organism, xic)%>%
                      gather(Organism, xic, 2:6)%>%
                      mutate(difference_from_water = xic - `Water control`)%>%
                      as.data.frame()%>%
                      dplyr::select(-`Water control`)%>%
                      unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)))%>%
  unnest(data)%>%
  unite(combined, c(combined, DayNight), sep = "_")

exudate <-organism_exudates%>%
  filter(Timepoint == "T0",
         difference_from_water > 0.00)

accumulites <-organism_exudates%>%
  filter(Timepoint == "TF",
         difference_from_water > 0.00)


# CLEANING -- feature ra data for graphing --------------------------------
feature_ra <- dorcierr_features_wdf%>%
  dplyr::select(c(feature_number, ends_with("_RA")))%>%
  gather(sample_ID, ra, 2:ncol(.))%>%
  mutate(sample_ID = gsub("_RA", "", sample_ID))%>%
  mutate(sample_ID = gsub("-", "_", sample_ID))%>%
  separate(sample_ID, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_")%>%
  filter(!Experiment %like% "%Blank%",
         !Organism %like% "%Blank",
         !Timepoint == "T0")%>%
  dplyr::mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       Experiment == "M" ~ "mordor",
                                       Experiment == "R" ~ "RR3",
                                       TRUE ~ as.character(Experiment)))%>%
  dplyr::mutate(Organism = case_when(Organism == "CC" ~ "CCA",
                                     Organism == "DT" ~ "Dictyota",
                                     Organism == "PL" ~ "Porites lobata",
                                     Organism == "PV" ~ "Pocillopora verrucosa",
                                     Organism == "TR" ~ "Turf",
                                     Organism == "WA" ~ "Water control",
                                     TRUE ~ as.character(Organism)))%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))


# CLEANING -- microbe RA data  --------------------------------------------
ra_bigger_TF <- microbe_combined%>%
  dplyr::select(-c(reads, sum, asin, numOtus, sample_code))%>%
  spread(Organism, ra)%>%
  gather(Organism, ra, 13:17)%>%
  mutate(difference = ra - `Water control`)%>%
  filter(difference > 0)%>%
  dplyr::select(c(DayNight, OFGO, Organism))
    

# PRE-STATS CLEAINING -- FCM Stats prep---------------------------------------------------------------------
## This should calculate mean cells per hour per µL for each organism
## Just looks at TF - mean(T0)
fcm_t0 <- fcm_wdf%>%
  filter(Timepoint == "T0",
         !Organism == 'Influent',
         !Organism == 'Offshore')%>%
  dplyr::select(c(4:6, 11))%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  rename(`T0` = 2)

fcm_t0_th <-fcm_wdf%>%
  filter(!Timepoint == "T5",
         !Timepoint == "T6",
         !Timepoint == "TF",
         !Organism == 'Influent',
         !Organism == 'Offshore')%>%
  dplyr::select(c(4:7, 11))%>%
  add_column(log10_cells = log10(.$`Cells µL-1`))

fcm_th <- fcm_t0_th%>%
  filter(Timepoint == 'T4')%>%
  dplyr::select(-log10_cells)%>%
  spread(Timepoint, "Cells µL-1")

fcm_t7 <- fcm_wdf%>%
  dplyr::select(c(4:7, 11))%>%
  filter(Timepoint == c('TF'),
         !Organism == 'Influent')%>%
  spread(Timepoint, "Cells µL-1")

fcm_rate_t7_t0 <- left_join(fcm_t7, fcm_t0, by = "Organism")%>%
  add_column(log10_TF = log10(.$TF))%>%
  add_column(change_per_hour = (.$TF - .$T0)/48)%>%
  add_column(log10_change_per_hour = log10(.$change_per_hour))

fcm_rate_th_t0 <- left_join(fcm_th, fcm_t0, by = "Organism")%>%
  add_column(change_per_hour = (.$T4 - .$T0)/24)%>%
  add_column(log_change_per_hour = log10(.$change_per_hour))

## This is the actual dataframe to use for running FCM stats
fcm_stats_df <- left_join(fcm_t7, fcm_rate_th_t0[c(2:3,7)], by = c("Organism", "Replicate"))%>%
  add_column(log10_TF = log10(.$TF), .before = 5)

# PRE-STATS CLEANING -- fDOM and DOC --------------------------------------------------------------
fdom_log10 <-
  as.data.frame(
    sapply(
      fdom_wdf[14:20],
      function(x) log10(x)))%>%
  add_column(sample_name = fdom_wdf$sample_name, .after = 1)%>%
  add_column(Replicate = fdom_wdf$Replicate, .after = 1)%>%
  add_column(Timepoint = fdom_wdf$Timepoint, .after = 1)%>%
  add_column(Organism = fdom_wdf$Organism, .after = 1)

doc_log10 <- moorea_doc%>%
  mutate_if(is.numeric, log10)%>%
  filter(sample_name != "D_OF_1_T0N",
         sample_name != "D_IN_2_T0N",
         sample_name != "D_PL_3_TFN",
         sample_name != "D_TR_1_T0N",
         sample_name != "D_WA_2_T0D",
         sample_name != "D_WA_1_T0D",
         sample_name != "D_CC_1_T0D",
         sample_name != "D_CC_2_T0D")

fdom_doc_log10 <- left_join(fdom_log10, doc_log10, by = "sample_name")%>%
  dplyr::select(-sample_name)

# PRE-STATS CLEANING -- Metabolomics, combining with fDOM DOC ------------------------------------
dom_stats_wdf<- dorc_wdf%>%
  full_join(., fdom_doc_log10, by = c("Organism", "Timepoint", "Replicate"))%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore")%>%
  gather(feature_number, asin, 6:ncol(.))

# PRE-STATS CLEANING -- 16s -----------------------------------------------
microbe_combined <- microbe_abundance_raw%>%
  dplyr::select(-1)%>%
  mutate(Group = case_when(Group == "Dorcierr_D_DT_1_TFD" ~ "D_DT_1_TFD",
                           Group == "DORCIERR_D_WA_2_TFN" ~ "D_WA_2_TFN",
                           Group == "D_PV_2_TFN_SA504_SC704" ~ "D_PV_2_TFN",
                           Group == "D_PV_2_TFN_SA503_SC704" ~ "D_PV_3_TFN",
                           Group == "D_WA_4_TFN_SA503_SC703" ~ "D_WA_4_TFN",
                           Group == "D_WA_4_TFN_SA504_SC703" ~ "D_WA_2_TFN",
                           TRUE ~ as.character(Group)))%>%
  filter(Group %like% "%D_%")%>%
  rename(sample_code = Group)%>%
  gather(OTU, reads, 3:ncol(.))%>%
  left_join(., microbe_taxonomy, by = "OTU")%>%
  dplyr::select(-Size)%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"), sep = ";")%>%
  mutate(Class = case_when(Class %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Class)),
         Order = case_when(Class == "unclassified" ~ "",
                                  TRUE ~ as.character(Order)),
         Family = case_when(Class == "unclassified" ~ "",
                            TRUE ~ as.character(Family)),
         Genus = case_when(Class == "unclassified" ~ "",
                           TRUE ~ as.character(Genus)),
         OTU = case_when(Class == "unclassified" ~ "",
                         TRUE ~ as.character(OTU)),
         Order = case_when(Order %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                           TRUE ~ as.character(Order)),
         Family = case_when(Order == "unclassified" ~ "",
                            TRUE ~ as.character(Family)),
         Genus = case_when(Order == "unclassified" ~ "",
                           TRUE ~ as.character(Genus)),
         OTU = case_when(Order == "unclassified" ~ "",
                         TRUE ~ as.character(OTU)),
         Family = case_when(Family %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                            TRUE ~ as.character(Family)),
         Genus = case_when(Family == "unclassified" ~ "",
                           TRUE ~ as.character(Genus)),
         OTU = case_when(Family == "unclassified" ~ "",
                         TRUE ~ as.character(OTU)),
         Genus = case_when(Genus %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                           TRUE ~ as.character(Genus)),
         OTU = case_when(Genus == "unclassified" ~ "",
                         OTU %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "sp",
                         TRUE ~ as.character(OTU)))%>%
  unite(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"), sep = ";")%>%
  group_by(sample_code, Taxonomy)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  separate(Taxonomy, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"), sep = ";")%>%
  separate(sample_code, c("Experiment", "Organism", "Replicate", "Timepoint"), sep = "_", remove = FALSE)%>%
  dplyr::mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       Experiment == "M" ~ "mordor",
                                       Experiment == "R" ~ "RR3",
                                       TRUE ~ as.character(Experiment)))%>%
  dplyr::mutate(Organism = case_when(Organism == "CC" ~ "CCA",
                                     Organism == "DT" ~ "Dictyota",
                                     Organism == "PL" ~ "Porites lobata",
                                     Organism == "PV" ~ "Pocillopora verrucosa",
                                     Organism == "TR" ~ "Turf",
                                     Organism == "WA" ~ "Water control",
                                     Organism == "IN" ~ "Influent",
                                     Organism == "OF" ~ "Offshore",
                                     TRUE ~ as.character(Organism)))%>%
  unite(OFGO, c("Order", "Family", "Genus", "OTU"), sep = ";", remove = FALSE)%>%
  separate(Timepoint, c("Timepoint", "DayNight"), sep = 2)%>%
  mutate(DayNight = case_when(DayNight == "D" ~ "Day",
                              TRUE ~ "Night"))%>%
  filter(Organism != "Offshore",
         Organism != "Influent")%>%
  group_by(sample_code)%>%
  mutate(sum = sum(reads),
         ra = reads/sum,
         asin = asin(sqrt(ra)))%>%
  ungroup()


# SET SEED ----------------------------------------------------------------
set.seed(2005)

# STATS ANOVA -- DOM TWO-WAY -------------------------------------------------
# This line makes the names of the rows which will be added into the pvalue table
twoway_anova_rows <- c("Organism", "Timepoint", "Organism*Timepoint")

# The actual Model and collection of f_values
aov_dom <- dom_stats_wdf%>%
  group_by(DayNight, feature_number)%>%
  mutate(sum = sum(asin))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  ungroup()%>%
  group_by(DayNight, feature_number)%>%
  nest()%>%
  mutate(anova = map(data, ~aov(asin ~ Organism*Timepoint, .x)%>%
                       tidy()%>%
                       filter(!term == "Residuals")%>%
                       dplyr::select(term, p.value)))%>%
  dplyr::select(-data)%>%
  unnest(anova)%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

# Save siginificant values from test to a vector to filter out for post-hoc tests
# For timepoint we care about the features which differ both between organisms and which all primary producers create or eat
significant_organism_dom_day <- aov_dom%>% filter(term == "Organism", DayNight == "Day")
significant_Timepoint_dom <- aov_dom%>% filter(!term == "Organism")

# STATS ANOVA -- FCM ONE-WAY -------------------------------------------------------
aov_fcm <-sapply(fcm_stats_df[c(4,6)], function(x) summary(
  aov(x ~ fcm_stats_df[["Organism"]]))[[1]][1,'Pr(>F)'])%>%
  as.data.frame()%>%
  rename(f_value = 1)%>%
  rownames_to_column(var = "test")%>%
  mutate(FDR = p.adjust(.$f_value, method = "BH"))



# STATS ANOVA -- Microbe TWO-Way ------------------------------------------
aov_microbe <- microbe_combined%>%
  filter(Timepoint == "TF")%>%
  group_by(OFGO)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(asin ~ Organism*DayNight, .x)%>%
                       tidy()%>%
                       filter(!term == "Residuals")%>%
                       dplyr::select(term, p.value)))%>%
  dplyr::select(-data)%>%
  unnest(anova)

anova_microbe_pvalues <- aov_microbe%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

organism_significant_microbes <- as.vector(anova_microbe_pvalues%>%
                                             filter(term == "Organism"))$OFGO

DayNight_significant_microbes <- as.vector(anova_microbe_pvalues%>%
                                             filter(term != "Organism"))$OFGO

# PRE-POST-HOC CLEANING -- Organism dunnetts ------------------------------------------------------
## Making Post-Hoc dataframes filtering for significant p_values from Two-Way Anova
dom_organism_post_hoc <- dom_stats_wdf%>%
  left_join(significant_organism_dom_day[2], ., by = "feature_number")

# PRE-POST-HOC CLEANING -- Timepoint anova --------------------------------------------------
dom_timepoint_post_hoc <- dom_stats_wdf%>%
  left_join(significant_Timepoint_dom[2], ., by = "feature_number")

dom_timepoint <- dom_timepoint_post_hoc%>%
  group_by(Organism, DayNight, feature_number)%>%
  mutate(sum = sum(asin))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  ungroup()



# PRE-POST-HOC CLEANING -- Microbe Dunnetts and DayNight anova -------------------------------
mic_organism_post_hoc <- microbe_combined%>%
  filter(OFGO %in% organism_significant_microbes)%>%
  filter(Timepoint == "TF")

daynight_microbe_post_hoc <- microbe_combined%>%
  filter(OFGO %in% DayNight_significant_microbes)%>%
  filter(Timepoint == "TF")

# STATS POST-HOC -- FCM Tukeys ----------------------------------------------------------
# Tukey growth rates for the first half of the expierment
tukey_model_fcm_th <-aov(fcm_rate_th_t0$log_change_per_hour ~ fcm_rate_th_t0$Organism, data = fcm_rate_th_t0)%>%
  TukeyHSD(p.adjust.methods = "BH")

p_values_tukey_fcm_th <- as.data.frame(tukey_model_fcm_th$`fcm_rate_th_t0$Organism`)%>%
  rownames_to_column(var = "Organism")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% "%p adj%")%>%
  filter(value < 0.05)

# Tukey TF FCM
tukey_model_fcm_TF <- aov(fcm_t7$TF ~ fcm_t7$Organism, data = fcm_t7)%>%
  TukeyHSD(p.adjust.methods = "BH")

p_values_tukey_fcm_TF <- as.data.frame(tukey_model_fcm_TF$`fcm_t7$Organism`)%>%
  rownames_to_column(var = "Organism")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% "%p adj%")%>%
  filter(value < 0.05)


# STATS POST-HOC — Timepoint anovas ---------------------------------------
aov_time <- dom_timepoint%>%
  group_by(Organism, feature_number, DayNight)%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(asin ~ Timepoint, .x)%>%
                       glance()%>%
                       dplyr::select(p.value)))%>%
  dplyr::select(-data)%>%
  unnest(anova)%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))

# STATS P-VALUE -- Timepoint anovas --------------------------------------------------------
time_aov_sigs <- aov_time%>%
  filter(FDR < 0.05)%>%
  unite(combined, c(Organism, feature_number, DayNight), sep = "_", remove = TRUE)%>%
  dplyr::select(-p.value)

time_aov_sig_features <- as.vector(time_aov_sigs$combined)


# STATS POST-HOC -- T0, TF, Day and Night Dunnetts -------------------------------------------------------------
organism_order <- as.factor(dom_organism_post_hoc$Organism)%>%
  relevel("Water control")%>%
  levels()%>%
  as.vector()

org_dunnetts_exudates <- dom_organism_post_hoc%>%
  group_by(Timepoint, DayNight, feature_number)%>%
  mutate(sum = sum(asin))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  ungroup()%>%
  mutate(Organism = factor(Organism))%>%
  mutate(Organism = fct_relevel(Organism, organism_order))%>%
  group_by(Timepoint, feature_number, DayNight)%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(asin ~ Organism, .x)%>%
                         glht(linfct = mcp(Organism = "Dunnett"))),
         dunnett_summary = map(dunnett, ~summary(.x)%>%
                                 tidy()))
         # tukey = map(data, ~ aov(asin ~ Organism, .x)%>%
         #             glht(linfct = mcp(Tension = "Tukey"))))

##Doing future_mapping across processors. Seems to take longer when you have a lot of map functions.
# org_dunnetts_exudates <- dom_organism_post_hoc%>%
#   group_by(Timepoint, DayNight, feature_number)%>%
#   mutate(sum = sum(asin))%>%
#   filter(!sum == 0)%>%
#   dplyr::select(-sum)%>%
#   ungroup()%>%
#   mutate(Organism = factor(Organism))%>%
#   mutate(Organism = fct_relevel(Organism, organism_order))%>%
#   group_by(Timepoint, feature_number, DayNight)%>%
#   nest()%>%
#   mutate(dunnett = future_map(data, ~ aov(asin ~ Organism, .x)%>%
#                                 glht(linfct = mcp(Organism = "Dunnett"))),
#          dunnett_summary = future_map(dunnett, ~ summary(.x)%>%
#                                 tidy()),
#          tukey = future_map(data, ~ aov(asin ~ Organism, .x)%>%
#                               glht(linfct = mcp(tension = "Tukey")))) #switched out map with future_map



# STATS P-VALUE -- T0 Dunnetts ----------------------------
dunnets_pvals <- org_dunnetts_exudates%>%
  dplyr::select(c(1:3, dunnett_summary))%>%
  unnest(dunnett_summary)%>%
  dplyr::select(-c(5:8))%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  mutate(lhs = gsub(" - Water control", "", lhs))%>%
  unite(combined, c(lhs, feature_number, DayNight), sep = "_")

significant_exudates <- dunnets_pvals%>%
  filter(Timepoint == "T0")

significant_accumulites <- dunnets_pvals%>%
  filter(Timepoint == "TF")

# STATS POST-HOC -- Day and Night Dunnetts TF -----------------------------
organism_order_micro <- as.factor(mic_organism_post_hoc$Organism)%>%
  relevel("Water control")%>%
  levels()%>%
  as.vector()

dunnett_microbe_pvals <- mic_organism_post_hoc%>%
  group_by(DayNight, OFGO)%>%
  mutate(sum = sum(asin))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  mutate(Organism = factor(Organism))%>%
  mutate(Organism = fct_relevel(Organism, organism_order_micro))%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(asin ~ Organism, .x)%>%
                         glht(linfct = mcp(Organism = "Dunnett"))),
         dunnett_summary = map(dunnett, ~summary(.x)%>%
                                 tidy()))%>%
  dplyr::select(-c(data,dunnett))%>%
  unnest(dunnett_summary)%>%
  dplyr::select(-c(4:7))%>%
  mutate(lhs = gsub(" - Water control", "", lhs))%>%
  rename("Organism" = "lhs")%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

# STATS POST-HOC -- DayNight anova grouped by organism --------------------
daynight_microbe_pvals <- mic_organism_post_hoc%>%
  group_by(Organism, OFGO)%>%
  mutate(sum = sum(asin))%>%
  filter(sum != 0)%>%
  dplyr::select(-sum)%>%
  nest()%>%
  mutate(data = map(data, ~ aov(asin ~ DayNight, .x)%>%
                       tidy()))%>%
  unnest(data)%>%
  dplyr::select(-c(4:7))%>%
  filter(term != "Residuals")%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))%>%
  filter(FDR < 0.05)

# META-STATS -- labile compounds --------------------------------------------------------
## These are defined as features which are present in T0 for any organism and then decrease from T0 -> TF
time_sigs_decrease <- inner_join(time_aov_sigs, decrease_over_time, by = "combined")
dunnett_exudates <- inner_join(exudate, significant_exudates, by = "combined")

labile_exudates <- inner_join(time_sigs_decrease, dunnett_exudates, by = "combined", suffix = c(".time", ".dunnetts"))%>%
  separate(combined, c("Organism", "feature_number", "DayNight"), sep = "_")%>%
  dplyr::select(1:3, FDR.dunnetts)%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~ spread(.x, Organism, FDR.dunnetts)%>%
                      add_column(number_exudate_organisms = rowSums(.[2:ncol(.)] >= 0, na.rm = TRUE))%>%
                      mutate(exudate_type = case_when(is.na(CCA) == FALSE & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) ~ "CCA",
                                                      is.na(CCA)  & 
                                                        is.na(Dictyota) == FALSE &
                                                        is.na(Turf) &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) ~ "Dictyota",
                                                      is.na(CCA)  & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) == FALSE &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) ~ "Turf",
                                                      is.na(CCA) & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) &
                                                        is.na(`Pocillopora verrucosa`) == FALSE &
                                                        is.na(`Porites lobata`) ~ "Pocillopora verrucosa",
                                                      is.na(CCA) & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) == FALSE ~ "Porites lobata",
                                                      is.na(`Pocillopora verrucosa`) == FALSE &
                                                        is.na(CCA) == FALSE & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf)  |
                                                        is.na(`Porites lobata`) == FALSE &
                                                        is.na(CCA) == FALSE &
                                                        is.na(Dictyota) &
                                                        is.na(Turf) ~ "Corraline",
                                                      is.na(`Pocillopora verrucosa`) &
                                                        is.na(CCA) == FALSE & 
                                                        is.na(Dictyota) == FALSE &
                                                        is.na(`Porites lobata`) |
                                                        is.na(CCA) == FALSE &
                                                        is.na(Turf) == FALSE &
                                                        is.na(`Porites lobata`) &
                                                        is.na(`Pocillopora verrucosa`) ~ "Algae",
                                                      is.na(`Dictyota`) &
                                                        is.na(CCA)  &
                                                        is.na(Turf) &
                                                        is.na(`Porites lobata`) == FALSE &
                                                        is.na(`Pocillopora verrucosa`) == FALSE ~ "Coral",
                                                      is.na(`Dictyota`) == FALSE &
                                                        is.na(CCA)  &
                                                        is.na(Turf) == FALSE &
                                                        is.na(`Porites lobata`) &
                                                        is.na(`Pocillopora verrucosa`) ~ "Fleshy Algae",
                                                      is.na(CCA)  & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) == FALSE &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) ~ "Turf",
                                                      number_exudate_organisms > 3 ~ "Primary Producers",
                                                      TRUE ~ "Cosmo"))))%>%
  unnest(data)%>%
  add_column(lability = "labile")



combined_labile_compounds <- right_join(networking, labile_exudates, by = "feature_number")

# ##Grouping by the labile exudates canopus annotations to understand variation
# grouped_summed_compounds <- combined_labile_compounds%>%
#   dplyr::select(c(network, canopus_annotation, exudate_type))%>%
#   add_column(value = 1)%>%
#   rownames_to_column(var = "unique")%>%
#   spread(network, value)%>%
#   dplyr::select(-unique)
# 
# grouped_summed_compounds[is.na(grouped_summed_compounds)] <- 0
# 
# summation_labile <- grouped_summed_compounds%>%
#   group_by(exudate_type, canopus_annotation)%>%
#   summarize_if(is.numeric, sum)%>%
#   ungroup()%>%
#   gather(network, value, 3:ncol(.))%>%
#   spread(exudate_type, value)
# 
# summation_labile[is.na(summation_labile)] <- 0
# 
# summation_labile$sum <- apply(summation_labile[3:ncol(summation_labile)], 1, sum)
# 
# canopus_networks_sums <- summation_labile%>%
#   filter(!sum == 0)
# 
# canopus_no_networks <- canopus_networks_sums%>%
#   dplyr::select(-network)%>%
#   group_by(canopus_annotation)%>%
#   summarize_if(is.numeric, sum)%>%
#   gather(Organism, number, 2:ncol(.))
# 
# write_csv(canopus_no_networks, "./staring_at_data/canopus_no_networks.dat")


# META-STATS -- accumulating compounds -------------
## Binning compounds different in remins by produced compounds 
time_sigs_increase <- inner_join(time_aov_sigs, increase_over_time, by = "combined")
dunnett_higher_than_h20 <- inner_join(accumulites, significant_accumulites, by = "combined")

accumulating_exudates <- inner_join(time_sigs_increase, dunnett_higher_than_h20, by = "combined", suffix = c(".time", ".dunnetts"))%>%
  separate(combined, c("Organism", "feature_number", "DayNight"), sep = "_")%>%
  dplyr::select(1:3, FDR.dunnetts)%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(data = map(data, ~ spread(.x, Organism, FDR.dunnetts)%>%
                      add_column(number_exudate_organisms = rowSums(.[2:ncol(.)] >= 0, na.rm = TRUE))%>%
                      mutate(exudate_type = case_when(is.na(CCA) == FALSE & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) ~ "CCA",
                                                      is.na(CCA)  & 
                                                        is.na(Dictyota) == FALSE &
                                                        is.na(Turf) &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) ~ "Dictyota",
                                                      is.na(CCA)  & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) == FALSE &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) ~ "Turf",
                                                      is.na(CCA) & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) &
                                                        is.na(`Pocillopora verrucosa`) == FALSE &
                                                        is.na(`Porites lobata`) ~ "Pocillopora verrucosa",
                                                      is.na(CCA) & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) == FALSE ~ "Porites lobata",
                                                      is.na(`Pocillopora verrucosa`) == FALSE &
                                                        is.na(CCA) == FALSE & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf)  |
                                                        is.na(`Porites lobata`) == FALSE &
                                                        is.na(CCA) == FALSE &
                                                        is.na(Dictyota) &
                                                        is.na(Turf) ~ "Corraline",
                                                      is.na(`Pocillopora verrucosa`) &
                                                        is.na(CCA) == FALSE & 
                                                        is.na(Dictyota) == FALSE &
                                                        is.na(`Porites lobata`) |
                                                        is.na(CCA) == FALSE &
                                                        is.na(Turf) == FALSE &
                                                        is.na(`Porites lobata`) &
                                                        is.na(`Pocillopora verrucosa`) ~ "Algae",
                                                      is.na(`Dictyota`) &
                                                        is.na(CCA)  &
                                                        is.na(Turf) &
                                                        is.na(`Porites lobata`) == FALSE &
                                                        is.na(`Pocillopora verrucosa`) == FALSE ~ "Coral",
                                                      is.na(`Dictyota`) == FALSE &
                                                        is.na(CCA)  &
                                                        is.na(Turf) == FALSE &
                                                        is.na(`Porites lobata`) &
                                                        is.na(`Pocillopora verrucosa`) ~ "Fleshy Algae",
                                                      is.na(CCA)  & 
                                                        is.na(Dictyota) &
                                                        is.na(Turf) == FALSE &
                                                        is.na(`Pocillopora verrucosa`)  &
                                                        is.na(`Porites lobata`) ~ "Turf",
                                                      number_exudate_organisms > 3 ~ "Primary Producers",
                                                      TRUE ~ "Cosmo"))))%>%
  unnest(data)%>%
  add_column(lability = "accumulite")


combined_accumulite_compounds <- right_join(networking, accumulating_exudates, by = "feature_number")



# META-STATS —- Labile + Accumulating compounds ------------------------------
combined_labile_accumulites_compounds <- bind_rows(combined_accumulite_compounds,combined_labile_compounds)%>%
  dplyr::select(feature_number, exudate_type, lability, everything())

# META-STATS -- microbes --------------------------------------------------
dunnett_micro_analysis <- dunnett_microbe_pvals%>%
  dplyr::select(-p.value)%>%
  spread(Organism, FDR)%>%
  add_column(number_exudate_organisms = rowSums(.[3:ncol(.)] >= 0, na.rm = TRUE))%>%
  mutate(microbe_organism = case_when(is.na(CCA) == FALSE & 
                                    is.na(Dictyota) &
                                    is.na(Turf) &
                                    is.na(`Pocillopora verrucosa`)  &
                                    is.na(`Porites lobata`) ~ "CCA",
                                  is.na(CCA)  & 
                                    is.na(Dictyota) == FALSE &
                                    is.na(Turf) &
                                    is.na(`Pocillopora verrucosa`)  &
                                    is.na(`Porites lobata`) ~ "Dictyota",
                                  is.na(CCA)  & 
                                    is.na(Dictyota) &
                                    is.na(Turf) == FALSE &
                                    is.na(`Pocillopora verrucosa`)  &
                                    is.na(`Porites lobata`) ~ "Turf",
                                  is.na(CCA) & 
                                    is.na(Dictyota) &
                                    is.na(Turf) &
                                    is.na(`Pocillopora verrucosa`) == FALSE &
                                    is.na(`Porites lobata`) ~ "Pocillopora verrucosa",
                                  is.na(CCA) & 
                                    is.na(Dictyota) &
                                    is.na(Turf) &
                                    is.na(`Pocillopora verrucosa`)  &
                                    is.na(`Porites lobata`) == FALSE ~ "Porites lobata",
                                  is.na(`Pocillopora verrucosa`) == FALSE &
                                    is.na(CCA) == FALSE & 
                                    is.na(Dictyota) &
                                    is.na(Turf)  |
                                    is.na(`Porites lobata`) == FALSE &
                                    is.na(CCA) == FALSE &
                                    is.na(Dictyota) &
                                    is.na(Turf) ~ "Corraline",
                                  is.na(`Pocillopora verrucosa`) &
                                    is.na(CCA) == FALSE & 
                                    is.na(Dictyota) == FALSE &
                                    is.na(`Porites lobata`) |
                                    is.na(CCA) == FALSE &
                                    is.na(Turf) == FALSE &
                                    is.na(`Porites lobata`) &
                                    is.na(`Pocillopora verrucosa`) ~ "Algae",
                                  is.na(`Dictyota`) &
                                    is.na(CCA)  &
                                    is.na(Turf) &
                                    is.na(`Porites lobata`) == FALSE &
                                    is.na(`Pocillopora verrucosa`) == FALSE ~ "Coral",
                                  is.na(`Dictyota`) == FALSE &
                                    is.na(CCA)  &
                                    is.na(Turf) == FALSE &
                                    is.na(`Porites lobata`) &
                                    is.na(`Pocillopora verrucosa`) ~ "Fleshy Algae",
                                  is.na(CCA)  & 
                                    is.na(Dictyota) &
                                    is.na(Turf) == FALSE &
                                    is.na(`Pocillopora verrucosa`)  &
                                    is.na(`Porites lobata`) ~ "Turf",
                                  number_exudate_organisms > 3 ~ "Primary Producers",
                                  TRUE ~ "Cosmo"))



# META-STATS -- correlation matrix ----------------------------------------


# GRAPHING —- PCoAs --------------------------------------
dom_graphing <- dom_stats_wdf%>%
  spread(6,7)

dorc_all_pcoa <- dom_stats_wdf%>%
  spread(6,7)%>%
  dplyr::select(-c(1:5, DOC:ncol(.)))%>%
  vegdist("bray")%>%
  pcoa()
  

#This plots Eigenvalues
#They will allow you to choose the best axes to show how your data varies
ylimit = c(0, 1.1*max(dom_pcoa_df$values$Relative_eig))

Eigan <- barplot(dom_pcoa_df$values$Relative_eig[1:10], ylim= ylimit)
# Add values to the bars
text(x = Eigan, y = dom_pcoa_df$values$Relative_eig[1:10], label = dom_pcoa_df$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

# S3 method for pcoa
pco_scores <- as.data.frame(dorc_all_pcoa$vectors)%>%
  mutate(feature = "all")%>%
  add_column(Organism = dom_graphing$Organism,
             DayNight = dom_graphing$DayNight,
             Timepoint = dom_graphing$Timepoint, .before = 1)

dorc_split <- dom_stats_wdf%>%
  spread(6,7)%>%
  dplyr::select(-c(1:5, DOC:ncol(.)))%>%
  split(.$DayNight)%>%
  map(~ dplyr::select(., -c(1:5)))%>%
  map(~ vegdist(., "bray")%>%
        pcoa(.))


pco_both <- dom_daynight%>%
  dplyr::select(-data)%>%
  unnest(pcoa)%>%
  add_column(Organism = dom_graphing$Organism)%>%
  add_column(Timepoint = dom_graphing$Timepoint)%>%
  mutate(feature = DayNight)%>%
  unite(time, c(Timepoint, DayNight), sep = "_")

pco_graph <- pco_scores%>%
  unite(time, c(Timepoint, DayNight), sep = "_")

pco_all <- bind_rows(pco_scores, pco_graph)%>%
  mutate(dorc_shape = case_when(time == "T0_Day" ~ 0,
                                time == "TF_Day" ~ 1,
                                time == "T0_Night" ~ 15,
                                time == "TF_Night" ~ 19,
                                TRUE ~ 3))

pco_all%>%
  split(.$feature)%>%
  map(~ ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism, shape = time)) +
        geom_point(stat = "identity", aes(size = 0.2)) +
        scale_shape_manual(values = .$dorc_shape) +
        ggtitle(.$feature) +
        theme(
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
          panel.grid.major.x = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
          legend.text = element_text(face = "italic")) +
        xlab("Axis 1") +
        ylab("Axis 2"))

# GRAPHING — network141_grant ---------------------------------------------
# Not visualized in 141 == "11633",
net141 <- feature_RA%>%
  dplyr::select(c(1:5, "11430", "12392",  "11716", "13535", "13491",  "16123", "15438",  "15198"))%>%
  filter(DayNight == "Day")%>%
  gather(feature_number, RA, 6:ncol(.))

net141_networking <- left_join(net141, networking, by = "feature_number")%>%
  rename("Feature Identification" = "combined_ID")

net141_networking$`Feature Identification` <- factor(net141_networking$`Feature Identification`, 
                                                     levels = c("3-Iodo-L-tyrosine", "L-3,5-Diiodotyrosine",
                                                                "N-Acetyl-L-phenylalanyl-3,5-diiodo-L-tyrosine",
                                                                "p-Fluoro-L-phenylalanine", "2-Naphthyl-D-alanine"))

ggplot(net141_networking, aes(x = reorder(Organism, - RA), y= (RA*100), fill = `Feature Identification`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values= c("#CB9E23", "#F8DF4F", "#1DACE8", "#7496D2", "#1C366B")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  facet_grid(~ Timepoint) +
  xlab("Organism") + 
  ylab("Relative Abundance (percent)")

tiff("test.tiff", units="in", width=10, height=5, res=300)
ggplot(net141_networking, aes(x = reorder(Organism, - RA), y= (RA*100), fill = `Feature Identification`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values= c("#CB9E23", "#F8DF4F", "#1DACE8", "#7496D2", "#1C366B")) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  facet_grid(~ Timepoint) +
  xlab("Organism") + 
  ylab("Relative Abundance (percent)")
dev.off()

# GRAPHING -- relative abundances of labile exudates ---------------------
labile_RA <- feature_ra%>%
  left_join(combined_labile_compounds%>%
              rename(Organism = exudate_type)%>%
              filter(number_exudate_organisms == 1),
            ., by = c("feature_number", "DayNight", "Organism"))%>%
  left_join(., networking%>%
              dplyr::select(c(feature_number, simplified_makeup)), by = "feature_number")%>%
  filter(Timepoint == "T0")%>%
  rename(`chemical composition` = simplified_makeup)%>%
  add_column(color = .$`chemical composition`)%>%
  mutate(color = case_when(color == "CHO" ~ "darkblue",
                           color == "CHON" ~ "#3B9AB2",
                           color == "CHN" ~ "#63ADBE",
                           color == "CHNP" ~ "#9EBE91",
                           color == "CHONP" ~ "#D1C74C",
                           color == "CHOP" ~ "#E4B80E",
                           color == "uncharacterized" ~ "#F21A00",
                           color == "CH" ~ "#E67D00",
                           color == "CHOS" ~ "goldenrod3",
                           color == "CHONS" ~ "goldenrod4",
                           is.na(color) ~ "red",
                           TRUE ~ as.character(color)),
         Timepoint = case_when(Timepoint == "T0" ~ "Exudate"))

labile_RA$`chemical composition` <- factor(labile_RA$`chemical composition`, 
                                                levels = c("", "CH", "CHOP", "CHONP", 
                                                           "CHNP", "CHON", "CHN", "CHO", "CHOS", 
                                                           "CHONS", "uncharacterized"))

lcolors <- labile_RA$color

names(lcolors) <- labile_RA$`chemical composition`


# pdf("labile_compound_composition.pdf")
labile_RA%>%
  ggplot(., aes(x = Organism, y= (ra*100))) +
  geom_bar(aes(fill = `chemical composition`), stat = "summary", fun.y = "sum", position = "stack") +
  scale_fill_manual(values= lcolors) +
  scale_y_continuous(limits = c(0,17.5), breaks= c(2.5, 5, 7.5, 10, 12.5, 15, 17.5)) +
theme(
  axis.text.x = element_text(angle = 60, hjust = 1),
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major.x = element_blank(), # get rid of major grid
  panel.grid.major.y = element_line(colour = "grey"),
  panel.grid.minor = element_blank(), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
) +
  facet_wrap(~DayNight) +
  xlab("Organism") + 
  ylab("Relative Abundance (percent)")
# dev.off()


# GRAPHING -- Labile accummulite ------------------------------------------
labile_accumulite <- as.vector(combined_labile_accumulites_compounds$feature_number)

labile_accumulite_RA <- feature_RA%>%
  dplyr::select(c(1:5, labile_accumulite))%>%
  gather(feature_number, RA, 6:ncol(.))%>%
  filter(DayNight == "Day",
         Timepoint == "T0")


labile_accumulite_RA_meta_full <- left_join(labile_accumulite_RA, combined_labile_accumulites_compounds%>%
                                              filter(!Organism == "algae",
                                                     !Organism == "fleshy",
                                                     !Organism == "primaryproducers")%>%
                                              unite(tag, c(Organism, Lability), sep = "_")%>%
                                              dplyr::select(1:29), 
                                            by = "feature_number")%>%
  rename(`chemical composition` = simplified_makeup)%>%
  add_column(color = .$`chemical composition`)%>%
  mutate(color = case_when(color == "CHO" ~ "darkblue",
                           color == "CHON" ~ "#3B9AB2",
                           color == "CHN" ~ "darkslatergray2", #Change color to something else
                           color == "CHNP" ~ "#9EBE91",
                           color == "CHONP" ~ "#D1C74C",
                           color == "CHOP" ~ "#E4B80E",
                           color == "" ~ "#F21A00",
                           color == "CH" ~ "#E67D00", 
                           TRUE ~ as.character(color)))%>%
  separate(tag, c("Organism-ish", "Lability"), sep = "_")

l_t0 <- labile_accumulite_RA_meta_full%>%
  filter(Timepoint == "T0" ,
         Lability == "labile")

a_tf <- labile_accumulite_RA_meta_full%>%
  filter(Timepoint == "T0",
         Lability == "accumulite")

labile_accumulite_RA_meta <- bind_rows(l_t0, a_tf)%>%
  mutate(Lability = case_when(Lability == "labile" ~ "Labile Exudate",
                              Lability == "accumulite" ~ "Microbial Accumulite",
                              TRUE ~ as.character(Lability)))

labile_accumulite_RA_meta$Lability <- factor(labile_accumulite_RA_meta$Lability, levels = c("Labile Exudate", "Microbial Accumulite"))

labile_accumulite_RA_meta$`chemical composition` <- factor(labile_accumulite_RA_meta$`chemical composition`, 
                                                           levels = c("", "CH", "CHOP", "CHONP", "CHNP", "CHON", "CHN", "CHO"))

lacolors <- labile_accumulite_RA_meta$color

names(lacolors) <- labile_accumulite_RA_meta$`chemical composition`

pdf("all_plots.pdf", height = 5, width = 7)
labile_accumulite_RA_meta%>%
  split(list(.$`Organism-ish`))%>%
  map(~ggplot(., aes(x = Organism, y= (RA*100), fill = `chemical composition`)) +
        geom_bar(stat = "summary", fun.y = "sum") +
        scale_fill_manual(values= lacolors)+
        ggtitle(unique(.$`Organism-ish`)) +
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
        ) +
        facet_wrap(~Lability) +
        xlab("Canopus Level 3") + 
        ylab("Relative Abundance (percent)"))
dev.off()



# GRAPHING -- Significant microbes ----------------------------------------
micro_sig_genera <- dunnett_micro_analysis%>%
  group_by(microbe_organism)%>%
  nest()%>%
  mutate(data = map(data, ~ left_join(.x, microbe_combined, by = c("OFGO", "DayNight"), suffix = c("_x", "_y"))%>%
                      inner_join(., ra_bigger_TF, by = c("OFGO", "DayNight", "Organism"))), 
         plots = map(data, ~ ggplot(.x, aes(Organism, ra, fill = OFGO)) +
                       geom_bar(stat = "summary", fun.y = "mean", position = "stack") +
                       ggtitle(microbe_organism) +
                       theme(
                         axis.text.x = element_text(angle = 60, hjust = 1),
                         panel.background = element_rect(fill = "transparent"), # bg of the panel
                         plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                         panel.grid.major = element_blank(), # get rid of major grid
                         panel.grid.minor = element_blank(), # get rid of minor grid
                         legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                         legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
                       ) +
                       facet_wrap(~ DayNight) +
                       ylab("Relative Abundance")
                       ))

pdf("~/Documents/SDSU/DORCIERR/Datasets/staring_at_data/microbes.pdf", height = 5, width = 6)
micro_sig_genera$plots
dev.off()

sig_genera <- dunnett_micro_analysis%>%
  filter(microbe_organism != "Primary Producers",
         microbe_organism != "Corraline",
         microbe_organism != "Fleshy Algae",
         microbe_organism != "Cosmo")%>%
  rename(Organism = microbe_organism)%>%
  left_join(., microbe_combined, by = c("OFGO", "DayNight", "Organism"), suffix = c("_x", "_y"))%>%
  inner_join(., ra_bigger_TF, by = c("OFGO", "DayNight", "Organism"))%>%
  separate(OFGO, c("Order", "Family", "Genus", "Species"), sep =";")%>%
  unite(genera, c("Family", "Genus"), sep = " ")%>%
  ggplot(., aes(Organism, ra, fill = genera)) +
  geom_bar(stat = "summary", fun.y = "mean", position = "stack") +
  theme(
    plot.margin = margin(2,.8,2,.8, "cm"),
    axis.text.x = element_text(angle = 60, hjust = 1),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  facet_wrap(~ DayNight) +
  ylab("Relative Abundance")

sig_genera

# WRITING -- Dataframe for Cytoscape ------------------------------
day_organism_mean <- feature_RA%>%
  filter(DayNight == "Day")%>%
  dplyr::select(-c("Replicate", "Experiment"))%>%
  unite(sample, c("Organism", "Timepoint"), sep = "_", remove = TRUE)%>%
  group_by(sample)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(feature_number, RA, 2:ncol(.))%>%
  spread(sample, RA)

organism_mean_networking <- right_join(networking, day_organism_mean, by = "feature_number")

microbial_accumulites_cyto <- microbial_accumulites%>%
  dplyr::select(1:3)%>%
  spread(Organism, p_value.time)

labile_exudates_cyto <- combined_labile_compounds%>%
  dplyr::select(-c(2:37))

meta_stats_cyto <- full_join(labile_exudates_cyto, microbial_accumulites_cyto, by = "feature_number", suffix = c(".labile", "accumulite"))

cyto_file<- left_join(organism_mean_networking, meta_stats_cyto, by = "feature_number")%>%
  rename(shared_name = feature_number)

write_csv(cyto_file, "./Cytoscape/Day_Dorc_cytoscape.csv")

# WRITING —- Two-Way ANOVA Significance table  -----------------------------
write_csv(two_way_signignificant, "~/Documents/SDSU/DORCIERR/Datasets/stats/Dorcierr_two_way_significant.csv")




# WRITING —- Metastats combined tables -------------------------------------
write_csv(combined_labile_compounds, "./staring_at_data/combined_labile_compounds.csv")

write_csv(combined_accumulite_compounds, "./staring_at_data/combined_accumulite_compounds.csv")

write_csv(combined_labile_accumulites_compounds, "./staring_at_data/combined_accumulating_labile_compounds.csv")

write_csv(cca_labile_accumulites, "./staring_at_data/cca_labile_accumulites.csv")

write_csv(level_2, "./Labile_table_figure/level_2.csv")

write_csv(organism_labile_level, "./Labile_table_figure/Organism_labile_table.csv")

# WRITING —- dataframes for graphing ---------------------------------------
write.csv(pco_scores, "./staring_at_data/PCo_scores.dat")

