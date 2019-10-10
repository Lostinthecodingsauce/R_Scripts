## Moorea2017Metabalomics removing blanks, running stats, etc.
## Written March 11th 2019 - March 27th 2019 for analyzing Mo'orea metabalomic data
## Editted specifically for RR3 September 23rd 2019 for re-anlysis and publication

# LOADING -- libraries -------------------------------------------------------
#Data mungering
library(multcomp)
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
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

# PRE-CLEANING -- SIRIUS_Zodiac elemental composition of molecular formulas -------------------------------------------
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

# PRE-CLEANING -- Canopus---------------------------

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



# PRE-CLEANING -- metadata and raw faeture table --------------------------
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


# BACKGROUND FEATURES -- flagging and removing ------------------------------------------------
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

# TRANSIENT FEATURES -- flagging and removing ---------------------------------------------
# Transient features are defined as features who's area under the curve is not more than 2E5 in at least 3 samples
# This was determined by comparing gap filled and non-gap filled data and selecting for the lowest peak area in non-gap filled data
# This gives the assumption that anything lower than 2E5 is noise. See supplemental figure
feature_table_no_background_trans_finder <- feature_table_dirty%>%
  gather(sample, xic, 2:ncol(.))%>%
  separate(sample, "experiment", sep = "_", remove = FALSE, extra = "drop")%>%
  filter(!experiment %like% "%Blank%")%>%
  mutate(experiment = case_when(experiment == "D" ~ "Dorcierr",
                                experiment == "M" ~ "Mordor",
                                experiment == "R" ~ "RR3",
                                TRUE ~ as.character(experiment)))%>%
  split(.$experiment)%>%
  map(~ spread(., sample, xic))%>%
  map(~ add_column(., trans_feature_finder = rowSums(.[3:ncol(.)] > 2E5), .before = 3))%>%
  map(~ mutate(., experiment = case_when(trans_feature_finder >= 3 ~ "real",
                                      TRUE ~ "transient")))

feature_table_no_back_trans_filter <- feature_table_dirty[1]%>%
  add_column(background_features = no_background$background, .after = 1)%>%
  add_column(SPIFFy_transient = feature_table_no_background_trans_finder$SPIFFy$experiment, .after = 2)%>%
  add_column(Mordor_transient = feature_table_no_background_trans_finder$Mordor$experiment, .after = 2)%>%
  add_column(RR3_transient = feature_table_no_background_trans_finder$RR3$experiment, .after = 2)%>%
  add_column(Dorcierr_transient = feature_table_no_background_trans_finder$Dorcierr$experiment, .after = 2)



# FILTERING -- out background and transient features ----------------------
rr3_real_features <- as.vector(feature_table_no_back_trans_filter%>%
                                 filter(background_features == "real")%>%
                                 filter(RR3_transient == "real"))$feature_number

feature_table_no_back_trans <- feature_table_dirty%>%
  gather(sample, val, 2:ncol(.))%>%
  spread(feature_number, val)%>%
  filter(sample %like any% c("R_%", "%Blank%"))%>%
  dplyr::select(c(1, rr3_real_features))%>%
  gather(feature_number, val, 2:ncol(.))%>%
  spread(sample, val)

# CLASSIFYING -- Ambient or Exudate features FOR RR3 ONLY ------------------------------
# Calculating exudates as TF/Water T0 > 1
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

ambient_exudate_features <- ambient_exudate_features_temp[3:ncol(ambient_exudate_features_temp)]%>%
  sapply(function(x) log2(x/ambient_exudate_features_temp$WA_T0))%>%
  as.data.frame()%>%
  add_column(feature_number = ambient_exudate_features_temp$feature_number, .before = 1)%>%
  add_column(DayNight = ambient_exudate_features_temp$DayNight, .before = 1)%>%
  add_column(max = apply(.[3:9],1, max))%>%
  dplyr::select(c(1:2, max))%>%
  mutate(exudate_behavior = case_when(max > 1 ~ "exudate",
                                      TRUE ~ "ambient"))%>%
  dplyr::select(-max)%>%
  spread(DayNight, exudate_behavior)%>%
  mutate(exudate_diel = case_when(D == "exudate" & N == "exudate" ~ "diel",
                                  N == "exudate" & D == "ambient" ~ "night_unique",
                                  D == "exudate" & N == "ambient" ~ "day_unique",
                                  TRUE ~ "ambient"))%>%
  mutate(exudate_behavior = case_when(exudate_diel == "ambient" ~ "ambient",
                                      TRUE ~ "exudate"))%>%
  dplyr::select(-c(D,N))


ambient_exudate_no_back_trans <- feature_table_no_back_trans%>%
  left_join(ambient_exudate_features, ., by = "feature_number")

# SUMMARY TABLE -- Background features ------------------------------------
sum_table_starting_values <- feature_table_dirty%>%
  gather(sample_name, Initial_TIC, 2:ncol(.))%>%
  add_column(Total_features = 1)%>%
  group_by(sample_name)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

sum_table_background_removal <- left_join(feature_table_no_back_trans_filter[1:2], feature_table_dirty,
                                          by = "feature_number")%>%
  filter(background_features == "real")%>%
  dplyr::select(-2)%>%
  gather(sample_name, sample_TIC_background_removal, 2:ncol(.))%>%
  add_column(background_removal_features = 1)%>%
  group_by(sample_name)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

not_background <- as.vector(no_background%>%
                              filter(background == "real"))$feature_number

sum_table_transient_removal <- feature_table_no_background_trans_finder%>%
  map(~ filter(., feature_number %in% not_background))%>%
  map(~ filter(., experiment == "real"))%>%
  map(~ dplyr::select(., -c(2:3)))%>%
  map(~ gather(., sample_name, sample_TIC, 2:ncol(.)))%>%
  map(~ add_column(., transient_removal_features = 1))%>%
  map(~ group_by(., sample_name))%>%
  map(~ summarize_if(., is.numeric, sum))%>%
  map(~ ungroup(.))

sum_table_all_transient <- bind_rows(sum_table_transient_removal$Dorcierr,
                                     sum_table_transient_removal$Mordor,
                                     sum_table_transient_removal$RR3,
                                     sum_table_transient_removal$SPIFFy)

sum_table_exudate_features <- exudate_no_back_trans%>%
  gather(sample_name, RR3_exudate_TIC, 2:ncol(.))%>%
  add_column(RR3_exudate__features = 1)%>%
  group_by(sample_name)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()


sum_table <- 
  left_join(
    left_join(
      left_join(
        sum_table_starting_values, sum_table_background_removal, by = "sample_name"),
      sum_table_all_transient, by = "sample_name"),
    sum_table_exudate_features, by = "sample_name")

write_csv(sum_table, "~/Documents/SDSU/Moorea_2017/cleaning_summary.csv")
# RELATIVIZATION AND NORMALIZATION -- grouped by ambient or exudate -----------------
feature_table_TIC <- ambient_exudate_no_back_trans%>%
  dplyr::select(-c(exudate_diel, exudate_behavior))%>%
  # group_by(exudate_behavior)%>%
  # nest()%>%
  # mutate(data = map(data, ~ 
  gather(sample_name, peak_area, 2:ncol(.))%>%
  spread(feature_number, peak_area)%>%
  add_column(TIC = apply(.[2:ncol(.)], 1, sum), .before = 2)

feature_table_relnorm <- feature_table_TIC%>%
         # transformations = map(data, ~ 
  gather(feature_number, xic, 3:(ncol(.)))%>%
  mutate(RA = .$xic/.$TIC)%>%
  mutate(asin = asin(sqrt(RA)))%>%
  dplyr::select(-TIC)%>%
  gather(transformation, values, xic:asin)%>%
  arrange(transformation)%>%
  unite(sample_transformed, c("sample_name", "transformation"), sep = "_")%>%
  spread(sample_transformed, values)%>%
  right_join(ambient_exudate_no_back_trans[1:3], ., by = "feature_number")

# PRE-CLEANING -- Build feature table working data frame ----------------------------------
# All three joined together (Peak area, RA, asin(sqrt))

exudate_feature_table_combined <- right_join(metadata, feature_table_relnorm, by = "feature_number")

exudate_table_wdf_temp <- exudate_feature_table_combined%>%
  dplyr::select(c(feature_number, exudate_diel, exudate_behavior, everything()))


# NORMALIZATION -- NOSC and energy to C -----------------------------------
# Percent is Reltaive Abundance * C content of feature Final equation is (RA*C) / (sum(RA*C)) * NOSC
carbon_normalized_ra_NOSC <- exudate_table_wdf_temp%>%
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
carbon_normalized_xic_NOSC <- exudate_table_wdf_temp%>%
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

# PRE-CLEANING -- adding carbon normalized values to wdf ------------------
exudate_table_wdf <- left_join(exudate_table_wdf_temp, carbon_normalized_NOSC, by = 'feature_number')

write_csv(exudate_table_wdf, "RR3_exudate_feature_table_master_post_filtered.csv")

# PRE-CLEANING -- Making Mo’orea working data frame for stats -----------------------------
rr3_transposed <- exudate_table_wdf%>%
  filter(exudate_behavior == "exudate")%>%
  dplyr::select(c(feature_number, ends_with("_asin")))%>%
  gather(sample_ID, angular, 2:ncol(.))%>%
  spread(feature_number, angular)

rr3_transposed$sample_ID <- rr3_transposed$sample_ID%>%
  gsub("_asin", "", .)%>%
  gsub("-", "_", .)

blanks_wdf <- rr3_transposed%>%
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

rr3_wdf <- rr3_transposed%>%
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


# SET -- SEED ----------------------------------------------------------------
set.seed(2005)

# STATS TWO-WAY ANOVA— model ---------------------------------------------------
# This line makes the names of the rows which will be added into the pvalue table
twoway_anova_rows_rr3 <- c("Organism", "DayNight", "Organism*DayNight")

# The actual Model and collection of f_values
aov_rr3 <- sapply(rr3_wdf[6:ncol(rr3_wdf)], 
                  function(x) summary(
                    aov(x ~ rr3_wdf[["Organism"]]*rr3_wdf[["DayNight"]]))[[1]][1:3,'Pr(>F)'])

two_way_rr3_tidy <- as.data.frame(aov_rr3)%>%
  add_column(anova_test = twoway_anova_rows_rr3, .before =1)%>%
  gather(feature, f_value, 2:ncol(.))%>%
  add_column(FDR = p.adjust(.$f_value, method = "BH"))%>%
  filter(FDR < 0.05)

# Significant 
significant_organism_rr3 <- as.vector(two_way_rr3_tidy%>% filter(!anova_test == "DayNight"))$feature
significant_Timepoint_rr3 <- as.vector(two_way_rr3_tidy%>% filter(!anova_test == "Organism"))$feature


# STATS POST-HOC-- dunnetts and Tukey Organism -------------------------------------
rr3_organism_post_hoc <- rr3_wdf%>%
  dplyr::select(c(1:5, significant_organism_rr3))

organism_order <- as.vector(levels(relevel(as.factor(rr3_organism_post_hoc$Organism),"Water control")))

rr3_organism_post_hoc <- rr3_organism_post_hoc%>%
  filter(Timepoint == "TF")%>%
  unite(sample, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_number, asin, 5:ncol(.))%>%
  dplyr::select(-c(Experiment, Timepoint))%>%
  group_by(DayNight, feature_number)%>%
  mutate(sum = sum(asin))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  separate(sample, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  mutate(Organism = fct_relevel(Organism, organism_order))%>%
  nest()%>%
  mutate(anova = map(data, ~ aov(asin ~ Organism, .x)),
         dunnett = map(anova, ~ glht(.x, linfct = mcp(Organism = "Dunnett"))),
         dunnett_summary = map(dunnett, ~summary(.x)%>%
                                 tidy()),
         tukey = map(anova, ~ TukeyHSD(.x, p.adjust.methods = "BH")))

# POST STATS -- Dunnetts summary table ------------------------------------
dunnets_pvalues <- rr3_organism_post_hoc%>%
  dplyr::select(c(DayNight, feature_number, dunnett_summary))%>%
  unnest(dunnett_summary)%>%
  add_column(FDR = p.adjust(.$p.value, method = "BH"))%>%
  filter(FDR < 0.05)%>%
  mutate(lhs = gsub(" - Water control", "_dunnetts", lhs))%>%
  dplyr::select(c(1,2,"lhs", "FDR"))%>%
  arrange(lhs)%>%
  group_by(DayNight)%>%
  spread(lhs, FDR)%>%
  ungroup()%>%
  group_by(DayNight, feature_number)%>%
  add_column(number_exudate_organisms = rowSums(.[3:ncol(.)] >= 0, na.rm = TRUE))%>%
  mutate(exudate_type = case_when(is.na(CCA_dunnetts) == FALSE & 
                                    is.na(Dictyota_dunnetts) &
                                    is.na(Turf_dunnetts) &
                                    is.na(`Pocillopora verrucosa_dunnetts`)  &
                                    is.na(`Porites lobata_dunnetts`) ~ "CCA",
                                  is.na(CCA_dunnetts)  & 
                                    is.na(Dictyota_dunnetts) == FALSE &
                                    is.na(Turf_dunnetts) &
                                    is.na(`Pocillopora verrucosa_dunnetts`)  &
                                    is.na(`Porites lobata_dunnetts`) ~ "Dictyota",
                                  is.na(CCA_dunnetts)  & 
                                    is.na(Dictyota_dunnetts) &
                                    is.na(Turf_dunnetts) == FALSE &
                                    is.na(`Pocillopora verrucosa_dunnetts`)  &
                                    is.na(`Porites lobata_dunnetts`) ~ "Turf",
                                  is.na(CCA_dunnetts) & 
                                    is.na(Dictyota_dunnetts) &
                                    is.na(Turf_dunnetts) &
                                    is.na(`Pocillopora verrucosa_dunnetts`) == FALSE &
                                    is.na(`Porites lobata_dunnetts`) ~ "Pocillopora verrucosa",
                                  is.na(CCA_dunnetts) & 
                                    is.na(Dictyota_dunnetts) &
                                    is.na(Turf_dunnetts) &
                                    is.na(`Pocillopora verrucosa_dunnetts`)  &
                                    is.na(`Porites lobata_dunnetts`) == FALSE ~ "Porites lobata",
                                  is.na(`Pocillopora verrucosa_dunnetts`) == FALSE &
                                    is.na(CCA_dunnetts) == FALSE & 
                                    is.na(Dictyota_dunnetts) &
                                    is.na(Turf_dunnetts)  |
                                    is.na(`Porites lobata_dunnetts`) == FALSE &
                                    is.na(CCA_dunnetts) == FALSE &
                                    is.na(Dictyota_dunnetts) &
                                    is.na(Turf_dunnetts) ~ "Corraline",
                                  is.na(`Pocillopora verrucosa_dunnetts`) &
                                    is.na(CCA_dunnetts) == FALSE & 
                                    is.na(Dictyota_dunnetts) == FALSE &
                                    is.na(`Porites lobata_dunnetts`) |
                                    is.na(CCA_dunnetts) == FALSE &
                                    is.na(Turf_dunnetts) == FALSE &
                                    is.na(`Porites lobata_dunnetts`) &
                                    is.na(`Pocillopora verrucosa_dunnetts`) ~ "Algae",
                                  is.na(`Dictyota_dunnetts`) &
                                    is.na(CCA_dunnetts)  &
                                    is.na(Turf_dunnetts) &
                                    is.na(`Porites lobata_dunnetts`) == FALSE &
                                    is.na(`Pocillopora verrucosa_dunnetts`) == FALSE ~ "Coral",
                                  is.na(`Dictyota_dunnetts`) == FALSE &
                                    is.na(CCA_dunnetts)  &
                                    is.na(Turf_dunnetts) == FALSE &
                                    is.na(`Porites lobata_dunnetts`) &
                                    is.na(`Pocillopora verrucosa_dunnetts`) ~ "Fleshy Algae",
                                  is.na(CCA_dunnetts)  & 
                                    is.na(Dictyota_dunnetts) &
                                    is.na(Turf_dunnetts) == FALSE &
                                    is.na(`Pocillopora verrucosa_dunnetts`)  &
                                    is.na(`Porites lobata_dunnetts`) ~ "Turf",
                                  number_exudate_organisms > 3 ~ "Primary Producers",
                                  TRUE ~ "Cosmo"))%>%
  ungroup()%>%
  gather(label, value, 3:ncol(.))%>%
  unite(label_b, c(label, DayNight), sep = "_")%>%
  spread(label_b, value)
  
  

# POST STATS -- Tulkey summary table --------------------------------------
## Not yet implimented, however the data has been analyzed
# tukey_pvalues <- rr3_organism_post_hoc%>%
#   dplyr::select(c(DayNight, feature_number, tukey))
#   unnest(tukey)
#   as.data.frame()
#   filter(FDR < 0.05)

# STATS -- PERMANOVA ------------------------------------------------------
permanovas_one_way <- rr3_wdf%>%
  filter(Timepoint == "TF")%>%
  group_by(DayNight)%>%
  nest()%>%
  mutate(permanova_by_DayNight = map(
    data, ~ adonis(.x[6:ncol(.x)] ~ Organism, .x, perm=1000, method="bray", p.adjust.methods = "BH")))

rr3_tf <- rr3_wdf%>%
  filter(Timepoint == "TF")

permanova_two_way <- adonis(rr3_tf[6:ncol(rr3_tf)] ~ Organism*DayNight, rr3_tf, perm = 1000,  p.adjust.methods = "BH")

# These will display the actual summary tables
permanova_two_way

permanovas_one_way$permanova_by_DayNight


# VISUALIZATION -- PCoA ---------------------------------------------------
# Ungrouped PCoA
rr3_pcoa <- rr3_wdf%>%
  filter(Timepoint == "TF")%>%
  dplyr::select(-c(1:5))%>%
  vegdist("bray")%>%
  pcoa()

#Eigen plot
ylimit = c(0, 1.1*max(rr3_pcoa$values$Relative_eig))
  
Eigan <- barplot(rr3_pcoa$values$Relative_eig[1:10], ylim= ylimit)
  
text(x = Eigan, y = rr3_pcoa$values$Relative_eig[1:10], label = rr3_pcoa$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

#PCoA plot
pco_scores_rr3 <- as.data.frame(rr3_pcoa$vectors)%>%
  add_column(Organism = filter(rr3_wdf, Timepoint == "TF")$Organism,
             test = "both",
             DayNight = filter(rr3_wdf, Timepoint == "TF")$DayNight,
             Replicate = filter(rr3_wdf, Timepoint == "TF")$Replicate , .before = 1)

# Effect of both Day and Night on organism (grouped)
rr3_pcoa_grouped <- rr3_wdf%>%
  filter(Timepoint == "TF")%>%
  split(.$DayNight)%>%
  map(~ dplyr::select(., -c(1:5)))%>%
  map(~ vegdist(., "bray")%>%
                     pcoa(.))

pco_night <- rr3_pcoa_grouped$Night$vectors%>%
  as.data.frame()%>%
  add_column(Organism = filter(rr3_wdf, Timepoint == "TF", DayNight == "Night")$Organism,
             test = filter(rr3_wdf, Timepoint == "TF", DayNight == "Night")$DayNight,
             DayNight = filter(rr3_wdf, Timepoint == "TF", DayNight == "Night")$DayNight,
             Replicate = filter(rr3_wdf, Timepoint == "TF", DayNight == "Night")$Replicate , .before = 1)

pco_day <- rr3_pcoa_grouped$Day$vectors%>%
  as.data.frame()%>%
  add_column(Organism = filter(rr3_wdf, Timepoint == "TF", DayNight == "Day")$Organism,
             test = filter(rr3_wdf, Timepoint == "TF", DayNight == "Day")$DayNight,
             DayNight = filter(rr3_wdf, Timepoint == "TF", DayNight == "Day")$DayNight,
             Replicate = filter(rr3_wdf, Timepoint == "TF", DayNight == "Day")$Replicate , .before = 1)


all_rr3_pcoa <- bind_rows(pco_scores_rr3, pco_day, pco_night)%>%
  mutate(rrshape = case_when(DayNight == "Day" ~ 16,
                             TRUE ~ (1)))


#Visualization 
pdf("all_plots.pdf", height = 5, width = 7)
all_rr3_pcoa%>%
  split(.$test)%>%
  map(~ ggplot(., aes(x = Axis.1, y = Axis.2, color = Organism, shape = DayNight)) +
        geom_point(stat = "identity", aes(size = 3)) +
        scale_shape_manual(values = .$rrshape) +
        ggtitle(.$test) +
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
dev.off()
  

# FINAL TABLES -- feature_table_post_stats------------------------------------------------
feature_table_post_stats <- left_join(exudate_table_wdf, dunnets_pvalues, by = "feature_number")

write_csv(feature_table_post_stats, 
          "~/Documents/SDSU/Moorea_2017/RR3/10022019_exudate_ambient/RR3_feature_table_post_stats.csv")
# VISUALIZATIONS -- XIC values ambient vs exudates ------------------------
rr3_xic <- left_join(exudates_column, feature_table_no_back_trans, by = "feature_number")%>%
  dplyr::select(-c(3:10, 53))%>%
  gather(sample_ID, XIC, 3:ncol(.))

rr3_xic_wdf <- rr3_xic%>%
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
  filter(Timepoint == "TF")



ggplot(rr3_xic_wdf, aes(x = Organism, y = log10(XIC)))+
  geom_bar(aes(fill= max),
               stat = "sum", position = "stack") +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    # panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',colour = "black"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.text.x = element_text(angle = 75, hjust = 1,face = "italic"),
    axis.title = element_text(face = "italic"),
    strip.text = element_text(face = "italic")
  ) +
  xlab("Organism") +
  ylab("log10(XIC)")


