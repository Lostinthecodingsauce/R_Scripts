## Script Written by Zach Quinlan 06/19/19
# Re-organization of DORCIERR_FCM_fDOM.R because it needs to be cleaner 07/15/2019
#  Only working on daytime exudation and remineralization

# LOADING -- packages -------------------------------------------------------
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

# LOADING -- dataframes  ------------------------------------------------------
## FCM and fDOM data
dorc_fcm_fdom <- read_xlsx("DORCIERR_fDOM_FCM.xlsx")%>%
  rename(sample_name =`Sample Name of DORCIERR_FCM_Final`)%>%
  rename('DayNight' = 'Experiment')%>%
  rename(Organism = 'Organsim')%>%
  filter(DayNight == "Day")%>%
  dplyr::select(-DayNight)%>%
  mutate(Organism = case_when(Organism == "Water" ~ "Water control",
                              TRUE ~ as.character(Organism)))

#This has samples as rows features as columns
dorc_metab <- read_csv("DORCIERR_metabolomics_wdf.csv")

#OTU table-esque
feature_metadata <- read_csv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/moorea_feature_table_master_post_filtered.csv")%>%
  dplyr::select(1:100)

#DOC data
moorea_doc <- read_xlsx("MO17_ExpSummary_DOC.2018.04.04.xlsx")%>%
  dplyr::select(1:2)%>%
  rename(sample_name = 1)%>%
  rename(DOC = 2)

#Canopus_annotations
chemont_anotations <- read_csv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/categories.canopus.strings.nelsonMarch2019.CSV")

#Relative abundance data for both graphing and for finding increasing/decreasing stuff
feature_RA <- read_csv("DORCIERR_RA_wdf.dat")

#Networking information for analyzing stats
network_id <- feature_metadata%>%
  dplyr::select('feature_number', 'network', 'LibraryID', 'Compound_NameAnalog_', 'ZodiacMF', 'ZodiacScore',
                'canopus_annotation', 'level', 'canopus_probability', 'CLASS_STRING')
network_id$feature_number <- as.character(network_id$feature_number)

networking_elements <- network_id%>%
  filter(!ZodiacMF == "not_explainable")%>%
  group_by(feature_number)%>% 
  do(., rownames_to_column(as.data.frame(makeup(.$ZodiacMF, multiplier = 1), var = "element")))%>%
  spread(rowname, 3)

networking_elements[is.na(networking_elements)] <- 0

networking_energy <- networking_elements%>%
  add_column(NOSC = (-((4*.$C + .$H - 3*.$N - 2*.$O + 5*.$P - 2*.$S)/.$C)+4))%>%
  add_column(cox_gibbs_energy = 60.3-28.5*.$NOSC)

networking_close <- left_join(network_id, networking_energy, by = "feature_number")%>%
  separate(CLASS_STRING, c('Level 1','Level 2','Level 3','Level 4','Level 5',
                           'Level 6','Level 7','Level 8'), sep = ";")%>%
  add_column(binary_ID = .$LibraryID, .after = 3)%>%
  mutate(binary_ID = case_when(binary_ID != "N/A" ~ "1",
                               Compound_NameAnalog_ != "NA" ~ "2",
                               TRUE ~ as.character(binary_ID)))%>%
  add_column(combined_ID = .$LibraryID, .after = 3)%>%
  mutate(combined_ID = case_when(binary_ID == "1" ~ LibraryID,
                                 binary_ID == "2" ~ Compound_NameAnalog_,
                                 binary_ID == "N/A" ~ canopus_annotation,
                               TRUE ~ as.character(binary_ID)))%>%
  dplyr::select(-c(LibraryID, Compound_NameAnalog_))

## building column for chemical diversity
cho <- networking_close%>%
  dplyr::select(c(feature_number, 17:23))%>%
  mutate(C = case_when(C > 0 ~ "C",
                       TRUE ~ ""),
         H = case_when(H > 0 ~ "H",
                       TRUE ~ ""),
         O = case_when(O > 0 ~ "O",
                       TRUE ~ ""),
         N = case_when(N > 0 ~ "N",
                       TRUE ~ ""),
         P = case_when(P > 0 ~ "P",
                       TRUE ~ ""))%>%
  unite(simplified_makeup, c(C,H,O,N,P), sep = "")

networking <- left_join(networking_close, cho, by = "feature_number")



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

# CLEANING -- Relative Abundance data ----------------------------------------
## Need to find difference between T0 and TF for each organism

# feature_fdom_RA <- left_join(feature_RA, fdom_wdf[c(2:4,14:20)], by = c("Organism", "Timepoint", "Replicate"))

day_time <- feature_RA%>%
  filter(DayNight == "Day")%>%
  gather(feature_number, RA, 6:ncol(.))%>%
  dplyr::select(-Replicate)%>%
  group_by(Organism, feature_number, Timepoint)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(Timepoint, RA)%>%
  add_column(difference = .$TF-.$T0)%>%
  dplyr::select(-c("TF", "T0"))

increase_over_time <- day_time%>%
  filter(difference > 0)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)

decrease_over_time <- day_time%>%
  filter(difference < 0)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)

##Finding Exudates
day_organism_exudates <- feature_RA%>%
  filter(DayNight == "Day",
         Timepoint == "T0")%>%
  gather(feature_number, RA, 6:ncol(.))%>%
  dplyr::select(-c(Replicate, Timepoint))%>%
  group_by(Organism, feature_number)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(Organism, RA)

difference_from_water <- as.data.frame(sapply(day_organism_exudates[2:6], function(x) x-day_organism_exudates$`Water control`))%>%
  add_column(feature_number = day_organism_exudates$feature_number, .before = 1)%>%
  gather(Organism, difference_water, 2:6)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)

exudate <- difference_from_water%>%
  filter(difference_water > 0.00)


##Finding microbial accumalites
day_organism_accumalites <- feature_RA%>%
  filter(DayNight == "Day",
         Timepoint == "TF")%>%
  gather(feature_number, RA, 6:ncol(.))%>%
  dplyr::select(-c(Replicate, Timepoint))%>%
  group_by(Organism, feature_number)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(Organism, RA)

difference_from_water_TF <- as.data.frame(sapply(day_organism_accumalites[2:6], function(x) x-day_organism_accumalites$`Water control`))%>%
  add_column(feature_number = day_organism_accumalites$feature_number, .before = 1)%>%
  gather(Organism, difference_water, 2:6)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)

accumulites <- difference_from_water_TF%>%
  filter(difference_water > 0.00)
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

# PRE-STATS CLEANING -- Metabolomics ------------------------------------
dorcierr_remins_day <- dorc_metab%>%
  filter(DayNight == "Day")

dorcierr_remins_cleaned <- dorcierr_remins_day%>%
  dplyr::select(-c("Experiment", "DayNight"))%>%
  unite(sample, c("Organism", "Timepoint", "Replicate" ), sep = "_", remove = TRUE)%>%
  gather(feature_name, asin, 2:ncol(.))%>%
  spread(sample, asin)%>%
  add_column(sum = apply(.[2:37], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(sample, asin, 2:37)%>%
  spread(feature_name, asin)%>%
  separate(sample, c("Organism", "Timepoint", "Replicate"), sep = "_", remove = TRUE)

dorcierr_exudates_day <- dorcierr_remins_cleaned%>%
  filter(Timepoint == "T0",
         !Replicate == 3,
         !Replicate == 4)


# PRE-STATS CLEANING -- fDOM and DOC --------------------------------------------------------------
fdom_log10 <-
  as.data.frame(
    sapply(
      fdom_wdf[14:20],
      function(x) log10(x)))%>%
  add_column(sample_name = fdom_wdf$sample_name, .before = 1)%>%
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

# PRE-STATS CLEANING -- All DOM types combined ------------------------------------------------
dom_stats_wdf <- full_join(left_join(fdom_log10, doc_log10, by = "sample_name"), dorcierr_remins_cleaned, by = c("Organism", "Timepoint", "Replicate"))%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore")


# SET SEED ----------------------------------------------------------------
set.seed(2005)

# STATS ANOVA -- DOM TWO-WAY -------------------------------------------------
# This line makes the names of the rows which will be added into the pvalue table
twoway_anova_rows <- c("Organism", "Timepoint", "Organism*Timepoint")

# The actual Model and collection of f_values
aov_dom <- sapply(dom_stats_wdf[5:ncol(dom_stats_wdf)], 
                  function(x) summary(
                    aov(x ~ dom_stats_wdf[["Organism"]]*dom_stats_wdf[["Timepoint"]]))[[1]][1:3,'Pr(>F)'])


two_way_dom <- as.data.frame(aov_dom)

two_way_dom$anova_test <- cbind(twoway_anova_rows)

two_way_tidy <- two_way_dom%>%
  gather(feature, f_value, 1:12367)

## Running Bejamini-Hochberg False Discovery Rate Corrections and filtering to significant values
two_way_tidy$FDR <- p.adjust(two_way_tidy$f_value, method = "BH")

two_way_signignificant <- two_way_tidy%>%
  filter(FDR < 0.05)

# Save siginificant values from test to a vector to filter out for post-hoc tests
# For timepoint we care about the features which differ both between organisms and which all primary producers create or eat
significant_organism_dom <- as.vector(two_way_signignificant%>% filter(anova_test == "Organism"))$feature
significant_Timepoint_dom <- as.vector(two_way_signignificant%>% filter(!anova_test == "Organism"))$feature

# STATS ANOVA -- FCM ONE-WAY -------------------------------------------------------
aov_fcm <-as.data.frame(
  sapply(fcm_stats_df[c(4,6)], 
         function(x) summary(
           aov(x ~ fcm_stats_df[["Organism"]]))[[1]][1,'Pr(>F)']))%>%
  rename(f_value = 1)%>%
  rownames_to_column(var = "test")

aov_fcm$FDR <- p.adjust(aov_fcm$f_value, method = "BH")


# PRE-POST-HOC CLEANING --Organism dunnetts ------------------------------------------------------
## Making Post-Hoc dataframes filtering for significant p_values from Two-Way Anova
dom_organism_post_hoc <- dom_stats_wdf%>%
  dplyr::select(c(1:4, significant_organism_dom))

dom_organism_exudates_ph <- dom_organism_post_hoc%>%
  filter(Timepoint == "T0")%>%
  dplyr::select(-sample_name)%>%
  unite(sample, c("Organism", "Timepoint", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_number, asin, 2:ncol(.))%>%
  spread(sample, asin)%>%
  add_column(sum = apply(.[2:ncol(.)], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(sample, asin, 2:ncol(.))%>%
  spread(feature_number, asin)%>%
  separate(sample, c("Organism", "Timepoint", "Replicate"), sep = "_", remove = TRUE)


dom_organism_accumulites <- dom_organism_post_hoc%>%
  filter(Timepoint == "TF")%>%
  dplyr::select(-sample_name)%>%
  unite(sample, c("Organism", "Timepoint", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_number, asin, 2:ncol(.))%>%
  spread(sample, asin)%>%
  add_column(sum = apply(.[2:ncol(.)], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(sample, asin, 2:ncol(.))%>%
  spread(feature_number, asin)%>%
  separate(sample, c("Organism", "Timepoint", "Replicate"), sep = "_", remove = TRUE)
# PRE-POST-HOC CLEANING -- Timepoint anova --------------------------------------------------
dom_timepoint_post_hoc <- dom_stats_wdf%>%
  dplyr::select(c(1:4, significant_Timepoint_dom))

dom_timepoint <- dom_timepoint_post_hoc%>%
  dplyr::select(-sample_name)%>%
  unite(sample, c("Organism", "Timepoint", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_name, asin, 2:ncol(.))%>%
  spread(sample, asin)%>%
  add_column(sum = apply(.[2:ncol(.)], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(sample, asin, 2:ncol(.))%>%
  spread(feature_name, asin)%>%
  separate(sample, c("Organism", "Timepoint", "Replicate"), sep = "_", remove = TRUE)

## Timepoint grouped by organisms and one-way anova run
timepoint_poc <- dom_timepoint%>%
  filter(Organism == "Pocillopora verrucosa")

timepoint_por <- dom_timepoint%>%
  filter(Organism == "Porites lobata")

timepoint_dic <- dom_timepoint%>%
  filter(Organism == "Dictyota")

timepoint_trf <- dom_timepoint%>%
  filter(Organism == "Turf")

timepoint_cca <- dom_timepoint%>%
  filter(Organism == "CCA")

# STATS POST-HOC -- FCM Tukeys ----------------------------------------------------------
# Tukey growth rates for the first half of the expierment
tukey_model_fcm_th <- 
  TukeyHSD(
    aov(
      fcm_rate_th_t0$log_change_per_hour ~ fcm_rate_th_t0$Organism, data = fcm_rate_th_t0), 
    p.adjust.methods = "BH")

p_values_tukey_fcm_th <- as.data.frame(tukey_model_fcm_th$`fcm_rate_th_t0$Organism`)%>%
  rownames_to_column(var = "Organism")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% "%p adj%")%>%
  filter(value < 0.05)

# Tukey TF FCM
tukey_model_fcm_TF <- 
  TukeyHSD(
    aov(
      fcm_t7$TF ~ fcm_t7$Organism, data = fcm_t7), 
    p.adjust.methods = "BH")

p_values_tukey_fcm_TF <- as.data.frame(tukey_model_fcm_TF$`fcm_t7$Organism`)%>%
  rownames_to_column(var = "Organism")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% "%p adj%")%>%
  filter(value < 0.05)


# STATS POST-HOC — Timepoint anovas ---------------------------------------
aov_time_poc <- as.data.frame(sapply(timepoint_poc[4:ncol(timepoint_poc)], 
                                     function(x) summary(
                                       aov(x ~ timepoint_poc[["Timepoint"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Pocillopora verrucosa", .before = 1)

aov_time_por <- as.data.frame(sapply(timepoint_por[4:ncol(timepoint_por)], 
                                     function(x) summary(
                                       aov(x ~ timepoint_por[["Timepoint"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Porites lobata", .before = 1)

aov_time_dic <- as.data.frame(sapply(timepoint_dic[4:ncol(timepoint_dic)], 
                                     function(x) summary(
                                       aov(x ~ timepoint_dic[["Timepoint"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Dictyota", .before = 1)

aov_time_trf <- as.data.frame(sapply(timepoint_trf[4:ncol(timepoint_trf)], 
                                     function(x) summary(
                                       aov(x ~ timepoint_trf[["Timepoint"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Turf", .before = 1)

aov_time_cca <- as.data.frame(sapply(timepoint_cca[4:ncol(timepoint_cca)], 
                                     function(x) summary(
                                       aov(x ~ timepoint_cca[["Timepoint"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "CCA", .before = 1)

time_aovs <- bind_rows(aov_time_cca, aov_time_dic, aov_time_poc, aov_time_por, aov_time_trf)

time_aovs$FDR <- p.adjust(time_aovs$p_value, method = "BH")

# STATS P-VALUE -- Timepoint anovas --------------------------------------------------------
time_aov_sigs <- time_aovs%>%
  filter(FDR < 0.05)%>%
  # dplyr::select(-FDR)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)%>%
  dplyr::select(-p_value)%>%
  rename(p_value = FDR)



time_aov_sig_features <- as.vector(time_aov_sigs$feature_number)


# STATS POST-HOC -- T0 Dunnetts -------------------------------------------------------------
dom_org_dunnetts_exudates <- lapply(dom_organism_exudates_ph[4:ncol(dom_organism_exudates_ph)], function(y, f)
  summary(glht(aov(y ~ f, dom_organism_exudates_ph), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(dom_organism_exudates_ph$Organism), "Water control")))

pvals_dunn_org_exudates_dom <- as.data.frame(
  sapply(
    dom_org_dunnetts_exudates,
    function(x) x$test$pvalues))

pvals_dunn_org_exudates_dom <- cbind("Organism" = c("CCA",
                                                    "Dictyota", 
                                                    "Pocillopora verrucosa",
                                                    "Porites lobata",
                                                    "Turf"),
                                     pvals_dunn_org_exudates_dom)

dom_dunnetts_FDR_exudates <- pvals_dunn_org_exudates_dom%>%
  gather(feature_number, p_value, 2:ncol(.))

dom_dunnetts_FDR_exudates$FDR_f <-p.adjust(dom_dunnetts_FDR_exudates$p_value, method = "BH")

# STATS P-VALUE -- T0 Dunnetts ----------------------------
dom_dunnett_sig_exudates <- dom_dunnetts_FDR_exudates%>%
  filter(FDR_f < 0.05)%>%
  # dplyr::select(-FDR_f)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)%>%
  dplyr::select(-p_value)%>%
  rename(p_value = FDR_f)


dunnett_sig_features_exudates <- as.vector(dom_dunnett_sig_exudates$feature_number)

# STATS POST-HOC -- TF Dunnetts ----------------------------------------------------
dom_org_dunnetts_remins <- lapply(dom_organism_accumulites[4:ncol(dom_organism_accumulites)], function(y, f)
  summary(glht(aov(y ~ f, dom_organism_accumulites), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(dom_organism_accumulites$Organism), "Water control")))

pvals_dunn_org_remins_dom <- as.data.frame(
  sapply(
    dom_org_dunnetts_remins,
    function(x) x$test$pvalues))

pvals_dunn_org_remins_dom <- cbind("Organism" = c("CCA",
                                                  "Dictyota", 
                                                  "Pocillopora verrucosa",
                                                  "Porites lobata",
                                                  "Turf"),
                                   pvals_dunn_org_remins_dom)

dom_dunnetts_FDR_remins <- pvals_dunn_org_remins_dom%>%
  gather(feature_number, p_value, 2:ncol(.))

dom_dunnetts_FDR_remins$FDR_f <-p.adjust(dom_dunnetts_FDR_remins$p_value, method = "BH")


# STATS P-VALUE -- TF Dunnetts  -------------------------------------
dom_dunnett_sig_remins <- dom_dunnetts_FDR_remins%>%
  filter(FDR_f < 0.05)%>%
  # dplyr::select(-FDR_f)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)%>%
  dplyr::select(-p_value)%>%
  rename(p_value = FDR_f)



dunnett_sig_features_remin <- as.vector(dom_dunnett_sig_remins$feature_number)




# META-STATS -- labile compounds --------------------------------------------------------
## These are defined as features which are present in T0 for any organism and then decrease from T0 -> TF
time_sigs_decrease <- inner_join(time_aov_sigs, decrease_over_time, by = "combined")
dunnett_exudates <- inner_join(exudate, dom_dunnett_sig_exudates, by = "combined")

labile_exudates <- inner_join(time_sigs_decrease, dunnett_exudates, by = "combined", suffix = c(".time", ".dunnetts"))%>%
  separate(combined, c("Organism", "feature_number"), sep = "_")

labile_exudate_filtered_RA <- labile_exudates%>%
  dplyr::select(-c(4,5))%>%
  gather(test, p_value, 3:4)%>%
  unite(name, c("Organism", "test"), sep = "_", remove = TRUE)%>%
  spread(name, p_value)

labile_exudate_filtered_RA[is.na(labile_exudate_filtered_RA)] <- 0

## specific groupings
#Single organisms
labile_poc <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` != 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` == 0,
         `Pocillopora verrucosa_p_value.time`  != 0)%>%
  add_column(tag = "poc_labile")

labile_por <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` != 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` == 0,
         `Porites lobata_p_value.time`  != 0)%>%
  add_column(tag = "por_labile")


labile_cca <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` != 0,
         `Turf_p_value.dunnetts` == 0,
         `CCA_p_value.time`  != 0)%>%
  add_column(tag = "cca_labile")


labile_dic <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` != 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` == 0,
         `Dictyota_p_value.time`  != 0)%>%
  add_column(tag = "dic_labile")


labile_turf <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` != 0,
         `Turf_p_value.time`  != 0)%>%
  add_column(tag = "trf_labile")



#Classes of organisms
labile_coral <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` != 0,
         `Porites lobata_p_value.dunnetts` != 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` == 0,
         `Pocillopora verrucosa_p_value.time`  != 0,
         `Porites lobata_p_value.time`  != 0)%>%
  add_column(tag = "coral_labile")

labile_coralline <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` != 0,
         `Porites lobata_p_value.dunnetts` != 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` != 0,
         `Turf_p_value.dunnetts` == 0,
         `Pocillopora verrucosa_p_value.time`  != 0,
         `Porites lobata_p_value.time`  != 0,
         `CCA_p_value.time` != 0)%>%
  add_column(tag = "coralline_labile")

labile_algae <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` != 0,
         `CCA_p_value.dunnetts` != 0,
         `Turf_p_value.dunnetts` != 0,
         `CCA_p_value.time`  != 0,
         `Dictyota_p_value.time`  != 0,
         `Turf_p_value.time` != 0)%>%
  add_column(tag = "algae_labile")

labile_fleshy <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` != 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` != 0,
         `Dictyota_p_value.time`  != 0,
         `Turf_p_value.time` != 0)%>%
  add_column(tag = "fleshy_labile")

labile_primary <-labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` != 0,
         `Porites lobata_p_value.dunnetts` != 0,
         `Dictyota_p_value.dunnetts` != 0,
         `CCA_p_value.dunnetts` != 0,
         `Turf_p_value.dunnetts` != 0,
         `Dictyota_p_value.time`  != 0,
         `Turf_p_value.time` != 0,
         `Pocillopora verrucosa_p_value.time`  != 0,
         `Porites lobata_p_value.time`  != 0,
         `CCA_p_value.time` != 0)%>%
  add_column(tag = "primaryproducers_labile")

##Calculating bond energies of individual labile compounds
bond_energies <- bind_rows(labile_poc, labile_por, labile_cca, labile_turf, labile_dic)

bond_networks <- right_join(networking, bond_energies, by = "feature_number")%>%
  dplyr::select('tag', 'cox_gibbs_energy')

bond_tukey <- TukeyHSD(aov(
  bond_networks$cox_gibbs_energy ~ bond_networks$tag, data = bond_networks), 
  p.adjust.methods = "BH")

bond_pvalues <- as.data.frame(bond_tukey$`bond_networks$tag`)%>%
  rownames_to_column(var = "organism")


#combined labile compounds  
combined_labile <- bind_rows(labile_poc, labile_por, labile_cca, labile_turf, labile_dic, labile_coral, 
                                       labile_algae, labile_fleshy, labile_primary)

combined_labile_compounds <- right_join(networking, combined_labile, by = "feature_number")

##Grouping by the labile exudates canopus annotations to understand variation
grouped_summed_compounds <- combined_labile_compounds%>%
  dplyr::select(c(network, canopus_annotation, tag))%>%
  add_column(value = 1)%>%
  rownames_to_column(var = "unique")%>%
  spread(network, value)%>%
  dplyr::select(-unique)

grouped_summed_compounds[is.na(grouped_summed_compounds)] <- 0

summation_labile <- grouped_summed_compounds%>%
  group_by(tag, canopus_annotation)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  gather(network, value, 3:ncol(.))%>%
  spread(tag, value)

summation_labile[is.na(summation_labile)] <- 0

summation_labile$sum <- apply(summation_labile[3:ncol(summation_labile)], 1, sum)

canopus_networks_sums <- summation_labile%>%
  filter(!sum == 0)%>%
  rename(CCA = cca_labile,
         Dictyota = dic_labile,
         Pocillopora = poc_labile,
         Porites = por_labile,
         Turf = trf_labile)

canopus_no_networks <- canopus_networks_sums%>%
  dplyr::select(-network)%>%
  group_by(canopus_annotation)%>%
  summarize_if(is.numeric, sum)%>%
  gather(Organism, number, 2:ncol(.))

write_csv(canopus_no_networks, "./staring_at_data/canopus_no_networks.dat")


# META-STATS -- accumulating compounds -------------
## Binning compounds different in remins by produced compounds 
time_sigs_increase <- inner_join(time_aov_sigs, increase_over_time, by = "combined")
dunnett_higher_than_h20 <- inner_join(accumulites, dom_dunnett_sig_remins, by = "combined")

microbial_accumulites <- inner_join(time_sigs_increase, dunnett_higher_than_h20, by = "combined", suffix = c(".time", ".dunnetts"))%>%
  separate(combined, c("Organism", "feature_number"), sep = "_")%>%
  dplyr::select(-c(4:5))%>%
  gather(test, p_value, 3:4)%>%
  unite(name, c(Organism, test), sep = "_", remove = TRUE)%>%
  spread(name, p_value)

microbial_accumulites[is.na(microbial_accumulites)] <- 0

## specific groupings
#Single organisms
accumulite_poc <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` != 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` == 0,
         `Pocillopora verrucosa_p_value.time`  != 0)%>%
  add_column(tag = "poc_accumulite")

accumulite_por <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` != 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` == 0,
         `Porites lobata_p_value.time`  != 0)%>%
  add_column(tag = "por_accumulite")


accumulite_cca <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` != 0,
         `Turf_p_value.dunnetts` == 0,
         `CCA_p_value.time`  != 0)%>%
  add_column(tag = "cca_accumulite")


accumulite_dic <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` != 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` == 0,
         `Dictyota_p_value.time`  != 0)%>%
  add_column(tag = "dic_accumulite")


accumulite_turf <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` != 0,
         `Turf_p_value.time`  != 0)%>%
  add_column(tag = "trf_accumulite")



#Classes of organisms
accumulite_coral <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` != 0,
         `Porites lobata_p_value.dunnetts` != 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` == 0,
         `Pocillopora verrucosa_p_value.time`  != 0,
         `Porites lobata_p_value.time`  != 0)%>%
  add_column(tag = "coral_accumulite")

accumulite_coralline <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` != 0,
         `Porites lobata_p_value.dunnetts` != 0,
         `Dictyota_p_value.dunnetts` == 0,
         `CCA_p_value.dunnetts` != 0,
         `Turf_p_value.dunnetts` == 0,
         `Pocillopora verrucosa_p_value.time`  != 0,
         `Porites lobata_p_value.time`  != 0,
         `CCA_p_value.time` != 0)%>%
  add_column(tag = "coralline_accumulite")

accumulite_algae <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` != 0,
         `CCA_p_value.dunnetts` != 0,
         `Turf_p_value.dunnetts` != 0,
         `CCA_p_value.time`  != 0,
         `Dictyota_p_value.time`  != 0,
         `Turf_p_value.time` != 0)%>%
  add_column(tag = "algae_accumulite")

accumulite_fleshy <- microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` == 0,
         `Porites lobata_p_value.dunnetts` == 0,
         `Dictyota_p_value.dunnetts` != 0,
         `CCA_p_value.dunnetts` == 0,
         `Turf_p_value.dunnetts` != 0,
         `Dictyota_p_value.time`  != 0,
         `Turf_p_value.time` != 0)%>%
  add_column(tag = "fleshy_accumulite")

accumulite_primary <-microbial_accumulites%>%
  filter(`Pocillopora verrucosa_p_value.dunnetts` != 0,
         `Porites lobata_p_value.dunnetts` != 0,
         `Dictyota_p_value.dunnetts` != 0,
         `CCA_p_value.dunnetts` != 0,
         `Turf_p_value.dunnetts` != 0,
         `Dictyota_p_value.time`  != 0,
         `Turf_p_value.time` != 0,
         `Pocillopora verrucosa_p_value.time`  != 0,
         `Porites lobata_p_value.time`  != 0,
         `CCA_p_value.time` != 0)%>%
  add_column(tag = "primaryproducers_accumulite")

##combined accumulite data frame
combined_accumulite_compounds <- right_join(networking, bind_rows(
  accumulite_poc, accumulite_por, accumulite_cca, accumulite_turf, accumulite_dic, accumulite_coral,
  accumulite_algae, accumulite_fleshy, accumulite_primary), by = "feature_number")


# META-STATS —- Labile + Accumulating compounds ------------------------------
combined_labile_accumulites_compounds <- bind_rows(combined_accumulite_compounds,combined_labile_compounds)%>%
  separate(tag, c("Organism", "Lability"), sep = '_')%>%
  dplyr::select(feature_number, Organism, Lability, everything())

cca_labile_accumulites <- bind_rows(right_join(networking, labile_cca, by = "feature_number"),
                                    right_join(networking, accumulite_cca, by = "feature_number"))%>%
  separate(tag, c("Organism", "Lability"), sep = '_')%>%
  dplyr::select(feature_number, Organism, Lability, everything())


# META-STATS — Microbial remineralized compounds --------------------------
comparing_features_remins <- dom_dunnett_sig_remins%>%
  separate(combined, c("Organism", "feature_number"), sep = "_")%>%
  spread(Organism, p_value)

comparing_features_remins[is.na(comparing_features_remins)] <- 0

## unique features by oragnism
reminned_poc <- comparing_features_remins%>%
  filter(`Pocillopora verrucosa` != 0,
         `Porites lobata` == 0,
         `Dictyota` == 0,
         `CCA` == 0,
         `Turf` == 0)%>%
  add_column(tag = "Pocillopora_reminned")

reminned_por <- comparing_features_remins%>%
  filter(`Pocillopora verrucosa` == 0,
         `Porites lobata` != 0,
         `Dictyota` == 0,
         `CCA` == 0,
         `Turf` == 0)%>%
  add_column(tag = "Porites_reminned")

reminned_coral <- comparing_features_remins%>%
  filter(`Pocillopora verrucosa` != 0,
         `Porites lobata` != 0,
         `Dictyota` == 0,
         `CCA` == 0,
         `Turf` == 0)%>%
  add_column(tag = "coral_reminned")

reminned_dic <- comparing_features_remins%>%
  filter(`Pocillopora verrucosa` == 0,
         `Porites lobata` == 0,
         `Dictyota` != 0,
         `CCA` == 0,
         `Turf` == 0)%>%
  add_column(tag = "Dictyota_reminned")

reminned_trf <- comparing_features_remins%>%
  filter(`Pocillopora verrucosa` == 0,
         `Porites lobata` == 0,
         `Dictyota` == 0,
         `CCA` == 0,
         `Turf` != 0)%>%
  add_column(tag = "Turf_reminned")

reminned_cca <- comparing_features_remins%>%
  filter(`Pocillopora verrucosa` == 0,
         `Porites lobata` == 0,
         `Dictyota` == 0,
         `CCA` != 0,
         `Turf` == 0)%>%
  add_column(tag = "CCA_reminned")

reminned_fleshy <- comparing_features_remins%>%
  filter(`Pocillopora verrucosa` == 0,
         `Porites lobata` == 0,
         `Dictyota` != 0,
         `CCA` == 0,
         `Turf` != 0)%>%
  add_column(tag = "fleshy_reminned")


reminned_algae <- comparing_features_remins%>%
  filter(`Pocillopora verrucosa` == 0,
         `Porites lobata` == 0,
         `Dictyota` != 0,
         `CCA` != 0,
         `Turf` != 0)%>%
  add_column(tag = "algae_reminned")


# FIGURES —- Labile Table by chemical class --------------------------------
Organism_labile_table <- combined_labile_compounds%>%
  dplyr::select(tag, canopus_annotation)%>%
  mutate(tag = case_when(tag == "poc_labile" ~ "Pocillopora verrucosa",
                         tag == "por_labile" ~ "Porites lobata",
                         tag == "trf_labile" ~ "Turf",
                         tag == "cca_labile" ~ "Crustose Corraline Algae",
                         tag == "dic_labile" ~ "Dictyota",
                         tag == "coral_labile" ~ "Coral",
                         tag == "algae_labile" ~ "Algae",
                         tag == "fleshy_labile" ~ "Fleshy Algae",
                         tag == "primaryproducers_labile" ~ "Primary Producer",
                         TRUE ~ as.character(tag)))%>%
  unite(combined, c("tag", "canopus_annotation"), sep = "_")%>%
  group_by(combined)%>%
  add_column(binary = 1)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()%>%
  separate(combined, c('tag', 'Chemical Classes'), sep = "_")%>%
  spread(tag, binary)

Organism_labile_table[is.na(Organism_labile_table)] <- 0

labile_table_canopus_rename <- Organism_labile_table%>%
  rename('name' = 'Chemical Classes')

level_by_string <- right_join(chemont_anotations[c(1,6)], labile_table_canopus_rename, by = "name")%>%
  separate(CLASS_STRING, c("level 1", "level 2", "level 3",
                           "level 4", "level 5", "level 6", "level 7", "level 8"), sep = ";")%>%
  rename('Chemical Classes' = "name")

level_2 <- level_by_string%>%
  dplyr::select(-c("Chemical Classes", "level 1", "level 3", "level 4", "level 5", "level 6", "level 7", "level 8"))%>%
  group_by(`level 2`)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup()

organism_labile_level <- right_join(level_by_string[c(1,3)], Organism_labile_table, by = "Chemical Classes")

# GRAPHING —- PCoA Labile_accumulites --------------------------------------
clac_features <- as.vector(combined_labile_accumulites_compounds$feature_number)

clac_sig_features_wdf <- dom_stats_wdf%>%
  dplyr::select(c(1:4, clac_features))

dom_pcoa_df <- dom_stats_wdf%>%
  dplyr::select(-c(1:4))

veg_bray <- vegdist(dom_pcoa_df, "bray") #Bray-curtis distances

pc_scores<-pcoa(veg_bray) #Calculating scores

#This plots Eigenvalues
#They will allow you to choose the best axes to show how your data varies
ylimit = c(0, 1.1*max(pc_scores$values$Relative_eig))

Eigan <- barplot(pc_scores$values$Relative_eig[1:10], ylim= ylimit)
# Add values to the bars
text(x = Eigan, y = pc_scores$values$Relative_eig[1:10], label = pc_scores$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

# S3 method for pcoa
biplot(pc_scores, Y=NULL, col = clac_sig_features_wdf$Organism, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1)

pco_scores <- as.data.frame(pc_scores$vectors)
pco_scores$Sample.Code <- clac_sig_features_wdf$sample_name     # This will add reference labels to the PCoA scores


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
labile_compound <- as.vector(combined_labile_compounds$feature_number)

labile_RA <- feature_RA%>%
  dplyr::select(c(1:5, labile_compound))%>%
  gather(feature_number, RA, 6:ncol(.))%>%
  filter(DayNight == "Day")


labile_RA_meta <- left_join(labile_RA, combined_labile_compounds%>%
                              dplyr::select(1:28, tag), by = "feature_number")%>%
  filter(!tag == "algae_labile",
         !tag == "fleshy_labile",
         !tag == "primaryproducers_labile")%>%
  add_column(color = .$simplified_makeup)%>%
  mutate(color = case_when(color == "CHO" ~ "darkblue",
                           color == "CHON" ~ "#3B9AB2",
                           color == "CHN" ~ "#63ADBE",
                           color == "CHNP" ~ "#9EBE91",
                           color == "CHONP" ~ "#D1C74C",
                           color == "CHOP" ~ "#E4B80E",
                           color == "" ~ "#F21A00",
                           color == "CH" ~ "#E67D00",
                           TRUE ~ as.character(color)))

labile_RA_meta$simplified_makeup <- factor(labile_RA_meta$simplified_makeup, 
                                           levels = c("", "CH", "CHOP", "CHONP", "CHNP", "CHON", "CHN", "CHO"))

lcolors <- labile_RA_meta$color

names(lcolors) <- labile_RA_meta$simplified_makeup

pdf("all_plots_test.pdf")
labile_RA_meta%>%
  split(list(.$tag))%>%
  map(~ggplot(., aes(x = Organism, y= (RA*100), fill = `simplified_makeup`)) +
        geom_bar(stat = "summary", fun.y = "sum") +
        scale_fill_manual(values= lcolors)+
        ggtitle(unique(.$tag)) +
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1),
          panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
        ) +
        facet_wrap(~Timepoint) +
        xlab("Canopus Level 3") + 
        ylab("Relative Abundance (percent)"))
dev.off()



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

