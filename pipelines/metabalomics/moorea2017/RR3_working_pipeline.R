## Written 07/25/2019 by Zach Quinlan

##RR3 working Pipeline

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
rr3_wdf <- read_csv("rr3_wdf.csv")%>%
  filter(Timepoint == "TF")

rr3_RA <- read_csv("rr3_RA.csv")%>%
  filter(Timepoint == "TF")

feature_metadata <- read_csv("~/Documents/SDSU/Moorea_2017/190312_new_fusion/moorea_feature_table_master_post_filtered.csv")%>%
  dplyr::select(1:100)

#Networking information for analyzing stats
network_id <- feature_metadata%>%
  dplyr::select('feature_number', 'network', 'LibraryID', 'Compound_NameAnalog_', 'ZodiacMF', 'ZodiacScore',
                'canopus_annotation', 'level', 'canopus_probability', 'CLASS_STRING')
network_id$feature_number <- as.character(network_id$feature_number)

networking_elements <- network_id%>%
  filter(!ZodiacMF == "not_explainable",
         !ZodiacMF == "NA")%>%
  group_by(feature_number)%>% 
  do(., rownames_to_column(as.data.frame(makeup(.$ZodiacMF, multiplier = 1), var = "element")))%>%
  spread(rowname, 3)

networking_elements[is.na(networking_elements)] <- 0

networking_energy <- networking_elements%>%
  add_column(NOSC = (-((4*.$C + .$H - 3*.$N - 2*.$O + 5*.$P - 2*.$S)/.$C)+4))%>%
  add_column(cox_gibbs_energy = 60.3-28.5*.$NOSC)

networking_close <-left_join(network_id, networking_energy, by = "feature_number")%>%
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
  dplyr::select(c(feature_number, 18:24))%>%
  dplyr::select(-Cl)%>%
  mutate(C = case_when(C > 0 ~ "C",
                       TRUE ~ ""),
         H = case_when(H > 0 ~ "H",
                       TRUE ~ ""),
         O = case_when(O > 0 ~ "O",
                       TRUE ~ ""),
         N = case_when(N > 0 ~ "N",
                       TRUE ~ ""),
         P = case_when(P > 0 ~ "P",
                       TRUE ~ ""),
         S = case_when(S > 0 ~ "S",
                       TRUE ~ ""))%>%
  unite(simplified_makeup, c(C,H,O,N,P,S), sep = "")

networking <- left_join(networking_close, cho, by = "feature_number")


# CLEANING — Relative Abundance data --------------------------------------
#feature_table_RA
rr3_feature_table_RA <- rr3_RA%>%
  dplyr::select(-c("Experiment","Timepoint", "Replicate"))%>%
  unite(sample_name, c("Organism", "DayNight"))%>%
  group_by(sample_name)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(feature_number, RA, 2:ncol(.))%>%
  spread(sample_name, RA)

## Finding Day Exudates
rr3_day_organism_exudates <- rr3_RA%>%
  filter(DayNight == "Day")%>%
  gather(feature_number, RA, 6:ncol(.))%>%
  dplyr::select(-c(Replicate, Timepoint, DayNight))%>%
  group_by(Organism, feature_number)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(Organism, RA)

difference_from_water_day <- as.data.frame(sapply(rr3_day_organism_exudates[2:6], function(x) x-rr3_day_organism_exudates$`Water control`))%>%
  add_column(feature_number = rr3_day_organism_exudates$feature_number, .before = 1)%>%
  gather(Organism, difference_water, 2:6)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)

exudate_day <- difference_from_water_day%>%
  filter(difference_water > 0.00)

## Finding Night Exudates
rr3_Night_organism_exudates <- rr3_RA%>%
  filter(DayNight == "Night")%>%
  gather(feature_number, RA, 6:ncol(.))%>%
  dplyr::select(-c(Replicate, Timepoint, DayNight))%>%
  group_by(Organism, feature_number)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  spread(Organism, RA)

difference_from_water_Night <- as.data.frame(sapply(rr3_Night_organism_exudates[2:6], function(x) x-rr3_Night_organism_exudates$`Water control`))%>%
  add_column(feature_number = rr3_Night_organism_exudates$feature_number, .before = 1)%>%
  gather(Organism, difference_water, 2:6)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)

exudate_Night <- difference_from_water_Night%>%
  filter(difference_water > 0.00)


##Difference between day and night
rr3_DayNight_difference_RA <- rr3_RA%>%
  filter(!Timepoint == "T0")%>%
  dplyr::select(-c(Experiment, Timepoint))%>%
  gather(feature_number, RA, 4:ncol(.))%>%
  spread(DayNight, RA)%>%
  unite(combined, c(Organism, feature_number), sep = "_")%>%
  dplyr::select(-Replicate)%>%
  group_by(combined)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  add_column(higher_during = .$Day-.$Night)%>%
  mutate(higher_during = case_when(higher_during < 0 ~ "Night",
                                   higher_during > 0 ~ "Day",
                                   TRUE ~ "equal"))%>%
  dplyr::select(-c(Day, Night))

# SET SEED ----------------------------------------------------------------
set.seed(2005)


# STATS — Two way ANOVA ---------------------------------------------------
# This line makes the names of the rows which will be added into the pvalue table
twoway_anova_rows_rr3 <- c("Organism", "DayNight", "Organism*DayNight")

# The actual Model and collection of f_values
aov_rr3 <- sapply(rr3_wdf[6:ncol(rr3_wdf)], 
                  function(x) summary(
                    aov(x ~ rr3_wdf[["Organism"]]*rr3_wdf[["DayNight"]]))[[1]][1:3,'Pr(>F)'])


two_way_rr3 <- as.data.frame(aov_rr3)

two_way_rr3$anova_test <- cbind(twoway_anova_rows_rr3)

two_way_tidy_rr3 <- two_way_rr3%>%
  gather(feature, f_value, 1:12367)

## Running Bejamini-Hochberg False Discovery Rate Corrections and filtering to significant values
two_way_tidy_rr3$FDR <- p.adjust(two_way_tidy_rr3$f_value, method = "BH")

two_way_signignificant_rr3 <- two_way_tidy_rr3%>%
  filter(FDR < 0.05)

# Save siginificant values from test to a vector to filter out for post-hoc tests
# For timepoint we care about the features which differ both between organisms and which all primary producers create or eat
significant_organism_rr3 <- as.vector(two_way_signignificant_rr3%>% filter(!anova_test == "DayNight"))$feature
significant_Timepoint_rr3 <- as.vector(two_way_signignificant_rr3%>% filter(!anova_test == "Organism"))$feature


# PRE-POST-HOC CLEANING -- Organism dunnetts ------------------------------------------------------
## Making Post-Hoc dataframes filtering for significant p_values from Two-Way Anova
rr3_organism_post_hoc <- rr3_wdf%>%
  dplyr::select(c(1:5, significant_organism_rr3))

rr3_organism_exudates_day <- rr3_organism_post_hoc%>%
  filter(DayNight == "Day")%>%
  dplyr::select(-c(Timepoint, DayNight, Experiment))%>%
  unite(sample, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_number, asin, 2:ncol(.))%>%
  spread(sample, asin)%>%
  add_column(sum = apply(.[2:ncol(.)], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(sample, asin, 2:ncol(.))%>%
  spread(feature_number, asin)%>%
  separate(sample, c("Organism", "Replicate"), sep = "_", remove = TRUE)


rr3_organism_exudates_night <- rr3_organism_post_hoc%>%
  filter(DayNight == "Night")%>%
  dplyr::select(-c(Timepoint, DayNight, Experiment))%>%
  unite(sample, c("Organism", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_number, asin, 2:ncol(.))%>%
  spread(sample, asin)%>%
  add_column(sum = apply(.[2:ncol(.)], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(sample, asin, 2:ncol(.))%>%
  spread(feature_number, asin)%>%
  separate(sample, c("Organism", "Replicate"), sep = "_", remove = TRUE)

# PRE-POST-HOC CLEANING -- Diel ANOVA ------------------------------------------------------
rr3_DayNight_post_hoc <- rr3_wdf%>%
  dplyr::select(c(1:5, significant_Timepoint_rr3))

rr3_Diel <- rr3_DayNight_post_hoc%>%
  dplyr::select(-c(Timepoint, Experiment))%>%
  unite(sample, c("Organism", "DayNight", "Replicate"), sep = "_", remove = TRUE)%>%
  gather(feature_number, asin, 2:ncol(.))%>%
  spread(sample, asin)%>%
  add_column(sum = apply(.[2:ncol(.)], 1, sum))%>%
  filter(!sum == 0)%>%
  dplyr::select(-sum)%>%
  gather(sample, asin, 2:ncol(.))%>%
  spread(feature_number, asin)%>%
  separate(sample, c("Organism", "DayNight", "Replicate"), sep = "_", remove = TRUE)

## subsetting Diel cycles to organisms
diel_poc <- rr3_Diel%>%
  filter(Organism == "Pocillopora verrucosa")

diel_por <- rr3_Diel%>%
  filter(Organism == "Porites lobata")

diel_dic <- rr3_Diel%>%
  filter(Organism == "Dictyota")

diel_trf <- rr3_Diel%>%
  filter(Organism == "Turf")

diel_cca <- rr3_Diel%>%
  filter(Organism == "CCA")

# STATS POST-HOC — Day organism dunnetts -------------------------------------------
dunnett_model_rrday <- lapply(rr3_organism_exudates_day[3:ncol(rr3_organism_exudates_day)], function(y, f)
  summary(glht(aov(y ~ f, rr3_organism_exudates_day), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(rr3_organism_exudates_day$Organism), "Water control")))

pvals_dunn_org_exudates_day <- as.data.frame(
  sapply(
    dunnett_model_rrday,
    function(x) x$test$pvalues))

pvals_dunn_org_exudates_day <- cbind("Organism" = c("CCA",
                                                    "Dictyota", 
                                                    "Pocillopora verrucosa",
                                                    "Porites lobata",
                                                    "Turf"),
                                     pvals_dunn_org_exudates_day)

rr3_day_dunnetts_FDR_exudates <- pvals_dunn_org_exudates_day%>%
  gather(feature_number, p_value, 2:ncol(.))

rr3_day_dunnetts_FDR_exudates$FDR_f <-p.adjust(rr3_day_dunnetts_FDR_exudates$p_value, method = "BH")

# STATS P-VALUE -- Day organism dunnetts ----------------------------
rr3_day_dunnett_sig_exudates <- rr3_day_dunnetts_FDR_exudates%>%
  filter(p_value < 0.1)%>%
  dplyr::select(-FDR_f)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)
# dplyr::select(-p_value)%>%
# rename(p_value = FDR_f)%>%


rr3_day_dunnett_sig_features_exudates <- as.vector(rr3_day_dunnett_sig_exudates$feature_number)


# STATS POST-HOC — Night organism dunnetts -------------------------------------------
dunnett_model_rrnight <- lapply(rr3_organism_exudates_night[3:ncol(rr3_organism_exudates_night)], function(y, f)
  summary(glht(aov(y ~ f, rr3_organism_exudates_night), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(rr3_organism_exudates_night$Organism), "Water control")))

pvals_dunn_org_exudates_night <- as.data.frame(
  sapply(
    dunnett_model_rrnight,
    function(x) x$test$pvalues))

pvals_dunn_org_exudates_night <- cbind("Organism" = c("CCA",
                                                    "Dictyota", 
                                                    "Pocillopora verrucosa",
                                                    "Porites lobata",
                                                    "Turf"),
                                     pvals_dunn_org_exudates_night)

rr3_night_dunnetts_FDR_exudates <- pvals_dunn_org_exudates_night%>%
  gather(feature_number, p_value, 2:ncol(.))

rr3_night_dunnetts_FDR_exudates$FDR_f <-p.adjust(rr3_night_dunnetts_FDR_exudates$p_value, method = "BH")

# STATS P-VALUE -- Night organism dunnetts ----------------------------
rr3_night_dunnett_sig_exudates <- rr3_night_dunnetts_FDR_exudates%>%
  filter(p_value < 0.1)%>%
  dplyr::select(-FDR_f)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)
# dplyr::select(-p_value)%>%
# rename(p_value = FDR_f)%>%


rr3_night_dunnett_sig_features_exudates <- as.vector(rr3_night_dunnett_sig_exudates$feature_number)


# STATS POST-HOC — Timepoint anovas ---------------------------------------
aov_diel_poc <- as.data.frame(sapply(diel_poc[4:ncol(diel_poc)], 
                                     function(x) summary(
                                       aov(x ~ diel_poc[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Pocillopora verrucosa", .before = 1)

aov_diel_por <- as.data.frame(sapply(diel_por[4:ncol(diel_por)], 
                                     function(x) summary(
                                       aov(x ~ diel_por[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Porites lobata", .before = 1)

aov_diel_dic <- as.data.frame(sapply(diel_dic[4:ncol(diel_dic)], 
                                     function(x) summary(
                                       aov(x ~ diel_dic[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Dictyota", .before = 1)

aov_diel_trf <- as.data.frame(sapply(diel_trf[4:ncol(diel_trf)], 
                                     function(x) summary(
                                       aov(x ~ diel_trf[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Turf", .before = 1)

aov_diel_cca <- as.data.frame(sapply(diel_cca[4:ncol(diel_cca)], 
                                     function(x) summary(
                                       aov(x ~ diel_cca[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "CCA", .before = 1)

diel_aovs <- bind_rows(aov_diel_cca, aov_diel_dic, aov_diel_poc, aov_diel_por, aov_diel_trf)

diel_aovs$FDR <- p.adjust(diel_aovs$p_value, method = "BH")

# STATS P-VALUE -- Timepoint anovas --------------------------------------------------------
diel_aov_sigs <- diel_aovs%>%
  filter(p_value < 0.1)%>%
  dplyr::select(-FDR)%>%
  unite(combined, c(Organism, feature_number), sep = "_", remove = TRUE)
# dplyr::select(-p_value)%>%
# rename(p_value = FDR_f)%>%

# METASTATS — Diel exudate differences ------------------------------------
## Fitler to exudates only
exudates_day_sig <- inner_join(rr3_day_dunnett_sig_exudates, exudate_day[1], by = "combined")

exudates_night_sig <- inner_join(rr3_night_dunnett_sig_exudates, exudate_Night[1], by = "combined")

## inner filtering to diel difference
day_only_exudates <- inner_join(exudates_day_sig, diel_aov_sigs, by = "combined", suffix = c(".dunnett", ".diel"))%>%
  add_column(diel_cycle = "Day")

night_only_exudates <- inner_join(exudates_night_sig, diel_aov_sigs, by = "combined", suffix = c(".dunnett", ".diel"))%>%
  add_column(diel_cycle = "Night")

##finding features which exude in both day night but differ between the two
day_unique_diel <- anti_join(day_only_exudates, night_only_exudates[1], by = "combined")

night_unique_diel <- anti_join(night_only_exudates, day_only_exudates[1], by = "combined")

not_unique_diel <- inner_join(day_only_exudates[c(1,2)], night_only_exudates[c(1,2)], by = "combined", suffix = c(".day", ".night"))%>%
  add_column(diel_cycle = .$p_value.dunnett.day - .$p_value.dunnett.night)%>%
  mutate(diel_cycle = case_when(diel_cycle > 0.00 ~ "preferentially day metabolite",
                                diel_cycle < 0.00 ~ "preferentially night metabolite",
                                TRUE ~ "equal?"))%>%
  add_column(p_value.dunnett = apply(.[2:3], 1, max))%>%
  dplyr::select(-c(2:3))



##combining dataframes
diel_combined <- left_join(bind_rows(day_unique_diel, night_unique_diel, not_unique_diel), rr3_DayNight_difference_RA, by = "combined")%>%
  separate(combined, c("Organism", "feature_number"), sep = "_")

diel_networking <- left_join(left_join(diel_combined, networking, by = "feature_number"), rr3_feature_table_RA, by = "feature_number")


# WRITING — Metastats -----------------------------------------------------
write_csv(diel_networking, "diel_exudates_means_rr3.csv")

