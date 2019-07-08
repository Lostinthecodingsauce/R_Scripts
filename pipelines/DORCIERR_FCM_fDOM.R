## Script Written by Zach Quinlan 06/19/19
#  Only working on daytime exudation and remineralization

# Loading libraries -------------------------------------------------------
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



# Loading in Library ------------------------------------------------------
dorc_fcm_fdom <- read_xlsx("DORCIERR_fDOM_FCM.xlsx")%>%
  rename(sample_name =`Sample Name of DORCIERR_FCM_Final`)%>%
  rename('DayNight' = 'Experiment')%>%
  rename(Organism = 'Organsim')%>%
  filter(DayNight == "Day")%>%
  dplyr::select(-DayNight)%>%
  mutate(Organism = case_when(Organism == "Water" ~ "Water control",
                              TRUE ~ as.character(Organism)))

dorc_metab <- read_csv("DORCIERR_metabolomics_wdf.csv")

feature_metadata <- read_csv("moorea_feature_table_master_post_filtered.csv")%>%
  dplyr::select(1:68)

moorea_doc <- read_xlsx("MO17_ExpSummary_DOC.2018.04.04.xlsx")%>%
  dplyr::select(1:2)%>%
  rename(sample_name = 1)%>%
  rename(DOC = 2)
  

feature_RA <- read_csv("DORCIERR_RA_wdf.dat")


network_id <- feature_metadata[c(1,17,23, 26, 33, 38, 29)]
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

networking <-left_join(network_id, networking_energy, by = "feature_number")

# Subsetting to FCM or fDOM -----------------------------------------------
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

# Relative Abundance data cleaning ----------------------------------------
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
  spread(Organism, difference)

decrease_over_time <- day_time%>%
  filter(difference < 0)%>%
  spread(Organism, difference)


# FCM Stats prep---------------------------------------------------------------------
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


# Metabolomics Cleaning and Subsetting ------------------------------------
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

# fDOM and DOC Stats prep --------------------------------------------------------------
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
  mutate_if(is.numeric, log10)


## FCM stats for One-Way Anovas -----------------------------------------------------------
#Combine dataframes into one
fcm_stats_df <- left_join(fcm_t7, fcm_rate_th_t0[c(2:3,7)], by = c("Organism", "Replicate"))%>%
  add_column(log10_TF = log10(.$TF), .before = 5)


# Two-way Anova Stats prep ------------------------------------------------
dom_stats_wdf <- full_join(left_join(fdom_log10, doc_log10, by = "sample_name"), dorcierr_remins_cleaned, by = c("Organism", "Timepoint", "Replicate"))%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore")



# Running TWO-Way Anova DOM -------------------------------------------------
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
# For timepoint only care about ones which differ by species and timepoint not just timepoint
significant_organism_dom <- as.vector(two_way_signignificant%>% filter(anova_test == "Organism"))$feature
significant_Timepoint_dom <- as.vector(two_way_signignificant%>% filter(anova_test == "Organism*Timepoint"))$feature


## Writing the actual significance table to a folder
write_csv(two_way_signignificant, "~/Documents/SDSU/DORCIERR/Datasets/stats/Dorcierr_two_way_significant.csv")


# DOM Organism Post-Hoc Tests ------------------------------------------------------
## Making Post-Hoc dataframes filtering for significant p_values from Two-Way Anova
dom_organism_post_hoc <- dom_stats_wdf%>%
  dplyr::select(c(1:4, significant_organism_dom))

dom_organism_exudates_ph <- dom_organism_post_hoc%>%
  filter(Timepoint == "T0")%>%
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

    
dom_organism_remineralized_ph <- dom_organism_post_hoc%>%
  filter(Timepoint == "TF")%>%
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


dom_timepoint_post_hoc <- dom_stats_wdf%>%
  dplyr::select(c(1:4, significant_organism_dom))


## Running Dunnetts only both T0 and TF
set.seed(2005)


# Exudates T0 (Dunnetts) -------------------------------------------------------------
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

dom_dunnett_sig_exudates <- dom_dunnetts_FDR_exudates%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)

dunnett_sig_features_exudates <- as.vector(dom_dunnett_sig_exudates$feature_number)

dom_dunnett_exudates_networking <- right_join(networking, dom_dunnett_sig_exudates, by = "feature_number")


# Different Exudates T0 ------------------------------------------------------
dom_organism_exudates_pre_tukey <- dom_organism_exudates_ph%>%
  dplyr::select(c(1:3, dunnett_sig_features_exudates))

tukey_exudates <- sapply(dom_organism_exudates_pre_tukey[4:ncol(dom_organism_exudates_pre_tukey)], function(x)
  TukeyHSD(aov(x ~ dom_organism_exudates_pre_tukey$Organism, data = dom_organism_exudates_pre_tukey), p.adjust.methods = "BH"))

tukey_exudates_p <- as.data.frame(tukey_exudates)%>%
  rownames_to_column(var = "variable")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% '%p.adj%')%>%
  filter(value, value < 0.05)%>%
  rename(feature_number = feature_info)

tukey_exudates_p$feature_number <- tukey_exudates_p$feature_number%>%
  gsub(".dom_organism_exudates_pre_tukey.Organism.p.adj", "",.)%>%
  gsub("X","",.)

tukey_exudates_networking <- left_join(tukey_exudates_p, networking, by = "feature_number")

# Finding unique exudates to organisms ------------------------------------
dom_dunnett_exudates_networking[is.na(dom_dunnett_exudates_networking)] <- 0

exudates <- dom_dunnett_exudates_networking[17:21]

exudates[exudates > 0] <- 1

exudates <- add_column(exudates, sum = apply(exudates, 1, sum))%>%
  add_column(feature_number = dom_dunnett_exudates_networking$feature_number)%>%
  filter(.$sum == 1)%>%
  dplyr::select(-sum)

exudates$max <- colnames(exudates)[max.col(exudates[1:5])]

poc_exudates_only <- exudates%>%
  filter(max == "Pocillopora verrucosa")

por_exudates_only <- exudates%>%
  filter(max == "Porites lobata")

cca_exudates_only <- exudates%>%
  filter(max == "CCA")

dic_exudates_only <- exudates%>%
  filter(max == "Dictyota")

trf_exudates_only <- exudates%>%
  filter(max == "Turf")


# Pocillopora exudates
poc_exudates <- right_join(dom_dunnett_exudates_networking, poc_exudates_only[6], by = "feature_number")%>%
  add_column(Organism = "Pocillopora Verrucosa", .before = 1)

poc_elements <- poc_exudates[8:14]
  # summarize_if(is.numeric, mean)%>%
  # gather(element, average, 1:ncol(.))

# CCA exudates
cca_exudates <- right_join(dom_dunnett_exudates_networking, cca_exudates_only[6], by = "feature_number")%>%
  add_column(Organism = "CCA", .before = 1)

cca_elements <- cca_exudates[8:14]
  # summarize_if(is.numeric, mean)%>%
  # gather(element, average, 1:ncol(.))

# Porites exudates
por_exudates <- right_join(dom_dunnett_exudates_networking, por_exudates_only[6], by = "feature_number")%>%
  add_column(Organism = "Porites lobata", .before = 1)

por_elements <- por_exudates[8:14]
  # summarize_if(is.numeric, mean)%>%
  # gather(element, average, 1:ncol(.))

# Dictyota Exudates
dic_exudates <- right_join(dom_dunnett_exudates_networking, dic_exudates_only[6], by = "feature_number")%>%
  add_column(Organism = "Dictyota", .before = 1)

dic_elements <- dic_exudates[8:14]
  # summarize_if(is.numeric, mean)%>%
  # gather(element, average, 1:ncol(.))

# Turf Exudates
trf_exudates <- right_join(dom_dunnett_exudates_networking, trf_exudates_only[6], by = "feature_number")%>%
  add_column(Organism = "Turf", .before = 1)

trf_elements <- trf_exudates[8:14]
  # summarize_if(is.numeric, mean)%>%
  # gather(element, average, 1:ncol(.))%>%
  

exudates_all <- bind_rows(cca_exudates, por_exudates, poc_exudates, dic_exudates, trf_exudates)

## Combining tables back toghether
mean_elemntal_composition <- bind_rows(poc_elements, por_elements, cca_elements, dic_elements, trf_elements)

aov_elements <- as.data.frame(
  sapply(mean_elemntal_composition[2:8], 
         function(x) summary(
           aov(x ~ mean_elemntal_composition[["Organism"]]))[[1]][1,'Pr(>F)']))%>%
  rename(f_value = 1)%>%
  rownames_to_column(var = "Organism")

aov_elements$FDR <- p.adjust(aov_elements$f_value, method = "BH")

## Tukey elements
tukey_elements <- sapply(mean_elemntal_composition[c(2,4,6,7)], function(x)
  TukeyHSD(aov(x ~ mean_elemntal_composition$Organism, data = mean_elemntal_composition), p.adjust.methods = "BH"))

p_values_tukey_elements <- as.data.frame(tukey_elements)%>%
  rownames_to_column(var = "variable")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% '%p.adj%')%>%
  filter(value, value < 0.05)

p_values_tukey_elements$feature_info <- gsub(".mean_elemntal_composition.Organism.p.adj", "",
                                             p_values_tukey_elements$feature_info)
## Graphing
graphing_elements <- mean_elemntal_composition%>%
  gather(element, number, c(2,4,6,7))

ggplot(graphing_elements, aes(x = Organism, y = number, fill = element))+
  geom_bar(stat = "summary", fun.y = "mean")+
  scale_color_manual(values=wes_palette(n=4, name="Darjeeling1"))
  

# Remineralized TF DOM ----------------------------------------------------
dom_org_dunnetts_remins <- lapply(dom_organism_remineralized_ph[4:ncol(dom_organism_remineralized_ph)], function(y, f)
  summary(glht(aov(y ~ f, dom_organism_remineralized_ph), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(dom_organism_remineralized_ph$Organism), "Water control")))

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

dom_dunnett_sig_remins <- dom_dunnetts_FDR_remins%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)


dunnett_sig_features_remin <- as.vector(dom_dunnett_sig_remins$feature_number)

dom_dunnett_remins_networking <- right_join(networking, dom_dunnett_sig_remins, by = "feature_number")


# writing significantly exuded compounds ----------------------------------
sig_compounds_wdf <- dom_organism_post_hoc%>%
  dplyr::select(c(1:3, dunnett_sig_features_exudates, dunnett_sig_features_remin))

write_csv(sig_compounds_wdf, "DORCIERR_Sig_compounds_dunnetts.dat")


# DOM Timepoint Post-hoc --------------------------------------------------
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

aov_time_poc <- as.data.frame(sapply(timepoint_poc[4:ncol(timepoint_poc)], 
                                     function(x) summary(
                                       aov(x ~ timepoint_poc[["Timepoint"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Pocillopora", .before = 1)

aov_time_por <- as.data.frame(sapply(timepoint_por[4:ncol(timepoint_por)], 
                       function(x) summary(
                         aov(x ~ timepoint_por[["Timepoint"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)%>%
  add_column(Organism = "Porites", .before = 1)

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

time_aov_sigs <- time_aovs%>%
  filter(FDR < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR)

time_aov_sig_features <- as.vector(time_aov_sigs$feature_number)

# labile compounds --------------------------------------------------------
## These are defined as features which are present in T0 for any organism and then decrease from T0 -> TF

time_differ_exudate_filtered <- dom_dunnett_exudates_networking%>%
  filter(feature_number %like any% time_aov_sig_features)

labile_exudate_filtered_RA <- inner_join(time_aov_sigs, 
                                         inner_join(time_differ_exudate_filtered, decrease_over_time, 
                                             by = "feature_number", suffix = c(".p_exudation", ".RA")),
                                             by = "feature_number", suffix = c(".p_timepoint", ""))


## specific groupings
#Single organisms
poc_labile <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` != 0,
         `Porites lobata.p_exudation` == 0,
         `Dictyota.p_exudation` == 0,
         `CCA.p_exudation` == 0,
         `Turf.p_exudation` == 0,
         Pocillopora != "NA")%>%
  add_column(tag = "poc_labile")

por_labile <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` == 0,
         `Porites lobata.p_exudation` != 0,
         `Dictyota.p_exudation` == 0,
         `CCA.p_exudation` == 0,
         `Turf.p_exudation` == 0,
         Porites != "NA")%>%
  add_column(tag = "por_labile")


cca_labile <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` == 0,
         `Porites lobata.p_exudation` == 0,
         `Dictyota.p_exudation` == 0,
         `CCA.p_exudation` != 0,
         `Turf.p_exudation` == 0,
         CCA != "NA")%>%
  add_column(tag = "cca_labile")


dic_labile <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` == 0,
         `Porites lobata.p_exudation` == 0,
         `Dictyota.p_exudation` != 0,
         `CCA.p_exudation` == 0,
         `Turf.p_exudation` == 0,
         Dictyota != "NA")%>%
  add_column(tag = "dic_labile")


trf_labile <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` == 0,
         `Porites lobata.p_exudation` == 0,
         `Dictyota.p_exudation` == 0,
         `CCA.p_exudation` == 0,
         `Turf.p_exudation` != 0,
         Turf != "NA")%>%
  add_column(tag = "trf_labile")



#Classes of organisms
coral_labile <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` != 0,
         `Porites lobata.p_exudation` != 0,
         `Dictyota.p_exudation` == 0,
         `CCA.p_exudation` == 0,
         `Turf.p_exudation` == 0,
         Pocillopora != "NA",
         Porites != "NA")

algae_labile <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` == 0,
         `Porites lobata.p_exudation` == 0,
         `Dictyota.p_exudation` != 0,
         `CCA.p_exudation` != 0,
         `Turf.p_exudation` != 0,
         CCA != "NA",
         Dictyota != "NA",
         Turf != "NA")

fleshy_labile <- labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` == 0,
         `Porites lobata.p_exudation` == 0,
         `Dictyota.p_exudation` != 0,
         `CCA.p_exudation` == 0,
         `Turf.p_exudation` != 0,
         Dictyota != "NA",
         Turf != "NA")

primary_labile <-labile_exudate_filtered_RA%>%
  filter(`Pocillopora verrucosa.p_exudation` != 0,
         `Porites lobata.p_exudation` != 0,
         `Dictyota.p_exudation` != 0,
         `CCA.p_exudation` != 0,
         `Turf.p_exudation` != 0,
         CCA != "NA",
         Dictyota != "NA",
         Turf != "NA",
         Pocillopora != "NA",
         Porites != "NA")

##Calculating bond energies of individual labile compounds
bond_energies <- bind_rows(poc_labile, por_labile, dic_labile, cca_labile, trf_labile)%>%
  dplyr::select(c(33, 21))

bond_tukey <- TukeyHSD(aov(
  bond_energies$cox_gibbs_energy ~ bond_energies$tag, data = bond_energies), 
  p.adjust.methods = "BH")

bond_pvalues <- as.data.frame(bond_tukey$`bond_energies$tag`)%>%
  rownames_to_column(var = "organism")


#combined labile compounds  
combined_labile_compounds <- bind_rows(poc_labile, por_labile, dic_labile, cca_labile, trf_labile)

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
  summarize_if(is.numeric, sum)

write_csv(canopus_no_networks, "./staring_at_data/canopus_no_networks.dat")


# Recalcitrant compounds -------------
## Binning compounds different in remins by produced compounds 
time_differ_produced_filtered <- dom_dunnett_remins_networking%>%
  filter(feature_number %like any% time_aov_sig_features)

recalcitrant_filtered_RA <- inner_join(time_aov_sigs, 
                                         inner_join(time_differ_produced_filtered, increase_over_time, 
                                                    by = "feature_number", suffix = c(".p_produced", ".RA")),
                                         by = "feature_number", suffix = c(".p_timepoint", ""))


# One Way Anova FCM -------------------------------------------------------
aov_fcm <-as.data.frame(
  sapply(fcm_stats_df[c(4,6)], 
         function(x) summary(
           aov(x ~ fcm_stats_df[["Organism"]]))[[1]][1,'Pr(>F)']))%>%
  rename(f_value = 1)%>%
  rownames_to_column(var = "test")

aov_fcm$FDR <- p.adjust(aov_fcm$f_value, method = "BH")

# FCM Post-Hoc Tests ----------------------------------------------------------
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


## Visualizations to help explain data --------------
ggplot(fcm_rate_t7_t0, aes(x = Organism, y = log10_change_per_hour))+
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)

fcm_graphing_wdf <- fcm_wdf%>%
  filter(!Organism == "Influent",
         !Organism == "Offshore")%>%
  add_column(log10_cells = log10(.$`Cells µL-1`))%>%
  dplyr::mutate(Timepoint = case_when(Timepoint == "T0" ~ 0,
                                      Timepoint == "T1" ~ 4,
                                      Timepoint == "T2" ~ 16,
                                      Timepoint == "T3" ~ 19,
                                      Timepoint == "T4" ~ 24,
                                      Timepoint == "T5" ~ 30,
                                      Timepoint == "T6" ~ 39,
                                      Timepoint == "TF" ~ 47,
                                      TRUE ~ as.numeric(Timepoint)))

ggplot(fcm_graphing_wdf, aes(x = Timepoint, y = log10_cells, col = Organism))+
  geom_point(stat = "summary", fun.y = "mean")+
  stat_summary(fun.y=mean, geom="line")

ggplot(bond_energies, aes(x = tag, y = bond_energy))+
  geom_point(stat = "summary", fun.y = "mean")


# PCoA --------------------------------------------------------------------

dom <- as.matrix(canopus_no_networks%>%
                   dplyr::select(-1)) ## No fDOM because it is negative

veg_bray_dom <- vegdist(dom, "bray") #Bray-curtis distances

pc_scores <- pcoa(veg_bray_dom) #Calculating scores

#This plots Eigenvalues
#They will allow you to choose the best axes to show how your data varies
ylimit = c(0, 1.1*max(pc_scores$values$Relative_eig))

Eigan <- barplot(pc_scores$values$Relative_eig[1:10], ylim= ylimit)
# Add values to the bars
text(x = Eigan, y = pc_scores$values$Relative_eig[1:10], label = pc_scores$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

# S3 method for pcoa
biplot(pc_scores, Y=NULL, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1)

pco_scores <- as.data.frame(pc_scores$vectors)
pco_scores$organism <- canopus_no_networks$canopus_annotation     # This will add reference labels to the PCoA scores

ggplot(pco_scores, mapping = aes(x = Axis.1, y = Axis.2))


write.csv(pco_scores, "PCoA_scores_labile_exudates.dat")

