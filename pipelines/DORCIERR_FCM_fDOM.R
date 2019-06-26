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


networking <- feature_metadata[c(1,17,23, 26, 33, 38, 29)]
networking$feature_number <- as.character(networking$feature_number)

# Subsetting to FCM or fDOM -----------------------------------------------
fdom_wdf <- dorc_fcm_fdom%>%
  dplyr::select(c(4:'PARAFAC6'))%>%
  filter(!Timepoint == "T1",
         !Timepoint == "T2",
         !Timepoint == "T3",
         !Timepoint == "T4",
         !Timepoint == "T5",
         !Timepoint == "T6")

fcm_wdf <- dorc_fcm_fdom%>%
  dplyr::select(c(1:7,33:ncol(.)))


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

# fDOM Stats prep --------------------------------------------------------------
fdom_log10 <-
  as.data.frame(
    sapply(
      fdom_wdf[14:20],
      function(x) log10(x)))%>%
  add_column(sample_name = fdom_wdf$sample_name, .before = 1)%>%
  add_column(Replicate = fdom_wdf$Replicate, .after = 1)%>%
  add_column(Timepoint = fdom_wdf$Timepoint, .after = 1)%>%
  add_column(Organism = fdom_wdf$Organism, .after = 1)


## FCM stats for One-Way Anovas -----------------------------------------------------------
#Combine dataframes into one
fcm_stats_df <- left_join(fcm_t7, fcm_rate_th_t0[c(2:3,7)], by = c("Organism", "Replicate"))%>%
  add_column(log10_TF = log10(.$TF), .before = 5)


# Two-way Anova Stats prep ------------------------------------------------
dom_stats_wdf <- full_join(fdom_log10, dorcierr_remins_cleaned, by = c("Organism", "Timepoint", "Replicate"))%>%
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

#Save siginificant values from test to a vector to filter out for post-hoc tests
significant_organism_dom <- as.vector(two_way_signignificant%>% filter(anova_test == "Organism"))$feature
significant_Timepoint_dom <- as.vector(two_way_signignificant%>% filter(anova_test == "Timepoint"))$feature


## Writing the actual significance table to a folder
write_csv(two_way_signignificant, "~/Documents/SDSU/DORCIERR/Datasets/stats/Dorcierr_two_way_significant.csv")


# DOM Post-Hoc Tests ------------------------------------------------------
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

#Exudates T0
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
                                                      "Porites Lobata",
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

# Remineralized TF
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
                                                    "Porites Lobata",
                                                    "Turf"),
                                     pvals_dunn_org_remins_dom)

dom_dunnetts_FDR_remins <- pvals_dunn_org_remins_dom%>%
  gather(feature_number, p_value, 2:ncol(.))

dom_dunnetts_FDR_remins$FDR_f <-p.adjust(dom_dunnetts_FDR_remins$p_value, method = "BH")

dom_dunnett_sig_remins <- dom_dunnetts_FDR_remins%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)%>%


dunnett_sig_features_remin <- as.vector(dom_dunnett_sig_remins$feature_number)

dom_dunnett_remins_networking <- right_join(networking, dom_dunnett_sig_remins, by = "feature_number")


# FCM Post-Hoc Tests ----------------------------------------------------------
#Tukey on the rate change
tukey_model_fcm <- 
  TukeyHSD(
    aov(
      fcm_rate_t7_t0$log10_change_per_hour ~ fcm_rate_t7_t0$Organism, data = fcm_rate_t7_t0), 
    p.adjust.methods = "BH")

p_values_tukey_fcm <- as.data.frame(tukey_model_fcm$`fcm_rate_t7_t0$Organism`)%>%
  rownames_to_column(var = "Organism")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% "%p adj%")%>%
  filter(value < 0.05)

# growth rates for the first half of the expierment
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

## Visualizations to help explain data
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




