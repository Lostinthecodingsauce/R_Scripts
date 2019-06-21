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
  dplyr::select(-DayNight)




# Subsetting to FCM or fDOM -----------------------------------------------
fdom_wdf <- dorc_fcm_fdom%>%
  dplyr::select(4:'PARAFAC6')%>%
  filter(Timepoint == c('T0', 'TF'))

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

## ANOVA on the rate change and TF TEMPORARY. WILL RUN TOGETHER WITH FDOM DOC in TWO-WAY
summary(aov(fcm_rate_t7_t0$log10_change_per_hour ~ fcm_rate_t7_t0$Organism))
summary(aov(fcm_rate_t7_t0$log10_TF ~ fcm_rate_t7_t0$Organism))

summary(aov(fcm_rate_th_t0$log_change_per_hour ~ fcm_rate_th_t0$Organism))

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


# Overall Stats -----------------------------------------------------------
#Combine dataframes into one
stats_df <- left_join



