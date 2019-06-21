## Script written for Colton Johnson by Zach Quinlan on June 18th 2019

library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(ape)
library(wesanderson)
library(readxl)


# Loading in ALL data frames ----------------------------------------------
master_functionalgroups_s1_raw <- read_excel("Master_FG_lvl1_Data_Raw.xlsx")
master_metadata <- read_excel("Master_FG_lvl1_Metadata_Raw.xlsx")%>%
  rename(sample_name = Sample)

master_bacteria_clades <- read_csv("All_Taxa_Leopard_Master_Redone_Aquarium.csv")
  

# Data Cleaning Functional Groups-----------------------------------------------------------
##  This is transposing your data. So now functional genes will be columns and samples will be rows
ls_tranfg <-master_functionalgroups_s1_raw %>%
  gather(sample_name, RA, 2:ncol(.))%>%
  spread('Subsystem 1', RA)%>%
  add_column(TR = apply(.[2:ncol(.)], 1, sum))


# Relative abundance
feature_relative_abundance <- 
  as.data.frame(
    sapply(
      ls_tranfg[2:ncol(ls_tranfg)],
      function(x) x/ls_tranfg$TR))%>%
  add_column(sample_name = ls_tranfg$sample_name, .before = 1)

# angular transformation of RA
feature_asin_sqrt <-
  as.data.frame(
    sapply(
      feature_relative_abundance[2:ncol(ls_tranfg)],
      function(x) asin(sqrt(x))))%>%
  add_column(sample_name = ls_tranfg$sample_name)


#Joining RA and angular transformed data for 
RA_asin_df <- left_join(feature_relative_abundance, feature_asin_sqrt, by = 'sample_name',
                        suffix = c(".RA", ".asin"))%>%
  dplyr::select(-c("TR.RA", "TR.asin"))

## Combining with Meta Data and filtering out the two bad samples to make a working data frame(wdf)
# This dataframe will have all of the samples for all types of sharks in all locations at subsystem 1
master_functional_s1_wdf <- left_join(master_metadata, RA_asin_df, by = "sample_name")%>%
  filter(!sample_name == "LS 2 2015 LJ")

# Subsetting Data Frames Functional Groups --------------------------------------------------
la_jolla_sharks_s1 <- master_functional_s1_wdf%>%
  filter(Location == "La Jolla",
         Group == "Leopard Shark")

ljs_meta <- 1:5
ljs_RA <-6:40
ljs_asin <- 41:75

# PERMANOVA  Functional Groups---------------------------------------------------------------
#Pairwise
pairwiseAdonis::pairwise.adonis(la_jolla_sharks_s1[ljs_asin], la_jolla_sharks$Year, p.adjust.m = "BH")

#Non-pairwise
adonis(la_jolla_sharks_s1[ljs_asin] ~ Year, la_jolla_sharks_s1, p.adjust.m = "BH")


# SIMPER  Functional Groups------------------------------------------------------------------
simper_ljs_s1 <- simper(la_jolla_sharks_s1[ljs_asin], la_jolla_sharks_s1$Year, permutations = 0, trace = FALSE)

summary(simper_ljs_s1)

# ANOVA  Functional Groups-------------------------------------------------------------------
p_values_oneway_sharks_functional <- 
  as.data.frame(
    sapply(
      la_jolla_sharks_s1[ljs_asin],
      function(x) summary(aov(x ~ la_jolla_sharks_s1[["Year"]]))[[1]][1,'Pr(>F)']))

sig_ljs_s1 <- p_values_oneway_sharks_functional%>%
  rownames_to_column(var = "Subsystem_1")%>%
  rename('p_value' = 2)%>%
  filter(p_value, p_value < 0.05)

sig_features_ljs_s1 <- as.vector(sig_ljs_s1$Subsystem_1)


# Tukeys  Functional Groups------------------------------------------------------------------
ljs_s1_tukey_cleaning <- la_jolla_sharks_s1%>%
  dplyr::select(c(sig_features_ljs_s1))%>%
  add_column(Year = as.factor(la_jolla_sharks_s1$Year), .before = 1)
  
tukey_model_ljs_s1<- sapply(ljs_s1_tukey_cleaning[2:ncol(ljs_s1_tukey_cleaning)], function(x)
  TukeyHSD(aov(x ~ ljs_s1_tukey_cleaning$Year, data = ljs_s1_tukey_cleaning), p.adjust.methods = "BH"))

p_values_tukey_functional_sharks <- as.data.frame(tukey_model_ljs_s1)%>%
  rownames_to_column(var = "subsystem_1")%>%
  gather(gene_info, value, 2:ncol(.))%>%
  filter(gene_info %like% '%p.adj%')%>%
  filter(value, value < 0.05)

p_values_tukey_functional_sharks$gene_info <- p_values_tukey_functional_sharks$gene_info%>%
  gsub(".ljs_s1_tukey_cleaning.Year.p.adj", "", .)


# PCoA Functional Groups--------------------------------------------------------------------
veg_bray <- vegdist(la_jolla_sharks_s1[ljs_asin], "bray") #Bray-curtis distances

pc_scores<-pcoa(veg_bray) #Calculating scores

#This plots Eigenvalues
#They will allow you to choose the best axes to show how your data varies
ylimit = c(0, 1.1*max(pc_scores$values$Relative_eig))

Eigan <- barplot(pc_scores$values$Relative_eig[1:10], ylim= ylimit)
# Add values to the bars
text(x = Eigan, y = pc_scores$values$Relative_eig[1:10], label = pc_scores$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

# S3 method for pcoa
biplot(pc_scores, Y=NULL, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1)

pco_scores_functional <- as.data.frame(pc_scores$vectors)
pco_scores_functional$Year <- as.factor(la_jolla_sharks_s1$Year)

ggplot(pco_scores_functional, mapping = aes(Axis.1, Axis.2, col = pco_scores_functional$Year)) +
  geom_point(stat = "identity") +
  scale_color_manual(values=wes_palette(n=4, name="Darjeeling1")) +
  xlab("PCoA 1") + 
  ylab("PCoA 2") +
  title("Sharks in La Jolla")


# Writing Functional Group Dataframes to csv -----------------------------------------------
write_csv(master_functional_s1_wdf, "Master_functionalgroups_RA_asin_subsystem1.csv")

write_csv(la_jolla_sharks_s1, "LaJolla_sharks_master_subsystem1.csv")

write_csv(p_values_tukey_functional_sharks, "LaJolla_sharks_tukey_pvalues.csv")

write_csv(simper_ljs_s1, "LaJolla_sharks_SimperValues.csv")


####### Cleaning Dataframes for Taxonomy ########----------------------------------------
ls_genera <- master_bacteria_clades%>%
  gather(sample_name, RA, 2:ncol(.))%>%
  spread(1, RA)

# angular transformation of RA THIS PART IS FUCKED CAUSE FOCUS DOES SOME TRANSFORMATION BUT I DO NOT KNOW WHAT IT IS.
# COLTON NEEDS TO CHECK THIS FOR ME.
genera_asin_sqrt <-
  as.data.frame(
    sapply(
      ls_genera[2:ncol(ls_genera)],
      function(x) asin(sqrt(x))))%>%
  add_column(sample_name = ls_genera$sample_name)

## Joining RA and Asin together
genera_RA_asin <- left_join(ls_genera, genera_asin_sqrt, by = "sample_name", suffix = c(".RA", ".asin"))



