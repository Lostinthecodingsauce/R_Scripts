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

true_hits <- read_tsv("true_hits.tsv")%>%
  rename(scan_number = 'Scan')
analog_hits <- read_tsv("analog_hits.tsv")%>%
  rename(scan_number = 'Scan')
node_info <- read_csv("node_info.csv")

canopus_anotations <- read_csv("SIRIUS_etc/converted/Canopus_classes.csv")

## join library hits, analog hits, and feature table
master_df_workinprogress <- full_join(node_info, 
                                      full_join(true_hits, analog_hits, by = scan_number, suffix = c("true_hits", "analog_hits")), 
                                      by = scan_number)
## samples
ions_samples <- 20:50

## 7 different blanks
ions_blanks <- 51:58

master_df_workinprogress$max_blanks <- apply(master_df_workinprogress[ions_blanks], 1, max)

master_df_workinprogress$mean_samples <- apply(master_df_workinprogress[ion_samples], 1, mean)

## highest probability canopus annotation
canopus_best_match <- as.data.frame(canopus_anotations$name)%>%
  rename('scan_number' = 1)
canopus_best_match$canopus_probability <- apply(canopus_anotations[2:ncol(canopus_anotations)], 1, max)
canopus_best_match$canopus_annotation <- colnames(canopus_anotations)[max.col(canopus_anotations[2:1288],
                                                                              ties.method="first")]

## flagging background features and subtraction from samples ------------------------------------------------
no_background <- apply(master_df_workinprogress[mean_samples], 1, function (x) filter(x > .5*max_blanks))





