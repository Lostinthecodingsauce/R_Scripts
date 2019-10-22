# Written by Zach Quinlan for Irina Koester for correlating microbe abundance to feature abundance
# written on 10/21/2019

# Read in libraries and dataframes ----------------------------------------

library(tidyverse)
library(d3heatmap)
library(Hmisc)
library(htmlwidgets)
library(pheatmap)
library(reshape)
library(RColorBrewer)


# Reading in dataframes ---------------------------------------------------
data <- read_csv("~/Downloads/Correlation_Pn_Ex2_t_Unfil_input_average.csv")


# Actually running the correlation test -----------------------------------

transpose <- data%>%
  gather(sample_name, ra, 2:ncol(.))%>%
  spread(SampleID, ra)%>%
  gather(otu, otu_ra, 4434:ncol(.))%>%
  gather(feature_number, feature_ra, 2:4433)%>%
  group_by(feature_number, otu)%>%
  nest()%>%
  mutate(corr = map(data, ~ cor.test(.x$feature_ra, .x$otu_ra, method = "pearson")%>%
                      broom::tidy()))

correlation_pvalues <- transpose%>%
  dplyr::select(-data)%>%
  unnest(corr)%>%
  mutate(FDR = p.adjust(p.value, method = "BH"))


