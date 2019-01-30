library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(magrittr)
library(readr)

crusade_df <- read.csv("Cr_OTU_tax_RA.dat")

#filter to order
data_sum_order <- crusade_df %>%
  dplyr::select(c(Taxonomy.3, Taxonomy.4), c(8:38))%>%
   # dplyr::filter(!str_detect(Taxonomy.3, "unclassified"),
   #        !str_detect(Taxonomy.3, "uncultured"),
   #        !str_detect(Taxonomy.3, "unidentified"))%>%
   dplyr::mutate(Taxonomy.4 = case_when(Taxonomy.4 %like any% c("%uncultured%", "%unclassified%", "%Subgroup%") ~ "unclassified",
                                 TRUE ~ as.character(Taxonomy.4))) %>%
  group_by(Taxonomy.3, Taxonomy.4) %>%
     summarize_if(is.numeric, sum)%>%
  ungroup(data_sum_order)
#  dplyr::filter(Taxonomy.6 == "Pseudoalteromonas")


#This next portion is simply for making pie charts. Used in CRUSADE to keep specific classes for coloring pie chart wedges
class_filter <- data_sum_order%>%
  mutate(Taxonomy.3 = case_when(Taxonomy.3 == "Alphaproteobacteria" ~ "Alphaproteobacteria",
                                  Taxonomy.3 == "Cyanobacteria" ~ "Cyanobacteria",
                                  Taxonomy.3 == "Cytophagia" ~ "Cytophagia",
                                  Taxonomy.3 == "Deltaproteobacteria" ~ "Deltaproteobacteria",
                                  Taxonomy.3 == "Flavobacteria" ~ "Flavobacteria",
                                  Taxonomy.3 == "Gammaproteobacteria" ~ "Gammaproteobacteria",
                                  Taxonomy.3 == "Planctomycetacia" ~ "Planctomycetacia",
                                  Taxonomy.3 == "Sphingobacteriia" ~ "Sphingobacteriia",
                                 TRUE ~ "Other"))

#write csv
# write.csv(transform_order, "transpose_order_T.dat")
# write.csv(transform_family, "transpose_family_T.dat")
# write.csv(data_sum_class, "transpose_class_NT.dat")
 
 