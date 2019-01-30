library(dplyr)
library(readr)
library(data.table)
library(DescTools)
library(tidyr)


crusade_df <- read_csv("data_sum_OTU.dat")

#data_sum_OTU <- crusade_df %>%
#  dplyr::select(c(Order, Family, Genus, OTU), c(8:38))



# Trying to replace All of this nonsense with a function ------------------


clean_test <- crusade_df%>%
  dplyr::mutate(Order = case_when(Order %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Order)))%>%
  dplyr::mutate(Family = case_when(Order == "unclassified" ~ "",
                                   TRUE ~ as.character(Family)))%>%
  dplyr::mutate(Genus = case_when(Order == "unclassified" ~ "",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Order == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(Family = case_when(Family %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                   TRUE ~ as.character(Family)))%>%
  dplyr::mutate(Genus = case_when(Family == "unclassified" ~ "",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Family == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(Genus = case_when(Genus %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Genus == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(OTU = case_when(OTU %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "sp",
                                TRUE ~ as.character(OTU)))


# Function Workshop -------------------------------------------------------


  
otu_OFGS <- function(data, Order, Family, Genus, OTU) {
  dplyr::mutate(data, Order = case_when(Order %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Order)))
  dplyr::mutate(data, Family = case_when(Order == "unclassified" ~ "",
                                     TRUE ~ as.character(Family)))
  dplyr::mutate(data, Genus = case_when(Order == "unclassified" ~ "",
                                    TRUE ~ as.character(Genus)))
  dplyr::mutate(data, OTU = case_when(Order == "unclassified" ~ "",
                                  TRUE ~ as.character(OTU)))
  dplyr::mutate(data, Family = case_when(Family %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                     TRUE ~ as.character(Family)))
  dplyr::mutate(data, Genus = case_when(Family == "unclassified" ~ "",
                                    TRUE ~ as.character(Genus)))
  dplyr::mutate(data, OTU = case_when(Family == "unclassified" ~ "",
                                  TRUE ~ as.character(OTU)))
  dplyr::mutate(data, Genus = case_when(Genus %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                    TRUE ~ as.character(Genus)))
  dplyr::mutate(data, OTU = case_when(Genus == "unclassified" ~ "",
                                  TRUE ~ as.character(OTU)))
  dplyr::mutate(data, OTU = case_when(OTU %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "sp",
                                  TRUE ~ as.character(OTU)))
}

otu_order <- function(data, taxon_c) {
  dplyr::select(data, -Order)
  dplyr::rename(data, "Order" = taxon_c)
}
  

cleanish <-otu_OFGS(crusade_df, crusade_df$Order, crusade_df$Family, crusade_df$Genus, crusade_df$OTU)


View(cleanish)
