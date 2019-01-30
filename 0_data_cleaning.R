# This script cleans up Zach's monstrosity of data and work flow
# 2/1/2018

#load libraries
library(dplyr)
library(stringr)
library(data.table)

# Load Data

data_raw <- read_csv("Cr_OTU_tax_RA.dat") 

# Rename columns
data_filt_family <- data_raw %>% 
  rename(kingdom = `Taxonomy 1`,
         phylum  = `Taxonomy 2`,
         class   = `Taxonomy 3`,
         order   = `Taxonomy 4`,
         family  = `Taxonomy 5`,
         genus   = `Taxonomy 6`,
         species = `Taxonomy 7`)

# group all the data unclassified at genus level into to one observation 
data_sum_genus <- data_filt_family %>%
  mutate(kingdom = case_when(kingdom %like% "uncultured"   ~ "",
                           kingdom %like% "unclassified" ~ "",
                           kingdom %like% "unidentified" ~ "",
                           TRUE                        ~ kingdom)) %>%
  mutate(phylum = case_when(phylum %like% "uncultured"   ~ "",
                           phylum %like% "unclassified" ~ "",
                           phylum %like% "unidentified" ~ "",
                           TRUE                        ~ phylum)) %>%
  mutate(class = case_when(class %like% "uncultured"   ~ "",
                           class %like% "unclassified" ~ "",
                           class %like% "unidentified" ~ "",
                           TRUE                        ~ class)) %>%
  mutate(family = case_when(family %like% "uncultured"   ~ "",
                            family %like% "unclassified" ~ "",
                            family %like% "unidentified" ~ "",
                            TRUE                        ~ family)) %>%
  mutate(order = case_when(order %like% "uncultured"   ~ "",
                           order %like% "unclassified" ~ "",
                           order %like% "unidentified" ~ "",
                           TRUE                        ~ order)) %>%
  mutate(genus = case_when(genus %like% "uncultured"   ~ "",
                           genus %like% "unclassified" ~ "",
                           genus %like% "unidentified" ~ "",
                           TRUE                        ~ genus)) %>% 
  mutate(species = case_when(species %like% "uncultured"   ~ "",
                             species %like% "unclassified" ~ "",
                             species %like% "unidentified" ~ "",
                             TRUE                          ~ species)) %>%
  mutate(species = case_when(species == "" & genus != ""  ~ "sp",
                             TRUE                          ~ species)) %>% 
  group_by_if(is.character) %>% 
  summarise_if(is.numeric, sum) 
 
# write data to csv

write_csv(data_sum_genus, "Crusade_RA_cleaned.dat")
