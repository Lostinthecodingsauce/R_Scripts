library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(magrittr)
library(readr)

crusade_df <- read.csv("Cr_OTU_tax_RA.dat")%>%
  rename(Kingdom = 'Taxonomy.1',
         Phyla   = 'Taxonomy.2',
         Class   =  'Taxonomy.3',
         Order   =  'Taxonomy.4',
         Family  =   'Taxonomy.5',
         Genus   =   'Taxonomy.6',
         OTU     =    'Taxonomy.7')

#filter to OTU
data_sum_OTU <- crusade_df %>%
  dplyr::select(c(Order, Family, Genus, OTU), c(8:38))%>%
  dplyr::mutate(Order = case_when(Order %like any% c("%uncultured%", "%unclassified%", "%Subgroup%") ~ "unclassified",
                                  TRUE ~ as.character(Order))) %>%
  dplyr::mutate(Family = case_when(Family %like any% c("%uncultured%", "%unclassified%", "%Subgroup%") ~ "unclassified",
                                   TRUE ~ as.character(Family))) %>%
  dplyr::mutate(Genus = case_when(Genus %like any% c("%uncultured%", "%unclassified%", "%Subgroup%") ~ "unclassified",
                                  TRUE ~ as.character(Genus))) %>%
  dplyr::mutate(OTU = case_when(OTU %like any% c("%uncultured%", "%unclassified%", "%Subgroup%") ~ "unclassified",
                                TRUE ~ as.character(OTU))) %>%
  group_by(Order, Family, Genus, OTU) %>%
  summarize_if(is.numeric, sum)%>%
  ungroup(data_sum_OTU)%>%
  unite(OFGO, c(Order, Family, Genus, OTU), sep = ";", remove = TRUE)
  
transform_OTU <-one_percent_OTU%>%
  mutate_if(is.numeric, sqrt)%>% #Normalizing data (arcsin(sqrt(x)))
  mutate_if(is.numeric, asin)

#write csv
write_csv(transform_OTU, "CR_10OTU_Clean_Transform.dat")

##transposed in JMP then put back in for data wrangling
transpose_otu <- read_csv("Transpose_CR_10OTU_Clean_Transform.dat")%>%
  separate(Label, paste("C", 1:4, sep = "."), sep = "_")%>%
  rename(Organism = `C.1`,
         Water = `C.2`,
         Timepoint = `C.3`,
         DNA.Source = `C.4`)%>%
  mutate(Organism = case_when(Organism == "CC" ~ "Calcium Carbonate Control",
                              Organism == "CH" ~ "Water Control",
                              Organism == "HR" ~ "Hydrolithon reinboldii",
                              Organism == "PO" ~ "Porolithon onkodes",
                              TRUE ~ "Start Water"))%>%
  mutate(Water = case_when(Water == "1F" ~ "Filtered",
                           Water == "2F" ~ "Filtered",
                           Water == "3F" ~ "Filtered",
                           TRUE ~ "Unfiltered"))%>%
  mutate(DNA.Source = case_when(DNA.Source == "CCA" ~ "CCA Microbiome",
                                TRUE ~ "Bacterioplankton"))%>%
  unite(Name, c(DNA.Source, Organism, Water, Timepoint), sep = "      ", remove = FALSE)


write_csv(transpose_otu, "CR_10OTU_WDF.dat")


##Getting Average relative abundance for HR, PO

report_table <-transpose_otu%>%
  group_by(Organism, Water, DNA.Source)%>%
  summarise_if(is.numeric, mean)%>%
  ungroup(transpose_otu)%>%
  dplyr::filter(!Organism == "Calcium Carbonate Control",
                !Organism == "Start Water",
                !Organism == "Water Control")


write_csv(report_table, "report_table_RA.dat")
