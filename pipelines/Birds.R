## Analyzing patch variability using vis_model
## written by Zach Quinlan for Ben Scott
## 10/22/2019

# READING -- libraries ------------------------------------------------------
library(tidyverse)
library(pavo)
library(phytools)
library(furrr)
library(future)


# SETTING CORES -- for parrelelizing code ---------------------------------
# sets the number of cores to one fewer than the max available
num_cores <- availableCores() -1
# actually sets the planning algoreithm for how the code will be mapped across the cores.
plan(multiprocess, workers = num_cores) 


# READING -- all CSV files in working directory ---------------------------
files <- dir(pattern = "*.csv")

reads <- files%>%
  purrr::map(read_csv)

overlap_patches <- c("back", "Mantle", "etc.")

vis_model <-reads%>%
  future_map(~ gather(., code, values, 2:ncol(.))%>%
               mutate(values = abs(values))%>%
               separate(code, c("patch", "speices_code"), sep = "_")%>%
               mutate(patch = case_when(patch == "Ventrail.tail" ~ "Ventral.tail",
                                        patch == "Undertail" ~ "Ventral.tail",
                                        patch == "PirFaries" ~ "Pirmaries",
                                        patch == "Wingbars" ~ "Wingbar",
                                        patch == "Beak.2" ~ "Back.2",
                                        patch == "Fantle" ~ "Mantle",
                                        patch == "FediuF.Coverts" ~ "Medium.Coverts",
                                        patch == "Mante.2" ~ "Mantle.2",
                                        patch == "RuFp" ~ "Rump",
                                        TRUE ~ as.character(patch)),
                      speices_code = case_when(speices_code == "GranSalF" ~ "GranSelF",
                                               TRUE ~ as.character(speices_code)))%>%
               dplyr::select(c(wl, overlap_patches))%>%
               separate(speices_code, c("species", "sex"), sep = -1)%>%
               spread(patch, values)%>%
               group_by(species, sex)%>%
               nest())%>%
  reduce(bind_rows)%>%
  mutate(data = future_map(data, ~ as.data.frame(.x)%>%
                                 vismodel(., visual = "avg.uv", achromatic = "bt.dc")))%>%
  unnest(data)

write_csv(vis_model, "quatnum.dat")
