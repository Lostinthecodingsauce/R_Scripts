## Plan-c analysis pipeline from Mo'orea 2019 thermal bleaching experiment
## Script created Oct. 24, 2019 - Zach Quinlan



# READING -- Libraries ----------------------------------------------------
#Data mungering
library(tidyverse)
library(data.table)
library(DescTools)
library(broom)
library(readxl)
library(multcomp)
library(CHNOSZ)
library(furrr)
library(future)
library(readxl)

#PCoA, PERMANOVA
library(vegan)
library(ape)
library(wesanderson)
library(RColorBrewer)

# SETTING CORES -- for parrelelizing code ---------------------------------
# sets the number of cores to one fewer than the max available
num_cores <- availableCores() -1
# actually sets the planning algoreithm for how the code will be mapped across the cores.
plan(multiprocess, workers = num_cores) 


# READING -- dataframes ---------------------------------------------------
# FCM data
fcm_files_raw <- dir(path = "~/Documents/SDSU/plan_c/FCM/", pattern = "*.csv")%>%
  map(read_csv)%>%
  map(~ dplyr::select(., c("Plate Name", "Sample Name", "Gate Name", "Event Count")))%>%
  reduce(bind_rows)%>%
  rename("plate" = "Plate Name",
         "Well" = "Sample Name")

fcm_names_raw <- dir(path = "~/Documents/SDSU/plan_c/FCM/", pattern = "*.xlsx")

fcm_names <- fcm_names_raw%>%
  map(read_xlsx)%>%
  set_names(fcm_names_raw)

## Surface area
surface_area_raw <- read_csv("~/Documents/SDSU/plan_c/Nubbin_sizing_notes.csv")

# CLEANING -- FCM ---------------------------------------------------------
fcm_raw <- fcm_names%>%
  map_df( ~as.data.frame(.), .id="plate")%>%
  mutate(plate = gsub(".Map.xlsx", "", plate),
         plate = gsub("OCT.2019", "OCT.18.2019", plate))%>%
  full_join(., fcm_files_raw, by = c("plate", "Well"))%>%
  separate(plate, c("experiment_plate", "dye", "month", "day", "year", "plate", "number"))%>%
  unite(date, c("month", "day", "year"), sep = "-")%>%
  unite(plate_numnber, c("plate", "number"), sep = "_")%>%
  separate(Sample, c("Experiment", "health", "Timepoint"), sep = "_")%>%
  separate(health, c("health", "replicate"), sep = 2)%>%
  separate(Timepoint, c("timepoint", "bath"), sep = 2)%>%
  mutate(health = case_when(health == "HE" ~ "healthy",
                            health == "PB" ~ "partially bleached",
                            health == "BL" ~ "bleached",
                            health == "WA" ~ "water control",
                            TRUE ~ as.character(health)),
         timepoint = gsub("T", "", timepoint),
         bath = case_when(bath == "C" ~ "control",
                          bath == "H" ~ "hot",
                          TRUE ~ as.character(bath)))

fcm_sybr <- fcm_raw%>%
  filter(Experiment == "PC",
         dye == "SYBR",
         `Gate Name` == "SYBR Polygon")%>%
  mutate(cells_µL = `Event Count`/0.75)%>%
  group_by(timepoint, health, bath)%>%
  mutate(std = sd(cells_µL))%>%
  ungroup()


# GRAPHING -- FCM ---------------------------------------------------------
ggplot(fcm_sybr, aes(x = timepoint, y = cells_µL, color = health))+
  geom_point(stat = "summary", fun.y = "mean") + 
  geom_line(aes(group = health), stat = "summary", fun.y = "mean") +
  facet_wrap(~bath)






