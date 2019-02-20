rm(list = ls())


library(tidyverse)
library(data.table)


analog_hits <- read_tsv("Analog_hits.tsv")%>%
  rename(Scan_Number = 'Scan_analog')
true_hits <-read_tsv("library_hits_true.tsv")%>%
  rename(Scan_Number = 'Scan')
old_hits <- read_csv("Moorea_SPEDOM_MSMS_Polished_April2018.csv")%>%
  rename(Scan_Number = 'row ID')
Node_info <- read_tsv("Node_Info.clustersummary.tsv")

new_run <- full_join(full_join(Node_info, true_hits, "Scan_Number"),analog_hits, "Scan_Number")

both_runs <- full_join(new_run, old_hits, "Scan_Number")

both_runs_by_old <- right_join(new_run, old_hits, "Scan_Number")


mass_check <- both_runs%>%
  group_by(Scan_Number)

scan_IDs <- as.data.frame(old_hits$Scan_Number)

write_csv(both_runs, "Moorea_2017_new_and_old_runs.csv")
write_delim(scan_IDs, "Feature_ID.txt", delim = " ", col_names = FALSE)

write_csv(both_runs_by_old, "new_by_old.csv")


samples <- both_runs_by_old%>%
  select(c(Feature_Name:D_DT_2_T0N_C18))%>%
  gather(sample_ID, RA, 2:318)%>%
  spread(Feature_Name, RA)

moorea_wdf <-samples%>%
  separate(sample_ID, paste("c", 1:5, sep = "."), sep = "_")%>%
  rename(Experiment = `c.1`,
         Organism = `c.2`,
         Replicate = `c.3`,
         Timepoint = `c.4`,
         DOM_source = `c.5`)%>%
  filter(!Experiment %like% "Blank")%>%
  dplyr::mutate(Experiment = case_when(Experiment == "D" ~ "dorcierr",
                                       Experiment == "M" ~ "mordor",
                                       Experiment == "R" ~ "RR3",
                                       TRUE ~ as.character(Experiment)))%>%
  dplyr::mutate(Organism = case_when(Organism == "D"))


moorea_subset <- moorea_wdf%>%
  select(1:6)
Experiment_codes <-moorea_subset%>%
  group_by(Experiment)%>%
  summarise_if(is.numeric, sum)
Organism_codes <- moorea_subset%>%
  group_by(Organism)%>%
  summarise_if(is.numeric, sum)
Timepoint_codes <- moorea_subset%>%
  group_by(Timepoint)%>%
  summarise_if(is.numeric, sum)
DOM_source_codes <- moorea_subset%>%
  group_by(DOM_source)%>%
  summarise_if(is.numeric, sum)

