library(dplyr)

cr_df<- read.csv("Old/Transpose_of_master_cleaned.csv")

cr_filt_rhodo <- cr_df %>%
  mutate()%>%
  dplyr::filter(genus =="Rhodobacter")

write.csv(cr_filt_rhodo, file ="crusade_Rhodo.csv")
