# Functions for Metabolomics


flag_background <- function(data, blank_columns, sample_columns) {
  require("tidyverse")
  data$max_blanks <- apply(data[blank_columns], 1, max)
  data$mean_samples <- apply(data[sample_columns], 1, mean)
  
  no_background <- data%>%
    mutate(background_features = case_when(mean_samples*0.5 > max_blanks ~ "real",
                                    TRUE ~ "background"))%>%
    dplyr::select(-c(max_blanks, mean_samples))
}

flag_transient <- function(data, sample_columns, noise_level = 2E5, replication_number = 3) {
  require("tidyverse")
  no_transient <- data%>%
    add_column(samples_over_noise = rowSums(.[sample_columns] > noise_level), .before = 2)%>%
    mutate(transient_features = case_when(samples_over_noise >= replication_number ~ "real",
                                          TRUE ~ "transient"))%>%
    dplyr::select(-samples_over_noise)
}

test <- feature_table_dirty%>%
  flag_background(ions_blanks, ions_samples)%>%
  flag_transient(ions_samples)


