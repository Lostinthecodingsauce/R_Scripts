## RR3 Stats
#Read it in
rr3_features <- read_tsv("AnnotatedRR3FeaturesForZach.txt")

moorea_binned <- moorea_post_stats[c(1, 948:963)]
  
rr3_features_check <- as.vector(rr3_features$feature_number)

moorea_binned_check <- inner_join(moorea_binned, rr3_features, by = "feature_number")

rr3_dunnetts <- left_join(rr3_features, moorea_binned, by = "feature_number")


exuded_day <- rr3_wdf%>%
  dplyr::select(c(2,3,5, dunnett_rr3_sig_features))%>%
  unite(org_rep, c("Organism", "Replicate", "DayNight"), sep = "_", remove = TRUE)%>%
  gather(feature_name, asin, 2:ncol(.))%>%
  spread(org_rep, asin)%>%
  gather(org_rep, asin, 2:ncol(.))%>%
  spread(feature_name, asin)%>%
  separate(org_rep, c("Organism", "Replicate", "DayNight"), sep = "_", remove = TRUE)


water_dvn <-exuded_day%>%
  filter(Organism == "Water control")

porites_dvn <- exuded_day%>%
  filter(Organism == "Porites lobata")

pocil_dvn <- exuded_day%>%
  filter(Organism == "Pocillopora verrucosa")

cca_dvn <- exuded_day%>%
  filter(Organism == "CCA")

turf_dvn <- exuded_day%>%
  filter(Organism == "Turf")

dictyota_dvn <- exuded_day%>%
  filter(Organism == "Dictyota")


water_day_v_night <- 
  as.data.frame(
    sapply(
      water_dvn[4:ncol(water_dvn)],
      function(x) summary(aov(x ~ water_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)

water_day_v_night$FDR_water <- p.adjust(water_day_v_night$p_value, method = "BH")


porites_day_v_night <- as.data.frame(sapply(porites_dvn[4:ncol(porites_dvn)],
      function(x) summary(aov(x ~ porites_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)

porites_day_v_night$FDR_porites <- p.adjust(porites_day_v_night$p_value, method = "BH")

pocil_day_v_night <- 
  as.data.frame(
    sapply(
      pocil_dvn[4:ncol(pocil_dvn)],
      function(x) summary(aov(x ~ pocil_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)

pocil_day_v_night$FDR_pocillopora <- p.adjust(pocil_day_v_night$p_value, method = "BH")

cca_day_v_night <- 
  as.data.frame(
    sapply(
      cca_dvn[4:ncol(cca_dvn)],
      function(x) summary(aov(x ~ cca_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)

cca_day_v_night$FDR_cca <- p.adjust(cca_day_v_night$p_value, method = "BH")

turf_day_v_night <- 
  as.data.frame(
    sapply(
      turf_dvn[4:ncol(turf_dvn)],
      function(x) summary(aov(x ~ turf_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value = 2)

turf_day_v_night$FDR_turf <- p.adjust(turf_day_v_night$p_value, method = "BH")

dictyota_day_v_night <- 
  as.data.frame(
    sapply(
      dictyota_dvn[4:ncol(dictyota_dvn)],
      function(x) summary(aov(x ~ dictyota_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value_dictyota = 2)

dictyota_day_v_night$FDR_dictyota <- p.adjust(dictyota_day_v_night$p_value, method = "BH")


p_values_combined <- 
  full_join(
    full_join(
      full_join(
        full_join(
          full_join(water_day_v_night, cca_day_v_night, by = "feature_number"
          ), porites_day_v_night, by = "feature_number"
        ), pocil_day_v_night, by = "feature_number"
      ), turf_day_v_night, by = "feature_number"
    ), dictyota_day_v_night, by = "feature_number"
  )


# 06/25/2019 --------------------------------------------------------------

rr3_new <- read_xlsx("MappedRRMetabolitesForCraig.xlsx")

feature_asin_sqrt <-
  as.data.frame(
    sapply(
      rr3_new[5:ncol(rr3_new)],
      function(x) asin(sqrt(x))))%>%
  add_column(level_4 = rr3_new$`Level 4`, .before = 1)

rr3_aov <- feature_asin_sqrt%>%
  gather(Organism, asin, 2:ncol(.))%>%
  spread(level_4, asin)

rr3_aov$Organism <- rr3_aov$Organism%>% 
  gsub("__1", "", .)%>%
  gsub("__2", "", .)%>%
  gsub("__3", "", .)%>%
  gsub("2", "", .)%>%
  gsub("3", "", .)


rr3_ra <- rr3_new%>%
  dplyr::select(-c(1:3))%>%
  gather(Organism, asin, 2:ncol(.))%>%
  spread(`Level 4`, asin)

rr3_ra$Organism <- rr3_ra$Organism%>% 
  gsub("__1", "", .)%>%
  gsub("__2", "", .)%>%
  gsub("__3", "", .)%>%
  gsub("2", "", .)%>%
  gsub("3", "", .)

#Anova
aov_rr3_model <- as.data.frame(
  sapply(
    rr3_aov[2:ncol(rr3_aov)],
    function(x) summary(aov(x ~ rr3_aov[["Organism"]]))[[1]][1,'Pr(>F)']))

aov_rr3_df <- aov_rr3_model%>%
  rownames_to_column(var = "feature_number")%>%
  rename("p_value" = 2)



## Running Bejamini-Hochberg False Discovery Rate Corrections and filtering to significant values
aov_rr3_df$FDR <- p.adjust(aov_rr3_df$p_value, method = "BH")

one_way_signignificant <- aov_rr3_df%>%
  filter(FDR < 0.05)

#Save siginificant values from test to a vector to filter out for post-hoc tests
significant_aov_rr3 <- as.vector(one_way_signignificant$feature)


## Dunnetts
rr3_dunnetts <- lapply(rr3_aov[2:ncol(rr3_aov)], function(y, f)
  summary(glht(aov(y ~ f, rr3_aov), 
               linfct = mcp(f = "Dunnett"))),
  f = 
    as.factor(
      relevel(
        factor(rr3_aov$Organism), "WAT")))

pvals_dunn_rr3 <- as.data.frame(
  sapply(
    rr3_dunnetts,
    function(x) x$test$pvalues))

pvals_dunn_rr3 <- cbind("Organism" = c("CCA",
                                       "Dictyota", 
                                       "Pocillopora verrucosa",
                                       "Porites Lobata",
                                       "Turf"),
                        pvals_dunn_rr3)

rr3_dunnetts_FDR <- pvals_dunn_rr3%>%
  gather(feature_number, p_value, 2:ncol(.))

rr3_dunnetts_FDR$FDR_f <-p.adjust(rr3_dunnetts_FDR$p_value, method = "BH")

rr3_dunnett_sigs <- rr3_dunnetts_FDR%>%
  filter(FDR_f, FDR_f < 0.05)%>%
  dplyr::select(-p_value)%>%
  spread(Organism, FDR_f)

## Tukeys
tukey_rr3_model <- sapply(rr3_aov[2:ncol(rr3_aov)], function(x)
  TukeyHSD(aov(x ~ rr3_aov$Organism, data = rr3_aov), p.adjust.methods = "BH"))

p_values_tukey_rr3 <- as.data.frame(tukey_rr3_model)%>%
  rownames_to_column(var = "variable")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% '%p.adj%')%>%
  filter(value, value < 0.05)

p_values_tukey_rr3$feature_info <- p_values_tukey_rr3$feature_info%>%
  gsub(".rr3_aov.Organism.p.adj","", .)

# Tukeys for diff columns
tukey_rr3_ra_model <- sapply(rr3_ra[2:ncol(rr3_ra)], function(x)
  TukeyHSD(aov(x ~ rr3_ra$Organism, data = rr3_ra), p.adjust.methods = "BH"))

diff_tukey_rr3 <- as.data.frame(tukey_rr3_ra_model)%>%
  rownames_to_column(var = "variable")%>%
  gather(feature_info, value, 2:ncol(.))%>%
  filter(feature_info %like% '%diff%')

diff_tukey_rr3$feature_info <- diff_tukey_rr3$feature_info%>%
  gsub(".rr3_ra.Organism.diff","", .)

p_values_diffs_tukey_rr3 <- left_join(p_values_tukey_rr3, diff_tukey_rr3, by = c("variable", "feature_info" ))%>%
  rename(p_value = 3, 
         diff = 4)

## writing the sig tables
write_csv(p_values_diffs_tukey_rr3, "rr3_mapped_tukey_significant.csv")
write_csv(rr3_dunnett_sigs, "rr3_mapped_dunnetts_significant.csv")




