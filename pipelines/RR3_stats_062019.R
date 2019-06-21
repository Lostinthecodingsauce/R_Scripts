## RR3 Stats
#Read it in
rr3_features <- read_tsv("AnnotatedRR3FeaturesForZach.txt")

moorea_binned <- moorea_post_stats[c(1, 948:963)]
  


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
  rename(p_value_water = 2)


porites_day_v_night <- 
  as.data.frame(
    sapply(
      porites_dvn[4:ncol(porites_dvn)],
      function(x) summary(aov(x ~ porites_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value_porities = 2)%>%
  

pocil_day_v_night <- 
  as.data.frame(
    sapply(
      pocil_dvn[4:ncol(pocil_dvn)],
      function(x) summary(aov(x ~ pocil_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value_pocillopora = 2)

cca_day_v_night <- 
  as.data.frame(
    sapply(
      cca_dvn[4:ncol(cca_dvn)],
      function(x) summary(aov(x ~ cca_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value_cca = 2)

turf_day_v_night <- 
  as.data.frame(
    sapply(
      turf_dvn[4:ncol(turf_dvn)],
      function(x) summary(aov(x ~ turf_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value_turf = 2)

dictyota_day_v_night <- 
  as.data.frame(
    sapply(
      dictyota_dvn[4:ncol(dictyota_dvn)],
      function(x) summary(aov(x ~ dictyota_dvn[["DayNight"]]))[[1]][1,'Pr(>F)']))%>%
  rownames_to_column(var = "feature_number")%>%
  rename(p_value_dictyota = 2)


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


