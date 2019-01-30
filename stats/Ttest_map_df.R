##TTest map script
# T_test using map_df function for all anova significant clades (MIGHT NOT HAVE TO RUN T-TESTS)
# T_test for interactions organism
# Filtered
occa_bp_filt <- occa_bp%>%
  dplyr::filter(Water == "Filtered")

occa_bp_unfi <- occa_bp%>%
  dplyr::filter(Water == "Unfiltered")

ttest_bp_org_filt <- occa_bp_filt%>%
  dplyr::select(asigs_bp_or_wa)%>%
  select_if(is.numeric)%>%
  map_df(~broom::tidy(t.test(. ~ occa_bp_filt$Organism)))

ttest_bp_org_filt$Clade <- anbp_Or_Wa$Clade
ttest_bp_org_filt$Water <- "Filtered"

#Unfiltered
ttest_bp_org_unfi <- occa_bp_unfi%>%
  dplyr::select(asigs_bp_or_wa)%>%
  select_if(is.numeric)%>%
  map_df(~broom::tidy(t.test(. ~ occa_bp_unfi$Organism)))

ttest_bp_org_unfi$Clade <- anbp_Or_Wa$Clade
ttest_bp_org_unfi$Water <- "Unfiltered"
