
# LOADING -- Libraries ----------------------------------------------------
library(tidyverse)
library(ggplot2)
library(wesanderson)
library(multcomp)

# LOADING -- Dataframes ---------------------------------------------------
raw_scoring_df <- read_csv("CRACK_Scoring.csv")%>%
  filter(!Sample_code == "PP_TI2_AP",
         !Sample_code == "PS_TI4_MF",
         !Sample_code == "PP_TI2_MF",
         !Sample_code == "HB_TI4_MF")

raw_scoring_df[is.na(raw_scoring_df)] <- 0


# CLEANING -- Scoring dataframe to only settled vs not settled-------------------------------------------
tidy_scoring <- raw_scoring_df%>%
  dplyr::select(-c(2:7))%>%
  add_column(percent = (.$Sum_settled/.$Alive)*100)%>%
  add_column(percent.asin = asin(sqrt((.$Sum_settled/.$Alive))))%>%
  separate(Sample_code, c("CCA_species", "metabolite_pool", "Coral_species"), sep = "_")%>%
  separate(metabolite_pool, c('metabolite_pool', 'replicate'), sep = 2)%>%
  mutate(CCA_species = case_when(CCA_species == "HB" ~ "Hydrolithon borgenesii",
                               CCA_species == "PS" ~ "Paragoniolithon solubile",
                               CCA_species == "PP" ~ "Porolithon pachydernum",
                               CCA_species == "WA" ~ "Water control",
                               TRUE ~ as.character(CCA_species)))%>%
  mutate(metabolite_pool = case_when(metabolite_pool == "EX" ~ "Exometabolites",
                                     metabolite_pool == "TI" ~ "Tissue metabolites",
                                     TRUE ~ as.character(metabolite_pool)))%>%
  mutate(Coral_species = case_when(Coral_species == "AP" ~ "Acropora palmata",
                                   Coral_species == "DL" ~ "Diploria labyrinthiformis",
                                   Coral_species == "MF" ~ "Montastrea faveolata",
                                   TRUE ~ as.character(Coral_species)))

tidy_scoring$CCA_species <- as.factor(relevel(factor(tidy_scoring$CCA_species), "Water control"))

# PRE-STATS -- Set seed ---------------------------------------------------
set.seed(29580)

# STATS -- TWO-Way ANOVA Scoring --------------------------------------------------
twoway_scoring_aov <- summary(aov(percent.asin ~ CCA_species*Coral_species, data = tidy_scoring))[[1]]['Pr(>F)']


# STATS -- Dunnetts scoring -----------------------------------------------
cca_species_vector <- as.vector(c("Hydrolithon borgenesii", "Paragoniolithon solubile", "Porolithon pachydernum"))

dunnett_model_scoring <- tidy_scoring%>%
  split(list(.$Coral_species, .$metabolite_pool), sep = "_")%>%
  map(~summary(glht(aov(percent.asin ~ CCA_species, .),
                    linfct = mcp(CCA_species = "Dunnett"))))

coral_labels <- as.vector(names(dunnett_model_scoring))

dunnett_pvalues <- as.data.frame(dunnett_model_scoring%>%
                                  map(~list(.$test$pvalues)),
                                row.names = cca_species_vector)

colnames(dunnett_pvalues) <- coral_labels

scoring_dunnett_pvalues <- dunnett_pvalues%>%
  rownames_to_column(var = "cca_species")%>%
  gather(spec_metab, p_value, 2:7)%>% 
  separate(spec_metab, c("coral_species", "metabolites_pool"), sep = "_")

# GRAPHING -- Scoring Bar Charts ------------------------------------------
ggplot(tidy_scoring, aes(x = CCA_species, y = percent))+
  geom_boxplot(aes(fill= metabolite_pool),
               stat = "boxplot", position = "dodge2") +
  geom_point(position=position_dodge(width=0.75), aes(group=metabolite_pool))+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major.y = element_line(size = 0.2, linetype = 'solid',colour = "gray"), # get rid of major grid
    # panel.grid.minor.x = element_line(size = 0.5, linetype = 'solid',colour = "black"), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
    axis.text.x = element_text(angle = 75, hjust = 1)
  ) +
  facet_wrap(~ Coral_species) +
  xlab("CCA species") +
  ylab("Percent Settled")


ggsave("percentsettlement.jpg", 
       plot = last_plot(), # or give ggplot object name as in myPlot,
       width = 6.5, height = 5, 
       units = "in", # other options c("in", "cm", "mm"), 
       dpi = 300)


# WRITING -- stats_dataframes ---------------------------------------------
write_csv(scoring_dunnett_pvalues, "settlement_dunnet_pvalues.csv")
