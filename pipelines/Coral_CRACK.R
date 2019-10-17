
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
  mutate(CCA_species = case_when(CCA_species == "HB" ~ "Hydrolithon borgesenii",
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
cca_species_vector <- as.vector(c("Hydrolithon borgesenii", "Paragoniolithon solubile", "Porolithon pachydernum"))

dunnett_model_scoring <- tidy_scoring%>%
  group_by(Coral_species, metabolite_pool)%>%
  nest()%>%
  mutate(dunnett = map(data, ~ aov(percent.asin ~ CCA_species, .x)%>%
                    glht(linfct = mcp(CCA_species = "Dunnett"))),
      dunnett_summary = map(dunnett, ~ summary(.x)%>%
                              tidy()))%>%
  dplyr::select(-c(data, dunnett))%>%
  unnest(dunnett_summary)%>%
  dplyr::select(-c(4:7))%>%
  mutate(lhs = gsub(" - Water control", "", lhs))

# GRAPHING -- Scoring Bar Charts ------------------------------------------
no_pachy <- tidy_scoring%>%
  filter(!CCA_species == "Porolithon pachydernum")

pdf("perent_settlement.pdf", height = 5, width = 8)
ggplot(no_pachy, aes(x = CCA_species, y = percent))+
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
    axis.text.x = element_text(angle = 75, hjust = 1,face = "italic"),
    axis.title = element_text(face = "italic"),
    strip.text = element_text(face = "italic")
  ) +
  facet_wrap(~ reorder(Coral_species, percent), nrow = 1) +
  xlab("CCA species") +
  ylab("Percent Settled")
dev.off()


# WRITING -- stats_dataframes ---------------------------------------------
write_csv(dunnett_model_scoring, "settlement_dunnet_pvalues.csv")
