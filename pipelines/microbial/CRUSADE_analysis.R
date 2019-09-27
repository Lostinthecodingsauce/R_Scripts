## Originally written by Zach Quinlan 12/17/2018
## Editted and re-written by Zach Quinlan (1/30_2019) for CRUSADE manuscript
## This should be able to take the raw outputs from mothur and give a working df

# LOADING -- libraries ----------------------------------------------------
# Libraries for cleaning
library(DescTools)
library(tidyverse)
library(data.table)

# Libraries for NMDS, PERMANOVA, PCoA
library(vegan)
library(MASS)
library(factoextra)
library(ape)

# LOADING -- Read in data frames  -----------------------------------------------------
# These are the raw outputs from mothur
taxon_file <- read_tsv("Undetermined_S0_L001_R1_001.trim.contigs.good.good.unique.filter.good.precluster.pick.pick.subsample.opti_mcc.0.03.cons.taxonomy.txt")
ra_file <- read_tsv("Undetermined_S0_L001_R1_001.trim.contigs.good.good.unique.filter.good.precluster.pick.pick.subsample.opti_mcc.relabund.txt")%>%
  gather(OTU, RA, 4:6522)%>%
  spread(Group, RA)

# DATAFRAME PREP -- Join two tables and seperate taxonomy -----------------------------------
crusade_df <-left_join(taxon_file, ra_file, by = "OTU")%>%
  separate(Taxonomy, paste("Taxonomy", 1:7, sep = "."), sep = ";")%>%
  dplyr::select(-c(OTU, Size, label, numOtus))%>%
  dplyr::rename("Kingdom" = "Taxonomy.1",
                "Phylum"  = "Taxonomy.2",
                "Class"   = "Taxonomy.3",
                "Order"   = "Taxonomy.4",
                "Family"  = "Taxonomy.5",
                "Genus"   = "Taxonomy.6",
                "OTU"     = "Taxonomy.7")

# DATAFRAME PREP -- Cleaning all naming columns ---------------------------------------------
## These are changing all unidentified classifiers to a single value so that 
## If there are two different unclassified/unidentified OTU of the same genus they will be combined into one unidentified genus
cleaned <- crusade_df%>%
  dplyr::mutate(Class = case_when(Class %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Class)))%>%
  dplyr::mutate(Order = case_when(Class == "unclassified" ~ "",
                                  TRUE ~ as.character(Order)))%>%
  dplyr::mutate(Family = case_when(Class == "unclassified" ~ "",
                                   TRUE ~ as.character(Family)))%>%
  dplyr::mutate(Genus = case_when(Class == "unclassified" ~ "",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Class == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(Order = case_when(Order %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Order)))%>%
  dplyr::mutate(Family = case_when(Order == "unclassified" ~ "",
                                   TRUE ~ as.character(Family)))%>%
  dplyr::mutate(Genus = case_when(Order == "unclassified" ~ "",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Order == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(Family = case_when(Family %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                   TRUE ~ as.character(Family)))%>%
  dplyr::mutate(Genus = case_when(Family == "unclassified" ~ "",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Family == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(Genus = case_when(Genus %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "unclassified",
                                  TRUE ~ as.character(Genus)))%>%
  dplyr::mutate(OTU = case_when(Genus == "unclassified" ~ "",
                                TRUE ~ as.character(OTU)))%>%
  dplyr::mutate(OTU = case_when(OTU %like any% c("%uncultured%", "%unclassified%", "%unidentified%") ~ "sp",
                                TRUE ~ as.character(OTU)))
## Making spread sheets from cleaned data. 
## Class_order for pie charts, OFGO for graphing/relative abun data, transformed for stats
class_order <-cleaned%>%
  group_by(Class, Order) %>%
  summarize_if(is.numeric, sum)%>%
  ungroup(cleaned)

OFGO <- cleaned%>%
  group_by(Order, Family, Genus, OTU) %>%
  summarize_if(is.numeric, sum)%>%
  ungroup(cleaned)%>%
  unite(OFGO, c(Order, Family, Genus, OTU), sep = ";", remove = TRUE)

transformed <- OFGO%>%  
  mutate_if(is.numeric, sqrt)%>% #Normalizing data (arcsin(sqrt(x)))
  mutate_if(is.numeric, asin)

## Transposing dataframes
transposed_transformed <-transformed%>%
  gather(sample_code, RA, 2:32)%>%
  spread(OFGO, RA)

transposed_graphing <- OFGO%>%
  gather(sample_code, RA, 2:32)%>%
  spread(OFGO, RA)

# DATAFRAME PREP -- Designing working df for stats and such ---------------------------------
working_df <- transposed_transformed%>%
  separate(sample_code, paste("C", 1:4, sep = "."), sep = "_", remove = FALSE)%>%
  rename(Organism = `C.1`,
         Water = `C.2`,
         Timepoint = `C.3`,
         DNA.Source = `C.4`)%>%
  mutate(Organism = case_when(Organism == "CC" ~ "Calcium Carbonate Control",
                              Organism == "CH" ~ "Water Control",
                              Organism == "HR" ~ "Hydrolithon reinboldii",
                              Organism == "PO" ~ "Porolithon onkodes",
                              TRUE ~ "Start Water"))%>%
  mutate(Water = case_when(Water == "1F" ~ "Filtered",
                           Water == "2F" ~ "Filtered",
                           Water == "3F" ~ "Filtered",
                           TRUE ~ "Unfiltered"))%>%
  mutate(DNA.Source = case_when(DNA.Source == "CCA" ~ "CCA Microbiome",
                                TRUE ~ "Bacterioplankton"))%>%
  unite(Name, c(Organism, Water, Timepoint, DNA.Source), sep = "      ", remove = FALSE)

graphing_df <- transposed_graphing%>%
  separate(sample_code, paste("C", 1:4, sep = "."), sep = "_", remove = FALSE)%>%
  rename(Organism = `C.1`,
         Water = `C.2`,
         Timepoint = `C.3`,
         DNA.Source = `C.4`)%>%
  mutate(Organism = case_when(Organism == "CC" ~ "Calcium Carbonate Control",
                              Organism == "CH" ~ "Water Control",
                              Organism == "HR" ~ "Hydrolithon reinboldii",
                              Organism == "PO" ~ "Porolithon onkodes",
                              TRUE ~ "Start Water"))%>%
  mutate(Water = case_when(Water == "1F" ~ "Filtered",
                           Water == "2F" ~ "Filtered",
                           Water == "3F" ~ "Filtered",
                           TRUE ~ "Unfiltered"))%>%
  mutate(DNA.Source = case_when(DNA.Source == "CCA" ~ "CCA Microbiome",
                                TRUE ~ "Bacterioplankton"))%>%
  unite(Name, c(Organism, Water, Timepoint, DNA.Source), sep = "      ", remove = FALSE)

write_csv(working_df, "Crusade_wdf.dat ")
write_csv(graphing_df, "Crusade_graphing_df.dat")
# DATAFRAME PREP -- Subsetted dataframes-------------------------------------------------------
# Bacterioplankton whole df
bacterioplankton_df <- working_df %>%
  dplyr::filter(Timepoint == "T3") %>%
  dplyr::filter(DNA.Source == "Bacterioplankton") %>%
  dplyr::filter(!Organism == "Start Water") %>%
  dplyr::select(-c(Name, Timepoint, DNA.Source))


# Bacterioplanknton unfiltered
bacterio_unfilt <- bacterioplankton_df %>%
  dplyr::filter(Water == "Unfiltered")

bacterio_un_abun <- bacterio_unfilt %>%
  dplyr::select(-c(Water, Organism, sample_code))


# Bacterioplankton Filtered
bacterio_filt <- bacterioplankton_df %>%
  dplyr::filter(!Water == "Unfiltered")

bacterio_f_abun <- bacterio_filt %>%
  dplyr::select(-c(Water,Organism, sample_code))

# CCA Microbiome
cca_microb <- working_df %>%
  dplyr::filter(Timepoint == "T3") %>%
  dplyr::filter(!DNA.Source == "Bacterioplankton") %>%
  dplyr::filter(!Organism == "Start Water") %>%
  dplyr::select(-c(Name, Timepoint, DNA.Source, Water))

cca_microb_abun <- cca_microb %>%
  dplyr::select(-c(Organism, sample_code))

# STATS -- PERMANOVAs ------------------------------------------------------
# Bacterioplankton unfiltered
adonis(bacterio_un_abun ~ Organism, bacterio_unfilt, perm=1000, method="bray", set.seed(100))

# Bacterioplankton filtered
adonis(bacterio_f_abun ~ Organism, bacterio_filt, perm=1000, method="bray", set.seed(100))

# CCA Microbiome
adonis(cca_microb_abun ~ Organism, cca_microb, perm=1000, method="bray", set.seed(100))

# STATS -- Two-Way ANOVAs OTU -----------------------------------------------------
only_cca <- working_df%>%
  dplyr::filter(!Organism == "Calcium Carbonate Control",
                !Organism == "Water Control",
                !Organism == "Start Water")

##  subset the data by bacterioplankton or microbiome
occa_mb <- only_cca%>%
  dplyr::filter(!DNA.Source == "Bacterioplankton")

occa_bp <- only_cca%>%
  dplyr::filter(DNA.Source == "Bacterioplankton")

##Running the two way anova for Bacterioplankton Organism*Water

# This line makes the names of the rows which will be added into the pvalue table
twoway_anova_rows <- c("Organism", "Water", "Organism*Water")

# This is the line which applies the anova over every single clade in the bacterioplankton and runs the two way anova
aov_bp <- sapply(occa_bp[8:667], function(x) summary(aov(x ~ occa_bp[["Organism"]]*occa_bp[["Water"]]))[[1]][1:3,'Pr(>F)'])

aov_mb <-sapply(occa_mb[8:667], function(x) summary(aov(x ~ occa_mb[["Organism"]]))[[1]][1,'Pr(>F)'])
# It comes out as a list so it has to be first converted to a data frame before we can add in the test names as a column
anova_bp <- as.data.frame(aov_bp)
anova_bp$Anova_test <- cbind(twoway_anova_rows)

anova_mb <- as.data.frame(aov_mb)
# Now the data is made Tidy and we filter to only significant values
anova_bp_tidy <- anova_bp%>%
  gather(Clade, F_value, 1:660)

anova_bp_tidy$Water <- "Bacterioplankton"

anova_mb_tidy <- anova_mb%>%
  rownames_to_column("Clade")%>%
  rename(F_value = 'aov_mb')

anova_mb_tidy$Water <- "Microbiome"

# Subset the data so we are only looking at interaction of Organism*Water
anbp_Or_Wa <- anova_bp_tidy%>%
  filter(Anova_test == "Organism*Water")

## Combining Pvalues, FDR correction and finding sigs
anova_all <- bind_rows(anova_bp_tidy, anova_mb_tidy)
anova_all$FDR_f <- p.adjust(anova_all$F_value, method = "BH")

anova_FDR <- anova_all%>%
  filter(FDR_f, FDR_f < 0.05)

## Printing Anova Sig tables
write_csv(anova_FDR, "ANOVA_FDRsigs_All.dat")

# META-STATS -- Mean RA’s of Significant OTUs -------------------------------------------
sig_otus <- as.vector(anova_FDR$Clade)

sig_otus_microbiome <- anova_FDR%>%
  filter(Water == "Microbiome")

sig_microbiome_clades <- as.vector(sig_otus_microbiome$Clade)

sig_otu_table <- graphing_df%>%
  dplyr::select(1:6,sig_otus)%>%
  dplyr::filter(!Organism == "Water Control",
                !Organism == "Calcium Carbonate Control",
                !Organism == "Start Water")

sig_otu_table_microbiome <- graphing_df%>%
  dplyr::select(1:6,sig_microbiome_clades)%>%
  dplyr::filter(!Organism == "Water Control",
                !Organism == "Calcium Carbonate Control",
                !Organism == "Start Water")



water_sig_otu <- sig_otu_table%>%
  filter(DNA.Source == "Bacterioplankton")%>%
  group_by(Water)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(Clade, RA, 2:28)%>%
  spread(Water, RA)
water_sig_otu$DNA <- "Bacterioplankton"

organism_sig_otu <- sig_otu_table%>%
  filter(DNA.Source == "Bacterioplankton")%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(Clade, RA, 2:28)%>%
  spread(Organism, RA)
organism_sig_otu$DNA <- "Bacterioplankton"

microbiome_sig_otu <- sig_otu_table_microbiome%>%
  filter(!DNA.Source == "Bacterioplankton")%>%
  group_by(Organism)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  gather(Clade, RA, 2:5)%>%
  spread(Organism, RA)
microbiome_sig_otu$DNA <- "microbiome"

mean_sigs_bacterioplaknton <- full_join(water_sig_otu, organism_sig_otu)

write_csv(mean_sigs_bacterioplaknton, "significant_OTU_means_bact.dat")
write_csv(microbiome_sig_otu, "significant_OTU_means_micro.dat")

# VISUALIZATIONS -- Making PCoA's-----------------------------------------------------------
bacterioplankton_noh20 <-bacterioplankton_df%>%
  filter(Organism != "Water Control")
#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- cca_microb%>%
  dplyr::select(-c(Organism, sample_code))
  

veg_bray <- vegdist(microb_f, "bray") #Bray-curtis distances

pc_scores<-pcoa(veg_bray) #Calculating scores

#This plots Eigenvalues
#They will allow you to choose the best axes to show how your data varies
ylimit = c(0, 1.1*max(pc_scores$values$Relative_eig))

Eigan <- barplot(pc_scores$values$Relative_eig[1:10], ylim= ylimit)
# Add values to the bars
text(x = Eigan, y = pc_scores$values$Relative_eig[1:10], label = pc_scores$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

# S3 method for pcoa
biplot(pc_scores, Y=NULL, col = microb_f_abun$Organism, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1)

pco_scores <- as.data.frame(pc_scores$vectors)
pco_scores$Sample.Code <- cca_microb$sample_code     # This will add reference labels to the PCoA scores
pco_scores$Organism <- cca_microb$Organism
# pco_scores$Water <- bacterioplankton_noh20$Water

ggplot(pco_scores, aes(x = Axis.1, y = Axis.2, color = Organism, shape = Water)) +
  geom_point(stat = "identity", aes(size = 10)) +
  scale_color_manual(values = c("darkgoldenrod3", 'deepskyblue', "lightcoral"))+ 
  scale_shape_manual(values = c(16,1))+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.line = element_line(color="black")
  ) +
  xlab("PCoA 1 (28.6%)") +
  ylab("PCoA 2 (17.3%)")

tiff("bacterioplankton.tiff", units="in", width=8, height=5, res=300)
ggplot(pco_scores, aes(x = Axis.1, y = Axis.2, color = Organism, shape = Water)) +
  geom_point(stat = "identity", aes(size = 10)) +
  scale_color_manual(values = c("darkgoldenrod3", 'deepskyblue', "lightcoral"))+ 
  scale_shape_manual(values = c(16,1))+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.x = element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.line = element_line(color="black")
  ) +
  xlab("PCoA 1 (28.6%)") +
  ylab("PCoA 2 (17.3%)")
dev.off()


write.csv(pco_scores, "PCo_scores.dat")           # This will print the Scores to remake the figure in JMP because JMP is bae
# DATAFRAME PREP -- Making Pie Charts data frames -------------------------------------------------------
class_renamed <-class_order%>%
  mutate(Class = case_when(Class == "Alphaproteobacteria" ~ "Alphaproteobacteria",
                           Class == "Cyanobacteria" ~ "Cyanobacteria",
                           Class == "Cytophagia" ~ "Cytophagia",
                           Class == "Deltaproteobacteria" ~ "Deltaproteobacteria",
                           Class == "Flavobacteria" ~ "Flavobacteria",
                           Class == "Gammaproteobacteria" ~ "Gammaproteobacteria",
                           Class == "Planctomycetacia" ~ "Planctomycetacia",
                           Class == "Sphingobacteriia" ~ "Sphingobacteriia",
                           TRUE ~ "Other"))

raw_pie <- class_renamed%>%
  gather(Code, RA, 3:33)%>%
  separate(Code, paste("C", 1:4, sep = "."), sep = "_")%>%   
  rename(Organism = `C.1`,
         Water = `C.2`,
         Timepoint = `C.3`,
         DNA.Source = `C.4`)%>%
  mutate(Organism = case_when(Organism == "CC" ~ "Calcium Carbonate Control",
                              Organism == "CH" ~ "Water Control",
                              Organism == "HR" ~ "Hydrolithon reinboldii",
                              Organism == "PO" ~ "Porolithon onkodes",
                              TRUE ~ "Start Water"))%>%
  mutate(Water = case_when(Water == "1F" ~ "R1_Filtered",
                           Water == "2F" ~ "R2_Filtered",
                           Water == "3F" ~ "R3_Filtered",
                           Water == "1U" ~ "R1_Unfiltered",
                           Water == "2U" ~ "R2_Unfiltered",
                           TRUE ~ "R3_Unfiltered"))%>%
  separate(Water, paste(c("Replicate", "Water")), sep = "_")%>%
  mutate(DNA.Source = case_when(DNA.Source == "CCA" ~ "CCA Microbiome",
                                TRUE ~ "Bacterioplankton"))

#For this next part I want to make a pie chart for each DNA source and for each organism
#Below is the exanple of how I began to get the dataframe for just bacterioplankton in filtered water treatments
#Filter by DNA source and group by organism
#Only from filtered bacterioplankton. No Water control or Start Water
cr_bact_filt_combine <- raw_pie%>%
  dplyr::filter(DNA.Source == "Bacterioplankton")%>%
  dplyr::filter(Water == "Filtered")%>%
  dplyr::filter(!Organism == "Start Water",
                !Organism == "Water Control")%>%
  group_by(Organism, Order, Class)%>%
  summarize_if(is.numeric, mean)%>% ##This will take the mean of CCA species orders 
  ungroup(cr_bact_filt_combine)

cr_CCA_combine <- raw_pie%>%
  dplyr::filter(!DNA.Source == "Bacterioplankton")%>%
  dplyr::filter(Water == "Filtered")%>%
  dplyr::filter(!Organism == "Start Water",
                !Organism == "Water Control",
                !Organism == "Calcium Carbonate Control")%>%
  group_by(Organism, Order, Class)%>%
  summarize_if(is.numeric, mean)%>% ##This will take the mean of CCA species orders 
  ungroup(raw_pie)

## Making individual data frames for each species treatment
PO_filt <- cr_bact_filt_combine%>%
  dplyr::filter(Organism == "Porolithon onkodes")

HR_filt <- cr_bact_filt_combine%>%
  dplyr::filter(Organism == "Hydrolithon reinboldii")

HR_CCA <-cr_CCA_combine%>%
  dplyr::filter(Organism == "Hydrolithon reinboldii")

PO_CCA <-cr_CCA_combine%>%
  dplyr::filter(Organism == "Porolithon onkodes")

# DATAFRAME PREP -- FINDING ONLY Classes WITH GREATER THAN 2% ABUNDNACE  -------------------
#Filtering to only X Rabundances, tried to limit to ~10 different classes
# This can be used to find the 
# bact_filt_02 <- cr_bact_filt_combine%>%
#   gather(Class, RA, 2:65)%>%
#   group_by(Organism, Class)%>%
#   dplyr::filter(RA, sum(RA[RA>0.02]))%>%
#   ungroup(bact_filt_02)

# CCA_02_a <- cr_CCA_combine%>%
#   group_by(Organism, Order)%>%
#   dplyr::filter(RA, sum(RA[RA>0.02]))%>%
#   ungroup(cr_CCA_combine)%>%
#   dplyr::filter(!Organism == "Calcium Carbonate Control")

# VIZUALIZATION -- PIE CHARTS ----------------------------------
#Visualize Pie Charts
#Filtered treatments
PO_filt_ggp <-ggplot(PO_filt, aes(x= "", y=RA, fill = Class, colour = pal)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "Porolithon onkodes filtered") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

HR_filt_ggp <-ggplot(HR_filt, aes(x= "", y=RA, fill = Class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "Hydrolithon reinboldii filtered") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

PO_CCA_ggp <-ggplot(PO_CCA, aes(x= "", y=RA, fill = Class)) +
  geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
  labs( x =NULL, y = NULL, title = "Hydrolithon reinboldii CCA") + 
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

PO_filt_pie <- PO_filt_ggp + coord_polar("y", start = 0)
HR_filt_pie <- HR_filt_ggp + coord_polar("y", start = 0)
HR_CCA_pie <- PO_CCA_ggp + coord_polar("y", start = 0)

## To visualize, uncomment the next line
PO_filt_pie
# HR_filt_pie

#Saving Pie Charts as pdfs to wd
# ggsave(filename = "PO_filt.pdf", plot= PO_filt_pie)
# ggsave(filename = "HR_filt.pdf", plot= HR_filt_pie)


# STATS -- ANOVA Order -------------------------------------------------------------
Order_filtered <-class_order %>%
  group_by(Order)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup(class_order)%>%
  mutate_if(is.numeric, sqrt)%>%
  mutate_if(is.numeric, asin)%>%
  gather(Sample_code, RA, 2:32)%>%
  spread(Order, RA)

Order_wdf <- Order_filtered%>%
  separate(Sample_code, paste("C", 1:4, sep = "."), sep = "_")%>%   
  rename(Organism = `C.1`,
         Water = `C.2`,
         Timepoint = `C.3`,
         DNA.Source = `C.4`)%>%
  mutate(Organism = case_when(Organism == "CC" ~ "Calcium Carbonate Control",
                              Organism == "CH" ~ "Water Control",
                              Organism == "HR" ~ "Hydrolithon reinboldii",
                              Organism == "PO" ~ "Porolithon onkodes",
                              TRUE ~ "Start Water"))%>%
  mutate(Water = case_when(Water == "1F" ~ "Filtered",
                           Water == "2F" ~ "Filtered",
                           Water == "3F" ~ "Filtered",
                           TRUE ~ "Unfiltered"))%>%
  mutate(DNA.Source = case_when(DNA.Source == "CCA" ~ "CCA Microbiome",
                                TRUE ~ "Bacterioplankton"))%>%
  dplyr::filter(!Organism == "Start Water",
                !Organism == "Water Control",
                !Organism == "Calcium Carbonate Control") ##NOTE WHAT IS FILTERED OUT

Order_filtered <-Order_wdf%>%
  dplyr::filter(DNA.Source == "Bacterioplankton")%>%
  dplyr::filter(Water == "Filtered")

Order_unfiltered <- Order_wdf%>%
  dplyr::filter(DNA.Source == "Bacterioplankton")%>%
  dplyr::filter(Water == "Unfiltered")

Order_CCA <- Order_wdf%>%
  dplyr::filter(DNA.Source == "CCA Microbiome")


# STATS -- ANOVAS ORDER summed data -----------------------------
##This appears breifly in the manuscript when I cite the one-way anovas for Orders in relation to the pie charts
orders <- c("Rhodospirillales", "Rhodobacterales", "Rhizobiales", "Alteromonadales",
            "Flavobacteriales", "Sphingobacteriales", 
            "Cytophagales", "Vibrionales", "Oceanospirillales", "SAR11_clade", 
            "SubsectionI", "Acidimicrobiales")

abun2_anova_filtered <- sapply(Order_filtered[orders], function(x) summary(aov(x ~ Order_filtered[["Organism"]]))[[1]][1,'Pr(>F)'])
ab2_an_f_matrix <- as.data.frame(abun2_anova_filtered)%>%
  rownames_to_column(var = "Clade")%>%
  rename(F_value = 'abun2_anova_filtered')
ab2_an_f_matrix$Water <- "Filtered"

abun2_anova_unfiltered <- sapply(Order_unfiltered[orders], function(x) summary(aov(x ~ Order_unfiltered[["Organism"]]))[[1]][1,'Pr(>F)'])
ab2_an_u_matrix <- as.data.frame(abun2_anova_unfiltered)%>%
  rownames_to_column(var = "Clade")%>%
  rename(F_value = 'abun2_anova_unfiltered')
ab2_an_u_matrix$Water <- "Unfiltered"

abun2_anova_cca<- sapply(Order_CCA[orders], function(x) summary(aov(x ~ Order_CCA[["Organism"]]))[[1]][1,'Pr(>F)'])
ab2_an_mi_matrix <- as.data.frame(abun2_anova_cca)%>%
  rownames_to_column(var = "Clade")%>%
  rename(F_value = 'abun2_anova_cca')
ab2_an_mi_matrix$Water <- "Microbiome"

Anova_all_orders <- bind_rows(ab2_an_f_matrix, ab2_an_u_matrix, ab2_an_mi_matrix)

Anova_all_orders$FDR <- p.adjust(Anova_all_orders$F_value, method = "BH")

anova_FDR_orders <- Anova_all_orders%>%
  filter(FDR, FDR < 0.05)

write_csv(Order_wdf, "CRUSADE_Order_WDF.dat")

# META-STATS -- Class Standard Deviations -----------------------------------------------

class_only <- raw_pie%>%
  dplyr::select(-c(Order, Timepoint))%>%
  group_by(Class, Replicate, Organism, DNA.Source, Water)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup(raw_pie)

class_only[is.na(class_only)] <- 0

means_class <- class_only%>%
  spread(Replicate, RA)%>%
  ungroup(class_only)%>%
  dplyr::filter(!Organism == "Water Control",
                !Organism == "Calcium Carbonate Control",
                !Organism == "Start Water")

means_bacter <- class_only%>%
  dplyr::filter(!Organism == "Water Control",
                !Organism == "Calcium Carbonate Control",
                !Organism == "Start Water")%>%
  dplyr::filter(DNA.Source == "Bacterioplankton")%>%
  unite(Name, c(Organism, Replicate), sep = "_")%>%
  spread(Name, RA)

means_bacter[is.na(means_bacter)] <- 0

means_CCA <- class_only%>%
  dplyr::filter(!Organism == "Water Control",
                !Organism == "Calcium Carbonate Control",
                !Organism == "Start Water")%>%
  dplyr::filter(!DNA.Source == "Bacterioplankton")%>%
  unite(Name, c(Organism, Replicate), sep = "_")%>%
  spread(Name, RA)

means_CCA[is.na(means_CCA)] <- 0

means_class$Mean <- apply(means_class[5:7], 1, FUN = mean)
means_class$std <- apply(means_class[5:7], 1, FUN = sd)

means_bacter$Mean <- apply(means_bacter[4:9], 1, FUN = mean)
means_bacter$std <- apply(means_bacter[4:9], 1, FUN = sd)

means_CCA$Mean <- apply(means_CCA[4:9], 1, FUN = mean)
means_CCA$std <- apply(means_CCA[4:9], 1, FUN = sd)


# SUPPLEMENTAL TABLE 1 ----------------------------------------------------
library(tidyverse)

fd_fc <- read_csv("CRUSADE_fDOM_FCM.csv")%>%
  dplyr::select(-c(7:47))

fd_fc_clean <- fd_fc%>%
  dplyr::select(-c(1,2,Replicate,14:ncol(.)))%>%
  filter(!Timepoint == "2",
         !Inhabitant == "Water Control")%>%
  gather(fluor, RA_cm2, 4:ncol(.))%>%
  unite(name, c(Inhabitant, Water, fluor, Timepoint), sep = "_")%>%
  group_by(name)%>%
  summarize_if(is.numeric, mean)%>%
  ungroup()%>%
  separate(name, c('Inhabitant', 'Water', 'fluor', 'Timepoint'), sep = "_")%>%
  mutate(Timepoint = case_when(Timepoint == "1" ~ "one",
                               Timepoint == "3" ~ "three",
                               TRUE ~ as.character(Timepoint)))%>%
  spread(Timepoint, RA_cm2)%>%
  add_column(rate = (.$three-.$one)/8)

supp_table <- fd_fc_clean%>%
  dplyr::select(-c("three","one"))%>%
  spread(fluor, rate)

ggplot(fd_fc_clean, aes(x = Inhabitant, y = rate, fill = Water))+
  geom_boxplot(stat = "boxplot", fun.y = "mean", position = "dodge2")+
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
  facet_wrap(~ fluor) +
  xlab("Organism") +
  ylab("fluorophore (R.U. cm-2 h-1")


write_csv(supp_table, "Supplemental_table1.csv")


