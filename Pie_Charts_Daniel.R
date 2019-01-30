## Okay so this script should make it so you can go from an unclean Relative abundance 16s dataframe to pretty graphs
## You need the taxonomy seperated out into its seven rows so that you can then filter to only include samples which are defined to the order level.
## Then because a lot of different people helped to make the SILVA files there are a lot of Uncultured, unclassified, crap.


#Loading in our nice little libraries
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(magrittr)
library(readr)
library(DescTools)

## Pie Chart Libraries
library(janitor)
library(ggplot2)
library(wesanderson)

## PCoA Libraries
library(MASS)
library(ape)
library(vegan)
library(factoextra)

####              DATA CLEANING           ####

crusade_raw_df <- read.csv("Cr_OTU_tax_RA.dat")%>%  ## This dataframe has your sample codes as variables and Taxonomies as observations
  rename(kingdom = `Taxonomy.1`,
         phylum  = `Taxonomy.2`,
         Class   = `Taxonomy.3`,
         Order   = `Taxonomy.4`,
         family  = `Taxonomy.5`,
         genus   = `Taxonomy.6`,
         otu     = `Taxonomy.7`)
##    This dataframe has samples now as observations and taxonomies as variables
##    It is classified to the OTU level and cleaned simiarly as I do below with order except with OTU then tranposed
Cr_PCoA_Data<-read.csv("CRUSADE_RA_NT.dat")

# filter to order and sum Unclassifieds together
# This data will be used for PIE CHARTS
data_sum_order <- crusade_raw_df %>%
  dplyr::select(c(Class, Order), c(8:38))%>%
  dplyr::mutate(Order = case_when(Order %like any% c("%uncultured%", "%unclassified%", "%Subgroup%") ~ "unclassified",
                                       TRUE ~ as.character(Order))) %>%
  group_by(Class, Order) %>%
  summarize_if(is.numeric, sum)%>%
  ungroup(data_sum_order)

##  Filter data for only the groups I want
##  Transforming to asin(sqrt)
##  This data will be used for PCoA's and NMDS
microb_f_abun <-Cr_PCoA_Data %>%
  mutate_if(is.numeric, sqrt)%>% #Normalizing data (arcsin(sqrt(x)))
  mutate_if(is.numeric, asin) %>%
  dplyr::filter(DNA.Source == "Bacterioplankton") %>%
  dplyr::filter(!Organism == "Water Control", !Organism == "State Water", !Organism == "Start Water")


####              PIE CHART DATA CLEANING          ####


##   This next portion is simply for making pie charts. Used in CRUSADE to keep specific classes for coloring pie chart wedges
##    Daniel, I would reccomend you find all classes which are 2% abundance or greater in your data which you will then color.
##    Everything else will be classified as "other".
##    Below I have a section where you will be able to do this and then you have to come back up to this section
##    The old section is commented out because once you find it you can just go forward with the script from here.


##    Once you have found these you can then change them out for the classes I have below.
class_filter <- data_sum_order%>%
  mutate(Class = case_when(Class == "Alphaproteobacteria" ~ "Alphaproteobacteria",
                           Class == "Cyanobacteria" ~ "Cyanobacteria",
                           Class == "Cytophagia" ~ "Cytophagia",
                           Class == "Deltaproteobacteria" ~ "Deltaproteobacteria",
                           Class == "Flavobacteria" ~ "Flavobacteria",
                           Class == "Gammaproteobacteria" ~ "Gammaproteobacteria",
                           Class == "Planctomycetacia" ~ "Planctomycetacia",
                           Class == "Sphingobacteriia" ~ "Sphingobacteriia",
                                TRUE ~ "Other"))

##  So here I am serparating out my sample code column and defining all the different codes
##  I am also making it so that there is only one column of relative abundances 
cr_raw_pie <- class_filter%>%
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
  mutate(Water = case_when(Water == "1F" ~ "Filtered",
                           Water == "2F" ~ "Filtered",
                           Water == "3F" ~ "Filtered",
                           TRUE ~ "Unfiltered"))%>%
  mutate(DNA.Source = case_when(DNA.Source == "CCA" ~ "CCA Microbiome",
                                TRUE ~ "Bacterioplankton"))

#For this next part I want to make a pie chart for each DNA source and for each organism
#Below is the exanple of how I began to get the dataframe for just bacterioplankton in filtered water treatments
#Filter by DNA source and group by organism
#Only from filtered bacterioplankton. No Water control or Start Water
cr_bact_filt_combine <-cr_raw_pie%>%
  dplyr::filter(DNA.Source == "Bacterioplankton")%>%
  dplyr::filter(Water == "Filtered")%>%
  dplyr::filter(!Organism == "Start Water",
                !Organism == "Water Control")%>%
  group_by(Organism, Order, Class)%>%
  summarize_if(is.numeric, mean)%>% ##This will take the mean of CCA species orders 
  ungroup(cr_bact_filt_combine)

## Making individual data frames for each species treatment
PO_filt <- cr_bact_filt_combine%>%
  dplyr::filter(Organism == "Porolithon onkodes")

HR_filt <- cr_bact_filt_combine%>%
  dplyr::filter(Organism == "Hydrolithon reinboldii")

####      FINDING ONLY Classes WITH GREATER THAN 2% ABUNDNACE      ####
#Filtering to only X Rabundances, tried to limit to ~10 different classes
# This can be used to find the 
# bact_filt_02 <- cr_bact_filt_combine%>%
#   gather(Class, RA, 2:65)%>%
#   group_by(Organism, Class)%>%
#   dplyr::filter(RA, sum(RA[RA>0.02]))%>%
#   ungroup(bact_filt_02)

####      ACTUAL VISUALIZATION OF THE PIE CHARTS    ####
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

PO_filt_pie <- PO_filt_ggp + coord_polar("y", start = 0)
HR_filt_pie <- HR_filt_ggp + coord_polar("y", start = 0)

## To visualize, uncomment the next line
##  PO_filt_pie

#Saving Pie Charts as pdfs to wd
# ggsave(filename = "PO_filt.pdf", plot= PO_filt_pie)
# ggsave(filename = "HR_filt.pdf", plot= HR_filt_pie)


####      PCoA MAKING   ####

#remove original dataframe
rm(Cr_PCoA_Data)

#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- as.matrix(microb_f_abun %>%
                        dplyr::select(-c(Label:DNA.Source)),
                      dinames = list(paste("Label", 1:31, sep = ",")))

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
pco_scores$Sample.Code <- microb_f_abun$Label     # This will add reference labels to the PCoA scores

write.csv(pco_scores, "PCo_scores.dat")           # This will print the Scores to remake the figure in JMP because JMP is bae


####          RUNNING A PERMANOVA TO CHECK FOR DIFFERENCES BETWEEN TREATMENTS       ####
#clean data?
#Need to have data completely cleaned of character data. Only have numerics.

crusade_abun <- microb_f_abun %>%
  dplyr::select(-c(Water,Organism, Label, DNA.Source, Timepoint))

##    This tests whether there is a significant difference between species using Bray-Curtis distances and 1000 permutations
adonis(crusade_abun ~ Organism, microb_f_abun, perm=1000, method="bray", set.seed(100))


