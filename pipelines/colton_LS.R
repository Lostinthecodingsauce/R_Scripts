# Loading In Libraries and Datasheets -------------------------------------

library(purrr)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr) 
library(janitor)
library(readr)
library(magrittr)
library(DescTools)
library(reshape2)
library(ggpubr)

##color pallette
library(wesanderson)

##PCoA
library(ggplot2)
library(MASS)
library(ape)
library(vegan)
library(factoextra)
library(schoolmath)

colton_raw <- read_csv("All Leopard Shark Data.csv")
colton_meta <- read_csv("LS_Metadata.csv")


# Data Cleaning -----------------------------------------------------------


##    You had some duplicate Genera
ls_tran <-colton_raw %>%
  gather(Sample, RA, 2:44)%>%
  group_by(Genus, Sample)%>%
  summarize_if(is.numeric, sum)%>%
  ungroup(coltron_raw)%>%
  spread(Genus, RA)

## Combining with Meta Data
ls_df <- left_join(colton_meta, ls_tran, by = "Sample")


# PCoA --------------------------------------------------------------------

#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- as.data.frame(ls_df %>%
                        dplyr::select(-c(Sample:Year)),
                      dinames = list(paste("", 1:43, sep = ",")))

veg_bray <- vegdist(microb_f, "bray") #Bray-curtis distances

pc_scores <- pcoa(veg_bray) #Calculating scores

#This plots Eigenvalues
#They will allow you to choose the best axes to show how your data varies
ylimit = c(0, 1.1*max(pc_scores$values$Relative_eig))

Eigan <- barplot(pc_scores$values$Relative_eig[1:10], ylim= ylimit)
# Add values to the bars
text(x = Eigan, y = pc_scores$values$Relative_eig[1:10], label = pc_scores$values$Relative_eig[1:10], pos = 4, cex = .7, col = "red")

# S3 method for pcoa
biplot(pc_scores, Y=NULL, col = ls_wdf$Location, plot.axes = c(1,2), dir.axis1=1,
       dir.axis2=1)

pco_scores <- as.data.frame(pc_scores$vectors)

pco_scores$Sample <- ls_wdf$Sample
pco_scores$Group <- ls_wdf$Group
pco_scores$Location <- ls_wdf$Location
pco_scores$Year <- ls_wdf$Year

ggplot(pco_scores, mapping = aes(Axis.1, Axis.2, col = ls_wdf$Location, shape = factor(ls_wdf$Year))) +
  geom_point(stat = "identity") +
  scale_color_manual(values=wes_palette(n=4, name="Darjeeling1")) +
  scale_shape_manual(values = c(0,1,2,3,5))

write_csv(pco_scores, "colton_ls_pcoa.dat")


# PERMANOVA ---------------------------------------------------------------
adonis(microb_f ~ Year, ls_wdf, perm=1000, method="bray", set.seed(100))


# ANOVA -------------------------------------------------------------------
ls_df[is.nan(ls_df)] <-  0
ls_df[is.na(ls_df)] <- 0

ls_wdf <-as.data.frame(ls_df)

twoway_anova_rows <- c("Loaction", "Group", "Location*Group")
aov_ls <- sapply(ls_wdf[5:297], function(x) summary(aov(x ~ ls_wdf[["Location"]]*ls_wdf[["Group"]]))[[1]][1:3,'Pr(>F)'])


anova_ls <- as.data.frame(aov_ls)
anova_ls$Anova_test <- cbind(twoway_anova_rows)

anova_ls_tidy <- anova_ls%>%
  gather(Clade, F_value, 1:293)
