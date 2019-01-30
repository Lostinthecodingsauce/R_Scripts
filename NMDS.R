# import libraries

library(dplyr)
library(ggplot2)
library(readr)
library(MASS)
library(vegan)
library(factoextra)

#Read in data
Cr_Data<-read.csv("CRUSADE_RA_NT.dat")

#filter data for only the groups I want
microb_f_abun <-Cr_Data %>%
  mutate_if(is.numeric, sqrt)%>% #Normalizing data (arcsin(sqrt(x)))
  mutate_if(is.numeric, asin) %>%
  dplyr::filter(!DNA.Source == "Bacterioplankton") %>%
  dplyr::filter(!Organism == "Water Control", !Organism == "State Water", !Organism == "Start Water")

#remove original dataframe
rm(Cr_Data)

#making abundance only matrix and saving columns with names/metadata into dinames
microb_f <- as.matrix(microb_f_abun %>%
                        dplyr::select(-c(Label:DNA.Source)),
                                      dinames = list(paste("Label", 1:31, sep = ",")))
# Run NmDS, provide community x species matrix, and set dimensions to 2
# outputs each iteration until a solution is reached (minimize stress after reconfiguration of points in 2D)

set.seed(1234)
microb_NMDS <- metaMDS(microb_f)

# produce Shepard plot
# large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. 
stressplot(microb_NMDS)

# plot the NMDS
# open circles are communities, red crosses are species
plot(microb_NMDS)

ordiplot(microb_NMDS, type = "n")
orditorp(microb_NMDS, display = "sites", air = 0.25)


# create a plot with convex hulls connecting vertices of the points made by a treatment

treatment <- c(microb_f_abun$Organism)
samples <- c(microb_f_abun$Label)
##Different Plotting options. Comment in to plot

# ordiplot(microb_NMDS, type = "n",
#          main = "NMDS of Incubation Vessels")
# ordihull(microb_NMDS, groups = treatment, draw = "polygon", label = F)
# #orditorp(microb_NMDS, display = "sites", labels = samples, air = 0.05)
# orditorp(microb_NMDS, display = "sites", air = 0.2)
# 
# 
# # spider plot
# ordiplot(microb_NMDS, type = "n",
#          main = "NMDS of Incubation Vessels")
# orditorp(microb_NMDS, display = "sites", air = 0.2)
# ordispider(microb_NMDS, groups = treatment)
# 
# 
# # ellipse plot
# ordiplot(microb_NMDS, type = "n",
#          main = "NMDS of Incubation Vessels")
# orditorp(microb_NMDS, display = "sites", air = 0.2)
# ordiellipse(microb_NMDS, groups = treatment)
# 
# 
# 
# # MST plot
# 
# dist_microb <- vegdist(microb)
# clust_microb <- hclust(dist_microb, method = "complete")
# 
# 
# ordiplot(microb_NMDS, type = "n",
#          main = "NMDS of Incubation Vessels")
# orditorp(microb_NMDS, display = "sites", air = 0.2)
# ordicluster(microb_NMDS, cluster = clust_microb)


# plot NMDS output in ggplot

df.scores <- as.data.frame(scores(microb_NMDS)) #Using the scores function from vegan to extract the site scores and convert to a data.frame

df.scores$site <- rownames(df.scores) # create a column of site names, from the rownames of data.scores

# add grouping variables
df.scores$sample <- microb_f_abun$Label
df.scores$group <- microb_f_abun$Water
df.scores$DNA_source <- microb_f_abun$DNA.Source
df.scores$Organism <- microb_f_abun$Organism


species.scores <- as.data.frame(scores(microb_NMDS, "species")) #Using the scores function from vegan to extract the species scores and convert to a data.frame

species.scores$species <- rownames(species.scores)



# plot NMDS

ggplot() +
  geom_point(data = df.scores, aes(x = NMDS1, y = NMDS2), alpha = 0.5) +
  coord_equal() +
  scale_shape_manual(values = c(0, 19)) +
  theme_bw()



# write csv for Zach to plot in JMP
NMDS_scores <- as.data.frame(scores(microb_NMDS, display = c("sites")))


# add grouping variables
#inputs back into your new table the metadata
NMDS_scores$Sample <- microb_f_abun$Label
NMDS_scores$Water_Treatment <- microb_f_abun$Water.Treatment
NMDS_scores$DNA_source <- microb_f_abun$DNA.Source
NMDS_scores$Organism <- microb_f_abun$Organism



write.csv(NMDS_scores, "Crusade_Micro_NMDS_scores.csv")