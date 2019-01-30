library(tidyverse)
library(cluster)
library(factoextra)
library(dendextend)
library(data.table)
library(RColorBrewer)


##      KEEPING THIS SCRIPT FOR POSTERITY BUT IT IS NOT TOTALLY COMPLETE AND DOES NOT LOOK GOOD   ##
##      JMP IS THE PLACE I SHOULD BE MAKING THE HEATMAPS IN; NOT R.     ##


##        METHODS USED    ##

##  AGGLOMERATIVE CLUSTERING: 
##    It’s also known as AGNES (Agglomerative Nesting). It works in a bottom-up manner. 
##    That is, each object is initially considered as a single-element cluster (leaf). 
##    At each step of the algorithm, the two clusters that are the most similar are combined into a new bigger cluster (nodes).
##    This procedure is iterated until all points are member of just one single big cluster (root) (see figure below). 
##    The result is a tree which can be plotted as a dendrogram.

##  WARD'S MINIMUM VARIANCE METHOD:
##    It minimizes the total within-cluster variance. At each step the pair of clusters with minimum between-cluster distance are merged
##    See bellow for method on how to select correct variance method (probably ward)

##      DATA PREPARATION    ##

##  1)  Rows are observations (individuals) and columns are variables
##  2)  Any missing value in the data must be removed or estimated.
##  3)  The data must be standardized (i.e., scaled) to make variables comparable. 

##  FOR CRUSADE I USED SUM TO MAKE THE TABLE FOR THIS SCRIPT  ##

cluster_raw_data <- read_csv("tp50_RA_cleaned.dat")

#Cleaning some MetaData and normalizing arcsin(sqrt)
normalized_cluster <- cluster_raw_data%>%
  dplyr::filter(!Label %like% "ST%")%>% ##  Successfully removed Start water using the "%like%" opperator
  mutate_if(is.numeric, sqrt)%>% #  Normalizing data (arcsin(sqrt(x)))
  mutate_if(is.numeric, asin)

Scnorm_cluster <- normalized_cluster%>%
  select(-c(1))%>%
  scale()

##      PICKING THE BEST METHOD TO USE    ##

## You can check for the high AC of different methods using the below code:
# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(Scnorm_cluster, method = x)$ac
}

map_dbl(m, ac)


##  compute differences with agnes
agnes_cluster <- agnes(Scnorm_cluster, method = "ward") ## I USED WARD HERE BECAUSE IT HAD THE HIGHEST AGG COEFFICIENT

## Agglomerative coefficient
agnes_cluster$ac

pltree(agnes_cluster, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 

hv <-heatmap(Scnorm_cluster, col = brewer.pal(n=11, name = "RdBu"), scale = "column", margins = c(5,10),
             xlab = "Test", ylab = "test2",
             main = "heatmap test")
