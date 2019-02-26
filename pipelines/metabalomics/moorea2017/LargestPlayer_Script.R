#RR3 metabolomes largest player script
#Made by Zach Quinlan 9/16/18
## A cleaner version of this script is in the Moorea2017Metabalomics.R pipeline under the Largest Player section.


#load in the libraries needed for this
library(dplyr)
library(tidyr) #I believe this one is needed but cannot remember. dplyr is almost exclusively used.


#Before importing into R I used jmp (could use excel or whatever) to take the mean values for each organism, calculated max values,
# and exported the dataframe as a .dat (csv). Taking means in r is annoying...
rr3_df <- read.csv("RR3_SignificantFeatures.dat") #typically I like .csv or .dat files

#getting rid of extra columns not needed for this data processing
rr3_means <- rr3_df %>%
  dplyr::select(-c(DPL1:DTR3)) %>% #This gets rid of anything which is not originally in 
  dplyr::arrange(Feature) #have to sort by feature because some colnames bellow automatically does this and so you need everything to line up


##Summarizing the four column's means into a new max column
##This part works however coerces the data into a data table not data frame which makes epxorting a pain.
##I just calculated mean using jmp instead of this. But you can probably pretty easily edit it to make it work.
#rr3_means$Max <- group_by(rr3_means, Feature) %>%
 #    summarize(max(Mean.PL. , Mean.DT. , Mean.CC. , Mean.PV. , Meant.TR.))

#Making data frame for the integration of max column names
#For colnames to work you have to only have the columns you want in it.
rr3_meansonly <-rr3_means %>%
  dplyr::select(-c(1:2)) #remove "feature" and "Cluster columns"
  dplyr::select(-c(Max)) #remove the Max column for this

rr3_means$MaxNames <- colnames(rr3_meansonly)[max.col(rr3_meansonly, ties.method="first")] #Finds column names for max values and makes a new column in old data frame

#writing csv
write.csv(rr3_means, file = "RR3_SignificantFeatures_Max.dat")
