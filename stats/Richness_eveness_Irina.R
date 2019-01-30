# Shannon and Pielous
library(vegan) #This library is for Shannon and Pielous
library(dplyr) #Organization

# read in data. Make sure:
# 1) Microbial community members as your response variables
# 2) Samples are your observations
original_df <- read.csv("Pn_Ex2_16S_counts_transposed.dat")

# #transform data arc(sqrt()) 
# #**I dont think this is entirely needed. But check your data and look at transformed as well.
# microb_df <- microb_NT%>%
#   mutate_if(is.numeric, sqrt)%>%
#   mutate_if(is.numeric, asin)

# edit so only microb reads
microb_abun <- original_df %>%
  dplyr::select(-c(Label))

# Richness index (Shannan-Weaver Index)
# Makes a new column in original data frame with Shannons Richness
original_df$Shannon <-diversity(microb_abun)

# Eveness index (Pielous Eveness Index)
# Makes a new column in original data frame with Pielous Eveness
# equation for Pielous == (Shannons/log(number of species))
original_df$Pielous <-diversity(microb_abun)/log(specnumber(microb_abun))


#Writing file
write.csv(original_df, "Pn_Ex2_16s_counts_transposed_Shannon_Pielous.dat")
