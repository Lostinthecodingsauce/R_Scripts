##This is a subset of my whole pipeline
# Running Two-Way ANOVAs OTU -----------------------------------------------------
cr_wdf <- read_csv("Crusade_wdf.dat ")
only_cca <- cr_wdf%>%
  dplyr::filter(!Organism == "Calcium Carbonate Control",
                !Organism == "Water Control",
                !Organism == "Start Water")

##  subset the data by bacterioplankton
occa_bp <- only_cca%>%
  dplyr::filter(DNA.Source == "Bacterioplankton")

##Running the two way anova for Bacterioplankton Organism*Water

# This line makes the names of the rows which will be added into the pvalue table
twoway_anova_rows <- c("Organism", "Water", "Organism*Water")

# This is the line which applies the anova over every single clade in the bacterioplankton and runs the two way anova
## This is R's faster version of a for loop. It will apply the function of summary(aov()) across all columns in the
## ONLY CCA BACTERIOPLANKTON (occa_bp) matrix which are relative abundance data
## The [[1]][1:3, 'Pr(>F)'] is pulling the F_values out of the summary tables and making a new summary data matrix called aov_bp
aov_bp <- sapply(occa_bp[8:667], function(x) summary(aov(x ~ occa_bp[["Organism"]]*occa_bp[["Water"]]))[[1]][1:3,'Pr(>F)'])


# It comes out as a list so it has to be first converted to a data frame before we can add in the test names as a column
anova_bp <- as.data.frame(aov_bp)
anova_bp$Anova_test <- cbind(twoway_anova_rows)

# Now the data is made Tidy and we filter to only significant values
anova_bp_tidy <- anova_bp%>%
  gather(Clade, F_value, 1:660)

anova_bp_tidy$Water <- "Bacterioplankton"

write_csv(anova_bp_tidy, "CRUSADE_TWANOVA_Bacterioplankton.dat")