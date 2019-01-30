# Load Packages
library(FactoMineR)
library(factoextra)
library(ggord)

#Transform data and create dataframe with categories
data.filt <- read.csv("CRUSADE_RA_NT.dat") %>%
  dplyr::filter(DNA.Source == "Bacterioplankton") %>%
  dplyr::filter(!Organism == "Water Control", !Organism == "State Water", !Organism == "Start Water")

data.asin<-data.filt%>%
  mutate_if(is.numeric, sqrt)%>% #Normalizing data (arcsin(sqrt(x)))
  mutate_if(is.numeric, asin) %>%
  dplyr::select(-c(Label,Organism,Water,Timepoint,DNA.Source))



#  Run PCA
res.pca<-PCA(data.asin)
# Check out your eigenvalues
eig.val<-get_eigenvalue(res.pca)
# Plot eigenvalues
fviz_eig(res.pca)
#Inspect the dimensions
res.desc<-dimdesc(res.pca,axes=c(1,2),proba=.05)
#Check out the variables and their influence
fviz_pca_var(res.pca, col.var='cos2', gradient.cols=c("Red","Orange","Blue"), repel = TRUE,select.var= list(cos2 = 10))
#Check out the Individuals and add ellipes
fviz_pca_ind(res.pca, geom.ind='point', palette=c('Red', 'Blue'), addEllipses=TRUE, mean.point=FALSE)
# Plot Variables and Individuals together only plotting the top 10 Individuals influencing vectors
plot<-fviz_pca_biplot(res.pca, geom.ind=c('point',"text"),
                      select.var=list(contrib=15),
                      mean.point=FALSE,
                      pch=19,
                      ggtheme = theme_bw(),
                      col.var = "black",
                      repel=TRUE,
                      vec_ext = 0, txt = NULL, arrow = 0)
#ggord is good for 95% confidence intervals for only one filter. Not trying to look at water and organism for example
ord<-ggord(res.pca, data.filt$Organism, grp_in = data.filt$Water, vec_ext = 0, txt = NULL, arrow = 0)
ord
#ggsave(filename = "CCA.PCA.HPC.pdf", plot = ord)
