#AYESHA HARGEY
#3650393
#PCA - Mite data

#Load libraries
library(tidyverse)
library(vegan)

mites <- read_delim("Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/mite.txt", 
                   "\t", escape_double = FALSE, trim_ws = TRUE)
View(mites)

mites.pca <- rda(mites, scale = TRUE)
mites.pca
summary(mites.pca)

eig <- sum(mites.pca$CA$eig) #inertia
#inertia tells us about spread of variation
#sum of the diagonal values

mites.pca.cor <- cor(mites)
mites.pca.cor

(mites.pca$CA$eig[1]/sum(mites.pca$CA$eig))*100

par(mfrow=c(1,2))
#biplot
biplot(mites.pca, scaling = 2, choices = c(1, 2), main = "PCA - Scaling 1")
biplot(mites.pca, scaling = 1, choices = c(1,2), main = "PCA - Scaling 2")

#cleanplot
cleanplot.pca(mites.pca, scaling=1, mar.percent=0.06)
cleanplot.pca(mites.pca, scaling=2, mar.percent=0.06)
