#AYESHA HARGEY
#3650393
#PCA

#Load libraries
library(tidyverse)
library(vegan)

env <- read_csv("~/Quant/Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsEnv.csv")
View(env)

#PCA
env.pca <- rda(env, scale = TRUE)
env.pca

cor.env <- cor(env) #unconstrained and inertia is the same because it's an unconstrained

summary(env.pca)

biplot(env.pca, scaling = 2, choices = c(1, 2), main = "PCA - scaling 1")
#shortest arrow has the lowest value eg. in this case, pca

#mites
#dunes
#do this on the doubs species data

(env.pca$CA$eig[1]/sum(env.pca$CA$eig))*100 #tells the percentage influence the first PCA has


  