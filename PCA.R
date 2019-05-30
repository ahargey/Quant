#AYESHA HARGEY
#3650393
#PCA

#Load libraries
library(tidyverse)
library(vegan)

spe <- DoubsSpe
spa <- DoubsSpa
env <- DoubsEnv

env <- read_csv("~/Quant/Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsEnv.csv")
View(env)
spe <- read_csv("~/Quant/Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpe.csv")
View(spe)

spe <- spe[-8]

source("evplot.R")
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")

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

spe.pca <- rda(spe, scale = TRUE)
spe.pca
summary(spe.pca)

eig <- sum(spe.pca$CA$eig) #inertia
#inertia tells us about spread of variation
#sum of the diagonal values

spe.pca.cor <- cor(spe)
spe.pca.cor

(spe.pca$CA$eig[1]/sum(spe.pca$CA$eig))*100

par(mfrow=c(1,2))
#biplot
biplot(spe.pca, scaling = 2, choices = c(1, 2), main = "PCA - Scaling 1")
biplot(spe.pca, scaling = 1, choices = c(1,2), main = "PCA - Scaling 2")

#cleanplot
cleanplot.pca(spe.pca, scaling=1, mar.percent=0.06)
cleanplot.pca(spe.pca, scaling=2, mar.percent=0.06)

# Site 12, 13, 18 have similar species
# Satr, Phph, Babl, Thth, Cogo, Teso
# 
# Site 19, site 20, site 21 have very similar species
# Alal, Gyce, Pato, Lete, Baba, Pefl, Gogo, Soce
# 
# Many sites such as site 1, 8, 23 have no species 
# 
# PC1 is responsible for 61.34% of the species spread (this was determined by eigenvalues divided by the inertia)
# 
# We are using two scalings to reiterate the values, and to see alternative interpretations. An alternative arrangement of species. 
# 
# Scaling 2 also suggests a clustering on the sites 20-29 with the exception of 23, 24, 25. 
# 2 ordination axes are meaningful
# Heuritistic approach  


spe.pca1 <- scores(spe.pca, display="species", scaling=1)
cleanplot.pca(spe.pca1, scaling=1, mar.percent=0.06)

# Hellinger pre-transformation of the species data
Doubsspe.h <- decostand(spe, "hellinger")
(spe.h.pca <- rda(Doubsspe.h))
spe.h.pca 

# Plot eigenvalues and % of variance for each axis
ev <- spe.h.pca$CA$eig
dev.new(title="PCA eigenvalues")
evplot(ev)

# PCA biplots
spe.pca.sc1 <- scores(spe.h.pca, display="species", scaling=1)
spe.pca.sc2 <- scores(spe.h.pca, display="species", scaling=2)


dev.new(title="PCA on fish species", width=12, height=6)
par(mfrow=c(1,2))



