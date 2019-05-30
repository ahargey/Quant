#AYESHA HARGEY
#3650393
#PCA

#Load libraries
library(tidyverse)
library(vegan)
library(ggpubr)

spe <- DoubsSpe
spa <- DoubsSpa
env <- DoubsEnv

env <- read_csv("~/Quant/Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsEnv.csv")
View(env)
spe <- read_csv("~/Quant/Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpe.csv")
View(spe)
spa <- read_csv("Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpa.csv")
View(spa)

#environment data
env.pca <- rda(env, scale = TRUE)
env.pca

cor.env <- cor(env) #unconstrained and inertia is the same because it's an unconstrained

summary(env.pca)

biplot(env.pca, scaling = 2, choices = c(1, 2), main = "PCA - scaling 2")

#species data
spe <- spe[-8] #cut out site 8

source("evplot.R")
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")

#PCA
spe_pca <- rda(spe, scale = TRUE)
spe_pca
summary(spe_pca) #Anan is highest (0.97) 
#total amount of species gives us total inertia 

eig <- sum(spe_pca$CA$eig) #inertia
eig
#inertia tells us about spread of variation
#sum of the diagonal values

spe_pca_cor <- cor(spe)
spe_pca_cor #association

(spe_pca$CA$eig[1]/sum(spe_pca$CA$eig))*100 #PCA1
(spe_pca$CA$eig[2]/sum(spe_pca$CA$eig))*100 #PCA2

par(mfrow=c(1,2))
#biplot
biplot(spe_pca, scaling = 2, choices = c(1, 2), main = "PCA - Scaling 2")
biplot(spe_pca, scaling = 1, choices = c(1,2), main = "PCA - Scaling 1")

#cleanplot
cleanplot_pca(spe_pca, scaling=1, mar.percent=0.06)
cleanplot_pca(spe_pca, scaling=2, mar.percent=0.06)

#environmental + species
one <- biplot(spe_pca, scaling = 2, choices = c(1, 2), main = "PCA - Scaling 2")
two <- biplot(env.pca, scaling = 2, choices = c(1, 2), main = "PCA - Scaling 1")

# Site 12, 13, 18 have similar species
# Satr, Phph, Babl, Thth, Cogo, Teso
# 
# Site 19, site 20, site 21 have very similar species
# Alal, Gyce, Pato, Lete, Baba, Pefl, Gogo, Soce
#
# Two major clusters of species, one on the left, one on the right. 
# Negative association between these two areas
#
# PC1 is responsible for 60.84% of the species spread (this was determined by eigenvalues divided by the inertia)
# 
# We are using two scalings to reiterate the values, and to see alternative interpretations. An alternative arrangement of species. 
# 
# Scaling 2 also suggests a clustering on the sites 20-29 with the exception of 23, 24, 25. 
# 2 ordination axes are meaningful, heuritistic approach
# There were four main clusters in environmental data, most notably of oxygen
# Suggesting there could be an oxygen gradient 

biplot(env.pca, scaling = 2, choices = c(1, 2), main = "PCA - scaling 1")
biplot(spe_pca, scaling = 1, choices = c(1, 2), main = "PCA - Scaling 1")

spa$X1 <- as.numeric(spa$X1)


plot(spa)
spa_plot <- ggplot(data = spa, aes(x = X, y = Y)) +
  geom_path(data = spa, aes(label = X1), size = 0.5, colour = "blue") +
  geom_text(data = spa, aes(label = X1), size = 4, colour = "red") +
  #geom_line() +
  #lines() +
  annotate("text", label = "Downstream", x = 5, y = 39, size = 4.0, angle = 0, colour = "red") +
  annotate("text", label = "Upstream", x = 90, y = 10, size = 4.0, angle = 0, colour = "red") +
  labs(x = "x coordinate (km)", y = "y coordinate (km)", title = "Site Locations") +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))
spa_plot

ggarrange(one, two, spa_plot, ncol = 2, nrow = 2)
