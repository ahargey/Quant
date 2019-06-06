#AYESHA HARGEY
#3650393
#Quantitative Ecology Practical

#Load libraries
library(tidyverse)
library(ade4)
library(adespatial)
library(vegan)
library(gclus)
library(cluster)
library(FD)
library(ggpubr)

#Loading the data 
#Environmental data had to be edited in Excel as it used comma decimals and not point
env <- read_delim("Final Env Comma Corrected.csv", 
                  ";", escape_double = FALSE, trim_ws = TRUE)
spa <- read_delim("Final Spatial.csv", ";", 
                  escape_double = FALSE, trim_ws = TRUE)
spe <- read_delim("Final Species.csv", ";", 
                  escape_double = FALSE, trim_ws = TRUE)

spe <- spe %>% mutate_all(funs(replace_na(.,0))) #turns NA values into 0 for the species table

#Exploring data
glimpse(spe) #overall preview of data, shows every column
head(spe) #first six rows
tail(spe) #last six rows
nrow(spe) #number of row
ncol(spe) #number of columns
summary(spe) #gives a summary of the mean, median, quartiles and min/max values

glimpse(env) 
head(env) 
tail(env) 
nrow(env) 
ncol(env) 
summary(env)

#Removing empty sites
spe <- spe[c(-5, -10, -16, -27), ]
env <- env[c(-5, -10, -16, -27), ]
spa <- spa[c(-5, -10, -16, -27), ]

env_test <- env %>% 
  select(-X1)

# Hellinger pre-transformation of the species data
spe_h <- spe %>% 
  select(-X1) #removing the named column
spe_h <- decostand(spe_h, "hellinger")
#Hellinger transformation done as this is species abundance data.
#It gives low weights to species with minimal observations.

#CCA
spe_cca <- cca(spe_h ~ env_test$dfs)
summary(spe_cca)

# The adjusted R2 --- the variance explained by the constrained axes
(spe_cca_R2a <- RsquareAdj(spe_cca)$adj.r.squared)

# Variance explained by full model
sum(spe_cca$CCA$eig) / spe_cca$tot.chi * 100

anova(spe_cca)

DoubsSpe.cca.axis.test <- anova(spe_cca, by = "term")

#spatial plot
spa_plot <- ggplot(data = spa, aes(x = X, y = Y)) +
  geom_path(data = spa, aes(label = X1), size = 0.5, colour = "blue") +
  geom_text(data = spa, aes(label = X1), size = 4, colour = "red") +
  annotate("text", label = "Entrance", x = 1, y = 0.5, size = 4.0, angle = 0, colour = "purple") +
  labs(x = "x coordinate (m)", y = "y coordinate (m)", title = "Quadrat Locations") +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
spa_plot
