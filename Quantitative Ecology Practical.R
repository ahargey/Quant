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
sum(spe == 0) #number of absences
sum(spe == 0) / (nrow(spe) * ncol(spe)) #proportion of zeros

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

env <- env %>% 
  select(-X1) #removing name column

spe <- spe %>% 
  select(-X1) #removing name column

#spatial plot of quadrats
spa_plot <- ggplot(data = spa, aes(x = X, y = Y)) + #assigning a name to plot 
  geom_path(data = spa, aes(label = X1), size = 0.5, colour = "blue") + #path
  geom_text(data = spa, aes(label = X1), size = 4, colour = "red") + #quadrat names
  annotate("text", label = "Entrance", x = 1, y = 0.5, size = 4.0, angle = 0, colour = "purple") +
  annotate("text", label = "Mountain", x = 7, y = 4.5, size = 4.0, angle = 0, colour = "green") +
  annotate("text", label = "Slope", x = 7, y = 7.8, size = 5.0, angle = 45, colour = "orange") +
  labs(x = "X coordinate (m)", y = "Y coordinate (m)", title = "Quadrat Locations") +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
spa_plot

#CCA analysis is done

spe.cca <- cca(spe ~ ., env) #the period uses all columns in the environmental dataframe
summary(spe.cca)

# The adjusted R2 --- the variance explained by the constrained axes
spe.cca.R2a <- RsquareAdj(spe.cca)$adj.r.squared
spe.cca.R2a

# Variance explained by full model
sum(spe.cca$CCA$eig) / spe.cca$tot.chi * 100

(spe.cca$CCA$eig[1]/sum(spe.cca$CCA$eig))*100 #the variation explained by the first scaling
(spe.cca$CCA$eig[2]/sum(spe.cca$CCA$eig))*100 #the variation explained by the second scaling

#Testing to see if this analysis is significant overall
anova(spe.cca, permutations = how(nperm = 999))
#It is not, p > 0.05
#Testing to see if the analysis is significant based on the environmental factors
anova(spe.cca, by = "axis", permutations = how(nperm = 999))
#It is not, p > 0.05

vif.cca(spe.cca) #no values over 10 so no redundant contraints 

#Graphs

par(mfrow = c(1, 2)) #parameters of graphs 
#Scaling 1: species are scored to the relative eigenvalues,
plot(spe.cca,
     scaling = 1,
     display = c("sp", "lc", "cn"),
     main = "Triplot CCA species and environmental data - scaling 1")

#Scaling 2: site scores scaled to the relative eigenvalues
plot(spe.cca,
     display = c("sp", "lc", "cn"),
     main = "Triplot CCA species and environmental data - scaling 2")