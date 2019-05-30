#AYESHA HARGEY
#3650393
#CCA

#Load libraries
library(tidyverse)
library(ade4)
library(vegan)
library(MASS)



#Load data
env <- read_csv("~/Quant/Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsEnv.csv")
View(env)
spe <- read_csv("~/Quant/Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpe.csv")
View(spe)
spa <- read_csv("~/Quant/Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpa.csv")
View(spe)

#Remove empty site 
spe <- spe[-8,-1]
env <- env[-8, ]
spa <- spa[-8, ]

#making column into a separate dataframe
dfs <- env[ ,1]

#removing the column from the original
env <- env[ ,-1]

#landscape feature/topography 
envtopo <- env[ ,c(1:2)]
names(envtopo)

#chemical variables only
envchem <- env[, c(4:10)]
names(envchem)

#Hellinger-transform the species dataset
#Suited to species abundance data. Low weight to species with low counts and many zeroes.
spe.hel <- decostand(spe, "hellinger")
spe.hel 

#Canonical correspondance analysis
#'~' as a function of
spe.cca <- cca(spe.hel ~., env) #dot means all the column names within the dataset
#Call: cca(formula = spe.hel ~ alt + slo + flo + pH + har + pho + nit +
#amm + oxy + bod, data = env)
summary(spe.cca)

#The community structure that we are going to determine is influenced by
#the nature of the environmental variables
#in which the species find themselves at the particular site
#i.e a CONSTRAINT

#sum of unconstrained axes makes up the unconstrained inertia

#in order to calculate proportion of constrained 

#Variance explained by full mode
sum(spe.cca$CCA$eig) / spe.cca$tot.chi * 100
#Adjusted R
(spe.cca.R2a <- RsquareAdj(spe.cca)$adj.r.squared)

#Variance explained by full model
#Null hypothesis: The model cannot explain anything about the data

anova(spe.cca) #ANOVA DONE to see if it can determine if it can probe the hypothesis
#P is 0.001
#Null hypothesis is rejected
#Hypothesis is accepted

#which of these axes can explain the observed variation?
#this anova is done
anova(spe.cca, by = "axis") #CCA1, CCA2, CCA3
#this analysis of variance determines which environmental factor has the most impact
#'by' is the terms used in the environmental data
#which of those explains the most variation
spe.cca.axis.test <- anova(spe.cca, by = "term")
#altitude, slope, hardness, phosphorous and oxygen

#3 different anovas answering 3 different facets
#is the model able to explain variation
#which axes explains the variation
#which of the terms explain the variation

#only significant environmental values
#selects rows that are less than 0.05
spe.cca.ax <- which(spe.cca.axis.test[ ,4] < 0.05)
env2 <- env[, spe.cca.ax] #of the total environmental dataset
#choose only the columns which are present in 'spe.cca.ax' i.e the significant columns

#cca on the edited dataset
spe.cca2 <- cca(spe.hel ~., env2)
summary(spe.cca2)
#same anova transformations as previously
anova(spe.cca2)
anova(spe.cca2, by = "axis")
anova(spe.cca2, by = "term") #hardness is now excluded

#variance inflation factors
vif.cca(spe.cca2) #determines what's co-linear so that the arrows all point in different directions
#becomes unambiguous 
#if it's two things pointing at the same direction it is ambiguous because
#you are unable to discern what is what
#runs the vif on the cca reduced model
#if any of these values are larger than 10, they would be colinear 

#ordiplots 
#Scaling 1

#plot(spe.cca2, scaling = "sites", choices 1:2, type = "none",
     main = "CCA (sites scaling)" 
