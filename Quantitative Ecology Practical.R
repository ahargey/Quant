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

summary(env)

# Hellinger pre-transformation of the species data
spe_h <- spe %>% 
  select(-X1)
spe_h <- decostand(spe_h, "hellinger")
