# QUANTITATIVE ECOLOGY
# Ayesha Hargey
# 3650393
# Textbook Data 
# 
#Load libraries
library(tidyverse)
library(ade4)
library(adespatial)
library(vegan)
library(gclus)
library(cluster)
library(FD)

source("coldiss.R")
source("panelutils.R")

env <- read_csv("Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsEnv.csv")
spa <- read_csv("Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpa.csv")
spe <- read_csv("Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpe.csv")

#remove empty site 8
spe <- spe %>% 
  select(-8)
env <- env %>% 
  select(-8)
spa <- spa %>% 
  select(-8)

dim(spe_db)

#env
glimpse(env) #overall preview of data, shows every column
head(env, n = 9) #first nine rows
tail(env, n = 9) #last nine rows
nrow(env) #number of row
ncol(env) #number of columns
any(is.na(env)) #is there any missing data?
summary(env) #summary of the data according to quartiles and min/max values

#spa
glimpse(spa) #overall preview of data, shows every column
head(spa, n = 9) #first nine rows
tail(spa, n = 9) #last nine rows
nrow(spa) #number of row
ncol(spa) #number of columns
any(is.na(spa)) #is there any missing data?
summary(spa) #summary of the data according to quartiles and min/max values

#spe
glimpse(spe) #overall preview of data, shows every column
head(spe, n = 9) #first nine rows
tail(spe, n = 9) #last nine rows
nrow(spe) #number of row
ncol(spe) #number of columns
any(is.na(spe)) #is there any missing data?
summary(spe) #summary of the data according to quartiles and min/max values

#bray curtis is being used because the species data contains the known abundance
#as opposed to jaccard which is when it's a simple presense/abscence

#Bray-Curtis dissimilarity on raw species data
spe_db <- vegdist(spe, diag = "TRUE", upper = "TRUE")
head(spe_db)

#Bray Curtis on log transformed data
spe_log <- vegdist(log1p(spe))
head(spe_log)

#Chord distance matrix
spe.dc <- dist.ldc(spe, "chord")

#Jaccard
spe.dj <- vegdist(spe, "jac", binary = TRUE)
head(spe.dj)
head(sqrt(spe.dj))

# Percentage difference (aka Bray-Curtis) dissimilarity matrix on
# raw species abundance data
coldiss(spe_db, byrank = FALSE, diag = TRUE)
# Same but on log-transformed data
coldiss(spe_log, byrank = FALSE, diag = TRUE)

hist(spe_db)

spe_db_dc <- decostand(spe_db, "jac", binary = TRUE)
