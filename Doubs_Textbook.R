# QUANTITATIVE ECOLOGY
# Ayesha Hargey
# 3650393
# Textbook Data 
# 
#Load libraries
library(tidyverse)

env <- read_csv("Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsEnv.csv")
spa <- read_csv("Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpa.csv")
spe <- read_csv("Textbook/NEwR_1st_edition_R332/NEwR_updated_material_R332-NEwR_ed1/DoubsSpe.csv")

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