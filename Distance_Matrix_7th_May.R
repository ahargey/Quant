#QUANTITATIVE ECOLOGY
#Ayesha Hargey
#3650393
#Distance Matrix

library(tidyverse)
library(vegan)
?vegdist

data("varespec")
data("varechem")
vegist <- vegdist(varespec, method = "bray")
vegist
