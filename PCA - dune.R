#AYESHA HARGEY
#3650393
#PCA - dune

??dune

dune <- data("dune")
dune_env <- data(dune.env)

dune.pca <- rda(dune, scale = TRUE)
dune.pca
summary(dune.pca)

eig <- sum(dune.pca$CA$eig) #inertia
#inertia tells us about spread of variation
#sum of the diagonal values

dune.pca.cor <- cor(dune)
dune.pca.cor

(dune.pca$CA$eig[1]/sum(dune.pca$CA$eig))*100

par(mfrow=c(1,2))
#biplot
biplot(dune.pca, scaling = 2, choices = c(1, 2), main = "PCA - Scaling 1")
biplot(dune.pca, scaling = 1, choices = c(1,2), main = "PCA - Scaling 2")

#cleanplot
cleanplot.pca(dune.pca, scaling=1, mar.percent=0.06)
cleanplot.pca(dune.pca, scaling=2, mar.percent=0.06)
