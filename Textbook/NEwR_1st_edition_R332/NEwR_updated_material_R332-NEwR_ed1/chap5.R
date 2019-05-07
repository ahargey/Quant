################################################################################
### CHAPTER 5: ORDINATION IN REDUCED SPACE
###
### Online supporting material for: 
### Borcard D., Gillet F. & Legendre P. Numerical Ecology with R, Springer, 2011
### Updated for R 3.3.2 : 3 February 2017
################################################################################

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(FactoMineR)

# Load additionnal functions
# (files must be in the working directory)
source("evplot.R")
source("cleanplot.pca.R")
source("PCA.newr.R")
source("CA.newr.R")

# Import the data from CSV files
# (files must be in the working directory)

# Doubs fish data

spe <- DoubsSpe
spa <- DoubsSpa
env <- DoubsEnv



# Remove empty site 8
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]

# decostand func
env_std <- decostand(env, method = "standardize")

cor(env)
cor(env_std)
# calculating the mean of each of the columns:
round(apply(env_std, 2, FUN = mean, 2)) # The first 2 in this formula is to calculate the mean for the column. if the value was 1 then it calculates the rows 
round(apply(env_std, 2, FUN = sd, 2))

env_pca <-  rda(DoubsEnv, scale = TRUE)
sum(env_pca$CA$eig) # Calculating the sum
# What does the correlation matrix mean?
# when you compare a variable with itself it equals 1
# what are eignen values - 
summary(env_pca)
# Species scores - eigan vectors - correltion between unscaled value and the new value
#

# Oribatid mite data
mite <- read.table("mite.txt")
mite.env <- read.table("mite_env.txt")


# PCA on the full environmental dataset
# *************************************

# A reminder of the content of the env dataset
summary(DoubsEnv)
?rda
# PCA based on a correlation matrix
# Argument scale=TRUE calls for a standardization of the variables
# Rda different from pca 
env.pca <- rda(DoubsEnv, scale=TRUE) ## Scale : Mean 0 and SD of 1- dont need decostand analysis when using the scale function
env.pca  # pC1 Has most amount of power when compared to all of the variables
summary(env.pca) # Default scaling 2
summary(env.pca, scaling=1)

# Examine and plot partial results from PCA output
?cca.object   # Explains how an ordination object
              # produced by vegan isstructured and how to
              # extract its results.

# Eigenvalues
(ev <- env.pca$CA$eig)
# Apply Kaiser-Guttman criterion to select axes
ev[ev > mean(ev)]

# Broken stick model
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n)
{
	bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
}
bsm$p <- 100*bsm$p/n
bsm

# Plot eigenvalues and % of variance for each axis
dev.new(title="PCA eigenvalues")
par(mfrow=c(2,1))
barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")		# average eigenvalue
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE, 
	main="% variance", col=c("bisque",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"), 
	pch=15, col=c("bisque",2), bty="n")

# Same plots using a single function: evplot()
# Plot eigenvalues and % of variance for each axis
evplot(ev)


# Two PCA biplots: scaling 1 and scaling 2
# ****************************************

# Plots using biplot.rda
dev.new(width=12, height=6, title="PCA biplots - environmental variables - biplot.rda")
par(mfrow=c(1,2))
biplot(env.pca, scaling=1, main="PCA - scaling 1")
biplot(env.pca, main="PCA - scaling 2")  # Default scaling 2
# Smaller the angle between the plot is the closer they are related to eachother



# Plots using cleanplot.pca - NEW VERSION OF THIS FUNCTION wich changed syntax
# A rectangular graphic window is needed for the two plots
dev.new(width=12, height=6, title="PCA biplots - environmental variables - cleanplot.pca")
par(mfrow=c(1,2))
cleanplot.pca(env.pca, scaling=1, mar.percent=0.08)
cleanplot.pca(env.pca, scaling=2, mar.percent=0.04)



# Combining clustering and ordination results
# *******************************************

# Clustering the objects using the environmental data: Euclidean 
# distance after standardizing the variables, followed by Ward
# clustering
env.w <- hclust(dist(scale(env)), "ward.D")
# Cut the dendrogram to yield 4 groups
gr <- cutree(env.w, k=4)
grl <- levels(factor(gr))

# Extract the site scores, scaling 1
sit.sc1 <- scores(env.pca, display="wa", scaling=1)

# Plot the sites with cluster symbols and colours (scaling 1)
dev.new(title="Ordination and clustering")
p <- plot(env.pca, display="wa", scaling=1, type="n", 
	main="PCA correlation + clusters")
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
for (i in 1:length(grl))
{
	points(sit.sc1[gr==i,], pch=(14+i), cex=2, col=i+1)
}
text(sit.sc1, row.names(env), cex=0.7, pos=3)
# Add the dendrogram
ordicluster(p, env.w, col="dark grey")
# Add legend interactively
legend(locator(1), paste("Cluster", c(1:length(grl))), pch=14+c(1:length(grl)), 
	col=1+c(1:length(grl)), pt.cex=2)


Find PCA reading explaination



# 
# PCA on species data(apply jaccard bindary = true)
# rda and say scaling = true

spec_jac <- round(vegdist(DoubsEnv, method = "jaccard"), 2)
spec.jac <- rda(spec_jac, scale=TRUE) ## Scale : Mean 0 and SD of 1- dont need decostand analysis when using the scale function
spec.jac 

# R-vegan PCA
# hypothesis testing 
# summary pca output and c how numbers relate
# total inertia(eigen value)
# Basic analysis within R discriptions and how to do PCA and understanding RDA in R 



# Plotting environment as a factor of species

# PCA on the fish abundance data
# ******************************

# Hellinger pre-transformation of the species data
Doubsspe.h <- decostand(DoubsSpe, "hellinger")
(spe.h.pca <- rda(Doubsspe.h))

# Plot eigenvalues and % of variance for each axis
ev <- spe.h.pca$CA$eig
dev.new(title="PCA eigenvalues")
evplot(ev)

# PCA biplots
spe.pca.sc1 <- scores(spe.h.pca, display="species", scaling=1)
spe.pca.sc2 <- scores(spe.h.pca, display="species", scaling=2)

dev.new(title="PCA on fish species", width=12, height=6)
par(mfrow=c(1,2))
cleanplot.pca(spe.h.pca, scaling=1, mar.percent=0.06)
cleanplot.pca(spe.h.pca, scaling=2, mar.percent=0.06)



# PCA on the environmental data using PCA.newr() and biplot.PCA.newr()
# ********************************************************************

# PCA; scaling 1 is the default for biplots
env.PCA.PL <- PCA.newr(env, stand=TRUE)
dev.new(title="PCA on environmental variables - scaling 1")
biplot.PCA.newr(env.PCA.PL)

# PCA; scaling 2 in the biplot
dev.new(title="PCA on environmental variables - scaling 2")
biplot.PCA.newr(env.PCA.PL, scaling=2)



# CA of the raw species dataset (original species abundances)
# ***********************************************************

# Compute CA
(spe.ca <- cca(spe))
summary(spe.ca)		# default scaling 2
summary(spe.ca, scaling=1)

# Plot eigenvalues and % of variance for each axis
(ev2 <- spe.ca$CA$eig)
dev.new(title="CA eigenvalues")
evplot(ev2)

# CA biplots
dev.new(title="CA biplots", width=14, height=7)
par(mfrow=c(1,2))
# Scaling 1: sites are centroids of species
plot(spe.ca, scaling=1, main="CA fish abundances - biplot scaling 1")
# Scaling 2 (default): species are centroids of sites
plot(spe.ca, main="CA fish abundances - biplot scaling 2")

# A posteriori projection of environmental variables in a CA
# The last plot produced (CA scaling 2) must be active
(spe.ca.env <- envfit(spe.ca, env, scaling=2)) # Scaling 2 is default
plot(spe.ca.env)
# Plot significant variables with a different colour
plot(spe.ca.env, p.max=0.05, col=3)

# Species data table ordered after the CA result
vegemite(spe, spe.ca)

##################################################
ENV data
spe <- spe[-8,]
env <- env[-8,]
spa <- spa[-8,]



# Compute CA
(env.ca <- cca(env))
summary(env.ca)		# default scaling 2
summary(env.ca, scaling=1)

# Plot eigenvalues and % of variance for each axis
(ev2 <- env.ca$CA$eig)
# dev.new(title="CA eigenvalues")
evplot(ev2)

# CA biplots
# dev.new(title="CA biplots", width=14, height=7)
par(mfrow=c(1,2))
# Scaling 1: sites are centroids of species
plot(env.ca, scaling=1, main="CA env variables - biplot scaling 1")
# Scaling 2 (default): species are centroids of sites
plot(env.ca, main="CA env variables - biplot scaling 2")

# A surface plot for the CA
require('viridis')
palette(viridis(8))
par(mar = c(4, 4, 0.9, 0.5) + .1, mfrow = c(2, 2))
with(env, tmp <- ordisurf(env.ca ~ alt, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Alt"))
abline(h = 0, v = 0, lty = 3)
with(env, tmp <- ordisurf(env.ca ~ flo, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Flo"))
abline(h = 0, v = 0, lty = 3)
with(env, tmp <- ordisurf(env.ca ~ bod, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Bod"))
abline(h = 0, v = 0, lty = 3)
with(env, tmp <- ordisurf(env.ca ~ amm, bubble = 3,
                          family = quasipoisson, knots = 2, col = 6,
                          display = "sites", main = "Amm"))
abline(h = 0, v = 0, lty = 3)

# A posteriori projection of environmental variables in a CA
# The last plot produced (CA scaling 2) must be active
(env.ca.env <- envfit(env.ca, env, scaling = 2)) # Scaling 2 is default
plot(env.ca.env)
# Plot significant variables with a different colour
plot(env.ca.env, p.max = 0.05, col = "red")

follows no linear distribution

        
        
        
        

# Species data table ordered after the CA result
vegemite(env, env.ca)







# CA using CA.newr() function
# ***************************

spe.CA.PL <- CA.newr(spe)
dev.new(title="Species CA with CA() function", width=14, height=8)
biplot.CA(spe.CA.PL, cex=1)

# Ordering of the data table following the first CA axis
# The table is transposed, as in vegemite() output
summary(spe.CA.PL)
t(spe[order(spe.CA.PL$F[,1]), order(spe.CA.PL$V[,1])])


###########################################################################################################33
# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix of 
# fish species
# *********************************************************************

spe.bray <- vegdist(spe)
spe.b.pcoa <- cmdscale(spe.bray, k=(nrow(spe)-1), eig=TRUE)
# Plot of the sites
dev.new(title="PCoA on fish species - Percentage difference")
ordiplot(scores(spe.b.pcoa, choices=c(1,2)), type="t", main="PCoA with species weighted averages")
abline(h=0, lty=3)
abline(v=0, lty=3)
# Add weighted average projection of species
spe.wa <- wascores(spe.b.pcoa$points[,1:2], spe)
text(spe.wa, rownames(spe.wa), cex=0.7, col="red")


# PCoA and projection of species vectors using function pcoa()
spe.h.pcoa <- pcoa(dist(spe.h))
# Biplots
dev.new(title="PCoA with species vectors", width=14, height=8)
par(mfrow=c(1,2))
# First biplot: Hellinger-transformed species data
biplot.pcoa(spe.h.pcoa, spe.h, dir.axis1=-1) 
abline(h=0, lty=3)
abline(v=0, lty=3)
# Second biplot: standardized Hellinger-transformed species data
spe.std <- scale(spe.h)
biplot.pcoa(spe.h.pcoa, spe.std, dir.axis1=-1) 
abline(h=0, lty=3)
abline(v=0, lty=3)


# Comparison of PCoA results with Euclidean and non-Euclidean
# dissimilarity matrices
# ***********************************************************

# PCoA on a Hellinger distance matrix
is.euclid(dist(spe.h))
summary(spe.h.pcoa) 
spe.h.pcoa$values

# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix
is.euclid(spe.bray)
spe.bray.pcoa <- pcoa(spe.bray) 
spe.bray.pcoa$values		# Observe eigenvalues 18 and following

# PCoA on the square root of a percentage difference (Bray-Curtis)
# dissimilarity matrix
is.euclid(sqrt(spe.bray))
spe.braysq.pcoa <- pcoa(sqrt(spe.bray))
spe.braysq.pcoa$values	# Observe the eigenvalues

# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix with
# Lingoes correction
spe.brayl.pcoa <- pcoa(spe.bray, correction="lingoes")
spe.brayl.pcoa$values		# Observe the eigenvalues, col. 1 and 2

# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix with
# Cailliez correction
spe.brayc.pcoa <- pcoa(spe.bray, correction="cailliez")
spe.brayc.pcoa$values		# Observe the eigenvalues, col. 1 and 2



# NMDS applied to the fish species - percentage difference
# dissimilarity matrix
# ********************************************************

spe.nmds <- metaMDS(spe, distance="bray")
spe.nmds
spe.nmds$stress
dev.new(title="NMDS on fish species - Percentage difference")
plot(spe.nmds, type="t", main=paste("NMDS/Percentage difference - Stress =", round(spe.nmds$stress,3)))

# Shepard plot and goodness of fit
dev.new(title="NMDS - Shepard plot", width=12, height=6)
par(mfrow=c(1,2))
stressplot(spe.nmds, main="Shepard plot")
gof <- goodness(spe.nmds)
plot(spe.nmds, type="t", main="Goodness of fit")
points(spe.nmds, display="sites", cex=gof*300)


# Add colours from a clustering result to an NMDS plot
# ****************************************************

# Ward clustering of percentage difference dissimilarity matrix
# and extraction of four groups
spe.bray.ward <- hclust(spe.bray, "ward.D") # In this case, better than ward.D2 for 4 groups
spe.bw.groups <- cutree(spe.bray.ward, k=4)
grp.lev <- levels(factor(spe.bw.groups))

# Combination with NMDS result
sit.sc <- scores(spe.nmds)
dev.new(title="NMDS plot with cluster colors")
p <- ordiplot(sit.sc, type="n", main="NMDS/Percentage difference + clusters Ward/Percentage difference")
for (i in 1:length(grp.lev))
{
	points(sit.sc[spe.bw.groups==i,], pch=(14+i), cex=2, col=i+1)
}
text(sit.sc, row.names(spe), pos=4, cex=0.7)
# Add the dendrogram
ordicluster(p, spe.bray.ward, col="dark grey")
# Add a legend interactively
legend(locator(1), paste("Group",c(1:length(grp.lev))), 
	pch=14+c(1:length(grp.lev)), col=1+c(1:length(grp.lev)), pt.cex=2)



# ------------------------------------------------------------------------------

# A simple function to perform PCA

myPCA <- function(Y)
{

	Y.mat <- as.matrix(Y)
	object.names <- rownames(Y)
	var.names <- colnames(Y)
	
	# Centre the data (will be needed to compute F)
	Y.cent <- scale(Y.mat, center=TRUE, scale=FALSE)
	
	# Covariance matrix S
	Y.cov <- cov(Y.cent)

	# Eigenvectors and eigenvalues of S (Legendre and Legendre 2012, 
	# eq. 9.1 and 9.2)
	Y.eig <- eigen(Y.cov)
	
	# Copy the eigenvectors to matrix U (used to represent variables
	# in scaling 1 biplots)
	U <- Y.eig$vectors
	rownames(U) <- var.names
	
	# Compute matrix F (used to represent objects in scaling 1 plots)
	F <- Y.cent%*%U			# eq. 9.4
	rownames(F) <- object.names

	# Compute matrix U2 (to represent variables in scaling 2 plots)
	# eq. 9.8
	U2 <- U%*%diag(Y.eig$values^0.5)
	rownames(U2) <- var.names

	# Compute matrix G (to represent objects in scaling 2 plots)
	# eq. 9.14
	G <- F%*%diag(Y.eig$values^0.5)
	rownames(G) <- object.names

	# Output of a list containing all the results
	result <- list(Y.eig$values,U,F,U2,G)
	names(result) <- c("eigenvalues","U", "F", "U2", "G")
	result
}

# ------------------------------------------------------------------------------

# PCA on fish species using hand-written function
fish.PCA <- myPCA(spe.h)
summary(fish.PCA)
# Eigenvalues
fish.PCA$eigenvalues
# Eigenvalues expressed as percentages
(pv <- round(100*fish.PCA$eigenvalues/sum(fish.PCA$eigenvalues),2))
# Alternate computation of total variation (denominator)
round(100*fish.PCA$eigenvalues/sum(diag(cov(spe.h))),2)
# Cumulative eigenvalues expressed as percentages
round(cumsum(100*fish.PCA$eigenvalues/sum(fish.PCA$eigenvalues)),2)

# Biplots
dev.new(title="PCA using homemade function", width=12, height=8)
par(mfrow=c(1,2))
# Scaling 1 biplot
biplot(fish.PCA$F, fish.PCA$U)
# Scaling 2 biplot
biplot(fish.PCA$G, fish.PCA$U2)

# Plots using generic R plot() function
dev.new(title="PCA using homemade function - generic plot function", width=12, height=8)
par(mfrow=c(1,2))
# Scaling 1
# Plot objects
plot(fish.PCA$F[,1], fish.PCA$F[,2], asp=1, main="PCA scaling 1", 
	xlab=paste("Axis 1 (", pv[1], "%)", sep=""), 
	ylab=paste("Axis 2 (", pv[2], "%)", sep=""))
# Plot variables
arrows(x0=0, y0=0, fish.PCA$U[,1], fish.PCA$U[,2], length=0.1, col="red")
# Add object numbers
text(fish.PCA$F[,1], fish.PCA$F[,2], labels=row.names(spe), pos=3, cex=0.8)
# Add variable names
text(fish.PCA$U[,1], fish.PCA$U[,2], labels=colnames(spe), adj=c(-0.2,0.2), 
	col="red", cex=0.8)
abline(h=0, lty=3)
abline(v=0, lty=3)
# Scaling 2
plot(fish.PCA$G[,1], fish.PCA$G[,2], asp=1, main="PCA scaling 2", 
	xlab=paste("Axis 1 (", pv[1], "%)", sep=""), 
	ylab=paste("Axis 2 (", pv[2], "%)", sep=""))
arrows(x0=0, y0=0, fish.PCA$U2[,1], fish.PCA$U2[,2], length=0.1, col="red")
text(fish.PCA$G[,1], fish.PCA$G[,2], labels=row.names(spe), pos=3, cex=0.8)
text(fish.PCA$U2[,1], fish.PCA$U2[,2], labels=colnames(spe), col="red", 
	adj=c(-0.2,0.2), cex=0.8)
abline(h=0, lty=3)
abline(v=0, lty=3)




str(dune)
data("dune")
data("dune.env")
dune.env
