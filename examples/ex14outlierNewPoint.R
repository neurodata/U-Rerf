source('../rfr_us.R')
library(ggplot2)
library(scatterplot3d)
library(MASS)


set.seed(10)
# number of trees for forest
numtrees <- 200
# number of samples in dataset
sizeD <- 2000
# the 'k' of k nearest neighbors
#k = 100 
#k = sqrt(sizeD) 
k =  5

X1 <- runif(sizeD, 0,1)
X2 <- runif(sizeD, 0,1)
X <- as.matrix(cbind(X1, X2))


tModel <- urerf(X, numtrees,K=k )

k = 5
distanceOldOnManifold <- newPointDist(X, tModel, k)


X1 <- runif(sizeD, 0,1)
X2 <- runif(sizeD, 0,1)
X <- as.matrix(cbind(X1, X2))

distanceNewOnManifold <- newPointDist(X, tModel, k)



X1a <- runif(sizeD/2, -2, 2)
X2b <- runif(sizeD/2, -2, 2)
X1b <- runif(sizeD/2, 0,1)
X2a <- runif(sizeD/2, 0,1)

X1 <- c(X1a, X1b)
X2 <- c(X2a, X2b)

#X1[which(X1 >= 0 & X1 <= 1)] <- X1[which(X1 >= 0 & X1 <= 1)]-1
#X2[which(X2 >= 0 & X2 <= 1)] <- X2[which(X2 >= 0 & X2 <= 1)]-1
X <- as.matrix(cbind(X1, X2))

distanceNewOffManifold <- newPointDist(X, tModel, k)

summary(distanceOldOnManifold)
summary(distanceNewOnManifold)
summary(distanceNewOffManifold)
