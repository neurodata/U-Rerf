
X <- as.matrix(iris[,1:4])
X <- cbind(X, runif(150, 0, 7.9))
X <- cbind(X, runif(150, 0, 7.9))

aN <- hclust(dist(X), method="average")
cutree(aN, k=3)
summary(as.factor(cutree(aN,k=3)[1:50]))
summary(as.factor(cutree(aN,k=3)[51:100]))
summary(as.factor(cutree(aN,k=3)[101:150]))


aN <- hclust(dist(X)), method="mcquitty")
cutree(aN, k=3)
summary(as.factor(cutree(aN,k=3)[1:50]))
summary(as.factor(cutree(aN,k=3)[51:100]))
summary(as.factor(cutree(aN,k=3)[101:150]))



source('../rfr_us.R')
library(ggplot2)

# number of trees for forest
numtrees <- 100
# number of dimensions in dataset
K <- 2
# max depth
maxDepth <- 8


# create a sizeD by m synthetic dataset

similarityMatrix <- 1-createSimilarityMatrix (X, numtrees, K)

aM <- hclust(as.dist(similarityMatrix), method="average")
cutree(aM, k=3)
summary(as.factor(cutree(aM,k=3)[1:50]))
summary(as.factor(cutree(aM,k=3)[51:100]))
summary(as.factor(cutree(aM,k=3)[101:150]))




aM <- hclust(as.dist(similarityMatrix), method="mcquitty")
cutree(aM, k=3)
summary(as.factor(cutree(aM,k=3)[1:50]))
summary(as.factor(cutree(aM,k=3)[51:100]))
summary(as.factor(cutree(aM,k=3)[101:150]))




X <- read.csv("yeast.data", sep="", header=FALSE)
X <- with(X, X[order(V10),])
Y <- X[,ncol(X)]
X <- as.matrix(X[,2:(ncol(X)-1)])

similarityMatrix <- 1-createSimilarityMatrix (X, numtrees, K)

aM <- hclust(as.dist(similarityMatrix), method="mcquitty")
cutree(aM, k=10)








