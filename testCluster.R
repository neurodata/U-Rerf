source('rfr_us.R')
library(ggplot2)

# number of trees for forest
numtrees <- 100
# number of dimensions in dataset
m <- 1
# number of samples in dataset
sizeD <- 1000
# the 'k' of k nearest neighbors
k = 3


# create a sizeD by m synthetic dataset
X <- as.matrix(iris[,1:4])

# create a similarity matrix using urerf
similarityMatrix <- createSimilarityMatrix (X, numtrees, k)

#print(findClusters(similarityMatrix, 3,k)) 
pa <- findClusters(similarityMatrix, 3,k)
print(pa)

