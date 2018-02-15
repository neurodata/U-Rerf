source('../rfr_us.R')

# number of trees for forest
numtrees <- 100
# the 'k' of k nearest neighbors
k <- 10


# create a sizeD by m synthetic dataset
X <- as.matrix(iris[,1:4])

# create a similarity matrix using urerf
similarityMatrix <- createSimilarityMatrix (X, numtrees, k)
distanceMatrix <- createDistanceMatrix (X, numtrees, k)


