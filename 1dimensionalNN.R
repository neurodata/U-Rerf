source('rfr_us.R')
library(ggplot2)
library(reshape2)
library(scatterplot3d)



# number of trees for forest
numtrees <- 100
# number of dimensions in dataset
m <- 1
# number of samples in dataset
sizeD <- 1000
# the 'k' of k nearest neighbors
k = 10


# create a sizeD by m synthetic dataset
X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)
AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

# find the actual euclidean distance between all samples of the synthetic dataset
for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}

# create a similarity matrix using urerf
similarityMatrix <- createSimilarityMatrix (X, numtrees, k)
nnzPts <- which(similarityMatrix != 0)

#create output
pdf(file="1dimDvN.pdf")
ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = similarityMatrix[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point() + labs(title="Distance vs Nearness\n1-D Line, n=1000, d=1, k=10, trees=100\nAll Samples (0 Nearness Omitted)")

nnzPts <- which(t(similarityMatrix[499:501,]) != 0)
ssd <- data.frame(Distance = t(AkNN[499:501,])[nnzPts], Nearness = t(similarityMatrix[499:501,])[nnzPts])
groupLabels <- c(rep("499",sizeD), rep("500", sizeD), rep("501",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\n1-D Line, n=1000, d=1, k=10, trees=100\nThree Samples (0 Nearness Omitted)")
dev.off()

