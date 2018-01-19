source('../rfr_us.R')
library(ggplot2)
library(scatterplot3d)

# number of trees for forest
numtrees <- 100
# number of dimensions in dataset
m <- 3
# number of samples in dataset
sizeD <- 1000
# the 'k' of k nearest neighbors
k = 10


# create a sizeD by m synthetic dataset
X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)
X <- (X + runif(m*sizeD)/10)
X[500,1] <-X[500,1] + .15

AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

# find the actual euclidean distance between all samples of the synthetic dataset
for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}

# create a similarity matrix using urerf
similarityMatrix <- createSimilarityMatrix (X, numtrees, k)
nnzPts <- which(similarityMatrix != 0)

#create output
pdf(file="../results/3dimDvNwO.pdf")
s3d <- scatterplot3d(X[c(1:499, 501:1000),2], X[c(1:499, 501:1000),3], X[c(1:499, 501:1000),1],        # x y and z axis
                 color="blue", pch=19, # filled blue circles
                 main="Noise 3-D Line w/ Outlier",
                 xlab="X",
                 ylab="Y",
                 zlab="Z")

s3d$points3d(X[500,2], X[500,3], X[500,1], color="green", pch=17)

ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = similarityMatrix[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point() + labs(title="Distance vs Nearness\nNoisy 3-D Line w/outlier, n=1000, d=3, k=10, trees=100\nThree Samples (0 Nearness Omitted)")

nnzPts <- which(t(similarityMatrix[499:501,]) != 0)
ssd <- data.frame(Distance = t(AkNN[499:501,])[nnzPts], Nearness = t(similarityMatrix[499:501,])[nnzPts])
groupLabels <- c(rep("499",sizeD), rep("500", sizeD), rep("501",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\nNoisy 3-D Line w/outlier, n=1000, d=3, k=10, trees=100\nThree Samples (0 Nearness Omitted)\n500 is outlier")+ geom_jitter()
dev.off()

