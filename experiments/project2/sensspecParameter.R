source('../../rfr_us.R')
library(ggplot2)
library(FNN)

# number of trees for forest
numtrees <- 100
# number of dimensions in dataset
m <- 10 
# number of samples in dataset
sizeD <- 4000
# the 'k' of k nearest neighbors
k = 10

depth = 3
topK <- 3

X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)
kNN <- as.matrix(dist(X))
Y <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)






X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)
kNN <- as.matrix(dist(X))


# create a similarity matrix using urerf
sM <- urerfDepth(X, numtrees, 8)
sM2 <- urerf(X, numtrees, k)

Y <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)

png(file="~/dropbox/results/simpleLines.png")
plot(Y)
points(Y, cex = 1.0, col = "green")
plot(X)
points(X, cex = .5, col = "dark red")
dev.off()

topK <- 3

NN <- (get.knnx(X,Y,k=topK,algorithm="brute"))$nn.index
apprNN <- ann(Y, sM2, 200)
apprNND <- ann(Y, sM, 200)

results <- data.frame(sensitivity=double(), specificity=double(), manifolds=integer())

	q<-1
for(i in topK:200){
	numRight <- 0
	numNegs <- 0
for(j in 1:nrow(Y)){
	right <- sum(NN[j,] %in% apprNND[j,1:i])
numRight <- numRight + right
numNegs <- numNegs + nrow(X)-i-(topK-right)
}
results[q,] <- c(numRight/(nrow(Y)*topK), numNegs/((nrow(X)-topK)*nrow(Y)),1)
q <- q+1
}

png(file="~/dropbox/results/sensSpec1Manifold.png")
p <- ggplot(aes(x = specificity, y = sensitivity), data = results) + geom_line() + labs(title="Sensitivity vs Specificity\n10-D Line, n=1000, d=10, k=10, trees=100")
p <- p + scale_x_reverse()
print(p)
dev.off()


sizeD <- 2 * sizeD

X <- matrix(sort(runif(m*sizeD/2)), nrow=sizeD/2, ncol=m)
X <- rbind(matrix(sort(runif(m*sizeD/2)+1), nrow=sizeD/2, ncol=m))
kNN <- as.matrix(dist(X))

# create a similarity matrix using urerf
sM <- urerfDepth(X, numtrees, 8)
sM2 <- urerf(X, numtrees, k)

Y <- matrix(sort(runif(m*sizeD/2)), nrow=sizeD/2, ncol=m)
Y <- rbind(matrix(sort(runif(m*sizeD/2)+1), nrow=sizeD/2, ncol=m))

NN <- (get.knnx(X,Y,k=topK,algorithm="brute"))$nn.index
apprNN <- ann(Y, sM2, 200)
apprNND <- ann(Y, sM, 200)

#results <- data.frame(sensitivity=double(), specificity=double())

#	q<-1
for(i in topK:200){
	numRight <- 0
	numNegs <- 0
for(j in 1:nrow(Y)){
	right <- sum(NN[j,] %in% apprNND[j,1:i])
numRight <- numRight + right
numNegs <- numNegs + nrow(X)-i-(topK-right)
}
results[q,] <- c(numRight/(nrow(Y)*topK), numNegs/((nrow(X)-topK)*nrow(Y)),2)
q <- q+1
}

png(file="~/dropbox/results/sensSpec2Manifold.png")
p <- ggplot(aes(x = specificity, y = sensitivity, colour=manifolds), data = results) + geom_line() + labs(title="Sensitivity vs Specificity\n10-D Line, n=1000, d=10, k=10, trees=100, 2 manifolds")
p <- p + scale_x_reverse()
print(p)
dev.off()


