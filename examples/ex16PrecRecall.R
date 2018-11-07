source('../rfr_us.R')
library(ggplot2)
library(FNN)

LoadAndPrintOnly <- FALSE
args <- commandArgs(trailingOnly = TRUE)
if(length(args)==1){
	if(args[1] == "print"){
		LoadAndPrintOnly <- TRUE
	}
}

if(!LoadAndPrintOnly){
# number of trees for forest
numtrees <- 100
# number of dimensions in dataset
m <- 10 
# number of samples in dataset
sizes <- c(100,1000)
#sizes <- c(100,1000,5000,10000)
# the 'k' of k nearest neighbors
#k = round(sizes^.5)

topK <- 3

results <- data.frame(recall=double(), precision=double(), datasetSize=integer())
	q<-1

for(sizeD in sizes){
	k <- round(sizeD^.5)
X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)
kNN <- as.matrix(dist(X))

# create a similarity matrix using urerf
#sM <- urerf(X, numtrees, depth=8)
sM <- urerf(X, numtrees, K=k)

Y <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)

	numToTest <- .2 * sizeD

NN <- (get.knnx(X,Y,k=topK,algorithm="brute"))$nn.index
apprNND <- ann(Y, sM, numToTest)

for(i in topK:numToTest){
	numRight <- 0
	numNegs <- 0
for(j in 1:nrow(Y)){
	right <- sum(NN[j,] %in% apprNND[j,1:i])
numRight <- numRight + right
#numNegs <- numNegs + nrow(X)-i-(topK-right)
}
results[q,] <- c(numRight/(nrow(Y)*topK), numRight/i,sizeD)
q <- q+1
}

}

results[,3] <- as.factor(results[,3])
save(results, file="resultsData/ex16_sensSpecSizes.Rdata")
}else{
load(file="resultsData/ex16_sensSpecSizes.Rdata")
}

png(file="results/ex16_sensSpecSizes.png")
p <- ggplot(aes(x = precision, y = recall, colour=datasetSize), data = results) + geom_line() + labs(title="Sensitivity vs Specificity\n10-D Line, n=100-10000, depth=8, trees=100")
p <- p + scale_x_reverse()
print(p)
dev.off()
