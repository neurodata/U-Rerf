source('../rfr_us.R')
library(ggplot2)

# number of trees for forest
numtrees <- 100
# number of dimensions in dataset
m <- 10
# number of samples in dataset
sizeD <- c(10,100,1000,10000, 20000)
maxDepth <- 8

testNum <- 1
testSize <- NA
testTime <- NA
testType <- factor(NA, levels=c("depth", "minK"))



for(numberOfSamples in sizeD){
# create a sizeD by m synthetic dataset
X <- matrix(runif(m*numberOfSamples), nrow=numberOfSamples, ncol=m)
# create a similarity matrix using urerf
start.time <- Sys.time()
similarityMatrix <- createSimilarityMatrixDepth (X, numtrees, maxDepth)
end.time <- Sys.time()
time.taken <- end.time-start.time
testSize[testNum] <- numberOfSamples
testTime[testNum] <- time.taken
testType[testNum] <- "depth"
testNum <- testNum + 1
print(time.taken)
}

K <- 10

for(numberOfSamples in sizeD){
# create a sizeD by m synthetic dataset
X <- matrix(runif(m*numberOfSamples), nrow=numberOfSamples, ncol=m)
# create a similarity matrix using urerf
start.time <- Sys.time()
similarityMatrix <- createSimilarityMatrix (X, numtrees, K)
end.time <- Sys.time()
time.taken <- end.time-start.time
testSize[testNum] <- numberOfSamples
testTime[testNum] <- time.taken
testType[testNum] <- "minK"
testNum <- testNum + 1

print(time.taken)
}

testSize <- as.factor(testSize)
df <- data.frame(dataSize=testSize,time=testTime,type=testType)
png(file="../results/timetest.png")
ggplot(aes(x = dataSize, y = time, colour=type), data = df) + geom_point()
dev.off()

