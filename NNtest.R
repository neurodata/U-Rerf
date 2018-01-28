source("rfr_us.R")
k <- 10
X <- as.matrix(iris[,1:4])
testSize <- 10
numSamples <- nrow(X)
testIndex <- sample(1:numSamples, testSize, replace=F)
testSamples <- X[testIndex,]
X <- X[-testIndex,]

AkNN <- matrix(0, nrow=testSize, ncol=(numSamples-testSize))
AppkNN <- matrix(0, nrow=testSize, ncol=(numSamples-testSize))

for(z in 1:testSize){
    nearK <- order(sqrt(rowSums(sweep(X,2,testSamples[z,         ])^2)))[1:k]
    AkNN[z,nearK] <- 1
}

forest <- invisible(rfrus(X,trees=200,MinParent=k))
    NN1 <- distNN(testSamples, forest, numSamples-testSize)

sensitivity <- rep(0,41) 
specificity <- NA
for(m in 0:40){
    numK <- testSize + m
for(j in 1:testSize){
    sensitivity[m+1] <- sensitivity[m+1] + sum(order(AkNN[j,], decreasing=T)[1:k] %in% order(NN1[j,], decreasing = T)[1:numK])
}
sensitivity[m+1] <- sensitivity[m+1]/(testSize*k)
specificity[m+1] <- m/(numSamples-testSize-k)
}
    
