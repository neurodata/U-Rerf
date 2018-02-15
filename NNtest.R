#library(ggplot)
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
numTrees <- c(25, 50, 100, 200, 400)
#numTrees <- c(25, 50)
numLen <- 40
results <- data.frame(sensitivity=double(), specificity=double(), lineNum=integer())
currentResult <- 1

for(numT in numTrees){
    forest <- invisible(rfrus(X,trees=numT,MinParent=k))
    NN1 <- distNN(testSamples, forest, numSamples-testSize)

    specificity <- NA
    for(m in 1:numLen){
        tempSensitivity <- 0
        numK <- testSize + m
        for(j in 1:testSize){
            tempSensitivity <- tempSensitivity + sum(order(AkNN[j,], decreasing=T)[1:k] %in% order(NN1[j,], decreasing = T)[1:numK])
        }
        results <- rbind(results, c(tempSensitivity/(testSize*k), 1-m/(numSamples-testSize-k), numT))
    }
}

colnames(results) <- c("specificity", "sensitivity", "trees")


png(filename="spec.png", height=500, width=500)
xRange <- rev(range(results[,2]))
yRange <- range(results[,1])
plot(xRange, yRange, type="n", xlab="specificity",
        ylab="sensitivity" ) 
colors <- rainbow(length(numTrees))
linetype <- c(1:length(numTrees)) 
plotchar <- seq(18,18+length(numTrees),1)

for (i in 1:length(numTrees)) { 
      res <- subset(results, trees==numTrees[i]) 
  #lines(res$specificity, res$sensitivity, type="b", lwd=1.5, lty=linetype[i], col=colors[i], pch=plotchar[i]) 
  lines(res$specificity, res$sensitivity, col=colors[i]) 
} 

title("Sensitivity vs Specificity of Iris Nearest Neighbors Using URerf")
legend(xRange[2], (yRange[1]+yRange[2])/2, numTrees, cex=0.8, col=colors, pch=plotchar, lty=linetype, title="Number of Trees")
dev.off()
    
