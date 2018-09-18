source('../rfr_us.R')

genCube <- function(nrPoints,nrDim,center=rep(0,nrDim), sideLength=1.0) {
	dat <- matrix(nrow=nrPoints,ncol=nrDim)
	r <- sideLength/2.0
	for(i in 1:nrPoints){
		dat[i,] <- runif(nrDim,-r-center,r+center)
	}
	dat
}

boundingSideLength <- 2

pointsInCube <- genCube(500,3)
sM <- urerf(pointsInCube, 100, 3)

testPoints <- genCube(500,3, sideLength=boundingSideLength)
testLabelsActual <- sapply(1:nrow(testPoints), function(x){
		 all(testPoints[x,] < .5) & all(testPoints[x,] > -.5)
})

testLabelsCalculated <- !is.outlier(testPoints, sM, 2.5)

positiveOutcomes <- which(testLabelsActual==TRUE)

sensitivity <- sum(testLabelsActual[positiveOutcomes] == testLabelsCalculated[positiveOutcomes])/length(positiveOutcomes)


specificity <- sum(testLabelsActual[-positiveOutcomes] == testLabelsCalculated[-positiveOutcomes])/length(testLabelsActual[-positiveOutcomes])
print(sensitivity)
print(specificity)
