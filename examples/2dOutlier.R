source('../rfr_us.R')
library(ggplot2)

genCube <- function(nrPoints,nrDim,center=rep(0,nrDim), sideLength=1.0) {
	dat <- matrix(nrow=nrPoints,ncol=nrDim)
	r <- sideLength/2.0
	for(i in 1:nrPoints){
		dat[i,] <- runif(nrDim,-r-center,r+center)
	}
	dat
}

boundingSideLength <- 2

pointsInCube <- genCube(500,2)
sM <- urerf(pointsInCube, 100, 3)

labelsTraining <- rep(500, TRUE)

training <- data.frame(Xpoints = pointsInCube[,1], Ypoints = pointsInCube[,2], labels=labelsTraining)

png(file="~/dropbox/results/outlierTraining.png")
p <- ggplot(training, aes(x=Xpoints, y=Ypoints, colour=labels)) + geom_point()
p <- p + xlim(-1,1) + ylim(-1,1) + guides(fill=FALSE)
p <- p + labs(title="Training for 1x1 Square")

print(p)
dev.off()


testPoints <- genCube(500,2, sideLength=boundingSideLength)

testLabelsActual <- sapply(1:nrow(testPoints), function(x){
		 all(testPoints[x,] < .5) & all(testPoints[x,] > -.5)
})


testLabelsCalculated <- is.outlier(testPoints, sM, 2.0)

results <- data.frame(Xpoints = testPoints[,1], Ypoints = testPoints[,2], labels=testLabelsCalculated)

png(file="~/dropbox/results/outlierResults.png")
p <- ggplot(results, aes(x=Xpoints, y=Ypoints, colour=labels)) + geom_point()
p <- p + labs(title="Is it in the Square?\nOutlier Detection For 1x1 Square")

print(p)
dev.off()

