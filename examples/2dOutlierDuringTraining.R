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

numInstances <- 500
numDims <- 2
boundingSideLength <- 2

pointsInCube <- genCube(numInstances,numDims)
pointsInCube <- rbind(pointsInCube, c(.75, .6), c(-.6,.3))
sM <- urerf(pointsInCube, 100, 3)

labelsTraining <- c(rep(TRUE,numInstances),FALSE,FALSE)

training <- data.frame(Xpoints = pointsInCube[,1], Ypoints = pointsInCube[,2], labels=labelsTraining)

png(file="~/dropbox/results/2dOutlierDuringTraining1.png")
p <- ggplot(training, aes(x=Xpoints, y=Ypoints, colour=labels)) + geom_point()
p <- p + xlim(-1,1) + ylim(-1,1) + guides(fill=FALSE)
p <- p + labs(title="Training for 1x1 Square")

print(p)
dev.off()

outliersFromData <- dataset.outlier(sM, 2)


labelsTraining <- c(rep(TRUE,numInstances+2))
labelsTraining[outliersFromData] <- FALSE

training <- data.frame(Xpoints = pointsInCube[,1], Ypoints = pointsInCube[,2], labels=labelsTraining)

png(file="~/dropbox/results/2dOutlierDuringTraining2.png")
p <- ggplot(training, aes(x=Xpoints, y=Ypoints, colour=labels)) + geom_point()
p <- p + xlim(-1,1) + ylim(-1,1) + guides(fill=FALSE)
p <- p + labs(title="Training for 1x1 Square")

print(p)
dev.off()

