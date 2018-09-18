source('../rfr_us.R')

genCube <- function(nrPoints,nrDim,center=rep(0,nrDim), sideLength=1.0) {
	dat <- matrix(nrow=nrPoints,ncol=nrDim)
	r <- sideLength/2.0
	for(i in 1:nrPoints){
		dat[i,] <- runif(nrDim,-r-center,r+center)
	}
	dat
}


	pointsInCube <- genCube(500,3)

	sM <- urerf(pointsInCube)

	#actual volume of sphere
	volumeOfCube <- 1

	#approximate volume of sphere
	mV <- manifoldVolume(sM, 500)

	print(paste("actual volume of Cube: ", volumeOfCube))
	print(paste("approximated volume of Cube: ", mV))

