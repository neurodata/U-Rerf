source('../rfr_us.R')

genSphereRejection <- function(nrPoints,nrDim,center=rep(0,nrDim),r=1.0) {
	dat <- matrix(nrow=nrPoints,ncol=nrDim)
	nAccept <- 0
	while (nAccept < nrPoints) {
		samp <- runif(nrDim,-r-center,r+center)
		if ( sqrt(sum( (samp-center)^2 )) < r ) {
			nAccept = nAccept + 1
			dat[nAccept,] <- samp
		}
	}
	dat
}


	pointsInSphere <- genSphereRejection(500,3)

	sM <- urerf(pointsInSphere)

	#actual volume of sphere
	volumeOfSphere <- 4.189

	#approximate volume of sphere
	mV <- manifoldVolume(sM, 25000)

	print(paste("actual volume of sphere: ", volumeOfSphere))
	print(paste("approximated volume of sphere: ", mV))
