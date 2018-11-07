source('../rfr_us.R')
library(ggplot2)
library(scatterplot3d)

# number of trees for forest
numtrees <- 100
# number of samples in dataset
sizeD <- 1000
# the 'k' of k nearest neighbors
k = 33
set.seed(22)


# create a sizeD by m synthetic dataset
X <- swissRoll(sizeD/2, size =1, dim3=T)
X <- as.matrix(X)[,2:4]
X[500,2] <- -1.0
X[500,3] <- -1.0
X[500,1] <- 1.0

#AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

# find the actual euclidean distance between all samples of the synthetic dataset
#for(z in 1:sizeD){
#AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
#}

# create a similarity matrix using urerf
#nnzPts <- which(similarityMatrix != 0)

#create output
png(file="~/dropbox/results/urerf/3dswissDvNwOutlier.png")
#png(file="../results/3dswissDvNwOutlier.png")
s3d <- scatterplot3d(X[c(1:499, 501:1000),2], X[c(1:499, 501:1000),3], X[c(1:499, 501:1000),1],        # x y and z axis
										                  color="blue", pch=19, # filled blue circles
																			                 main="",
																			                 #main="3-D Swiss Roll w/ Outlier",
																			                 xlab="X",
																											                  ylab="Y",
																											                  zlab="Z")

s3d$points3d(X[500,2], X[500,3], X[500,1], col="red", pch=17, cex=1)
#s3d$points3d(X[579,2], X[579,3], X[579,1], color="green", pch=17, cex=2)
dev.off()

####################################
######## Random Forest #############
library(randomForest)
proxMat <- randomForest::randomForest(X)
outliers <- randomForest::outlier(proxMat)

png(file="~/dropbox/results/urerf/rfTest.png")
plot(outliers, pch="*", cex=ifelse(names(outliers)=="500", 4, 2), main="Random Forest Outlier Scores", col=ifelse(names(outliers)=="500", "red", "black"))  #plot cooks distance;
abline(h = 3, col="red")  # add cutoff line
#text(x=1:length(outliers)+1, y=outliers, labels=ifelse(outliers>3,names(outliers),""), col="red") 
dev.off()


####################################
######## Outliers ##################
library(outliers)
oLs <- rowSums(outliers::outlier(X, logical=TRUE))
attr(oLs, "names") <- 1:length(oLs)

png(file="~/dropbox/results/urerf/outliersTest.png")
plot(oLs, pch="*", cex=ifelse(names(outliers)=="500", 4, 2), main="'Outliers' Outlier Scores", col=ifelse(names(outliers)=="500", "red", "black"))  #plot cook's distance;
#text(x=1:length(oLs)+1, y=oLs, labels=ifelse(oLs>0,names(oLs),""), col="red") 
dev.off()


####################################
######## U-RerF ####################
source('../rfr_us.R')
similarityMatrix <- urerf(X, numtrees, k)

outliersUrerf <- dataset.outlier(similarityMatrix)
which(order(outliersUrerf)==500)

png(file="~/dropbox/results/urerf/urerfTest.png")
plot(outliersUrerf, cex=ifelse(names(outliersUrerf)=="500", 1, 2), main="",xaxt='n', ylab="Outlier Score",xlab="", pch= ifelse(names(outliers)=="500", 17, 19), col=ifelse(names(outliers)=="500", "red", "blue"));
#abline(h = 3, col="red")  # add cutoff line
#text(x=1:length(outliers)+1, y=outliers, labels=ifelse(outliers>3,names(outliers),""), col="red") 
dev.off()
