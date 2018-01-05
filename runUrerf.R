#source('rfr_us_p.R')
source('rfr_us.R')
library(ggplot2)
library(reshape2)
library(scatterplot3d)
require(compiler)

numtrees <- 100
m <- 1
sizeD <- 10
k = 3

#prepare data
rfrus <- cmpfun(rfrus)
distNN <- cmpfun(distNN)





numtrees <- 100
m <- 1
sizeD <- 1000
k = 10

AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)
AppkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

## Remove samples at random for testing.
X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m)

for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}

#NNorder <- order(NN1, decreasing=TRUE)
forest <- invisible(rfrus(X,trees=numtrees, MinParent = k))
for(z in 1:sizeD){
  NN1 <- distNN(X[z,], X, forest)
  AppkNN[z,] = NN1
for(q in 1:sizeD){
if(NN1[q]==0){
  AppkNN[z,q]<-0
}
}
  }

nnzPts <- which(AppkNN != 0)

ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = AppkNN[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point() + labs(title="Distance vs Nearness\n1-D Line, n=1000, d=1, k=10, trees=100\nAll Samples (0 Nearness Omitted)")

nnzPts <- which(t(AppkNN[499:501,]) != 0)
ssd <- data.frame(Distance = t(AkNN[499:501,])[nnzPts], Nearness = t(AppkNN[499:501,])[nnzPts])
groupLabels <- c(rep("499",sizeD), rep("500", sizeD), rep("501",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\n1-D Line, n=1000, d=1, k=10, trees=100\nThree Samples (0 Nearness Omitted)")

rm(forest)









# 3 dimensional line

```{r 3D, cache = FALSE, echo=FALSE, message=FALSE}
numtrees <- 100
m <- 3
sizeD <- 1000
k = 10

AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)
AppkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

## Remove samples at random for testing.
X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m, byrow=TRUE)

#NNorder <- order(NN1, decreasing=TRUE)
#print(X)
# for(z in 1:sizeD){
# nearK <- order(sqrt(rowSums(sweep(X,2,X[z,])^2)))[1:k]
# AkNN[z,nearK] <- 1
# }

#NNorder <- order(NN1, decreasing=TRUE)
# forest <- rfrus(X,trees=numtrees, MinParent = k, bagging=0)
# for(z in 1:sizeD){
#   NN1 <- distNN(X[z,], X, forest)
# nearK <- order(NN1,decreasing = TRUE)[1:k]
# AppkNN[z,nearK] <- 1
# }

#NNorder <- order(NN1, decreasing=TRUE)
for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}

#NNorder <- order(NN1, decreasing=TRUE)
forest <- invisible(rfrus(X,trees=numtrees, MinParent = k))
for(z in 1:sizeD){
  NN1 <- distNN(X[z,], X, forest)
  AppkNN[z,] = NN1
}

# #NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN != 0)
# plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=3", font.main = 4)
# 
# #NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN[,500] != 0)
# plot(cbind(AppkNN[nnzPts,500],AkNN[nnzPts,500]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (1 sample, #500), n=1000, d=3", font.main = 4)

nnzPts <- which(AppkNN != 0)
ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = AppkNN[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point() + labs(title="Distance vs Nearness\n3-D Line, n=1000, d=3, k=10, trees=100\nThree Samples (0 Nearness Omitted)")

#plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
#title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=1", font.main = 4)

#NNorder <- order(NN1, decreasing=TRUE)
nnzPts <- which(t(AppkNN[499:501,]) != 0)
ssd <- data.frame(Distance = t(AkNN[499:501,])[nnzPts], Nearness = t(AppkNN[499:501,])[nnzPts])
groupLabels <- c(rep("499",sizeD), rep("500", sizeD), rep("501",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\n3-D Line, n=1000, d=3, k=10, trees=100\nThree Samples (0 Nearness Omitted)")+ geom_jitter()



rm(forest)
```


# 10 dimensional line

```{r 10D, cache = FALSE, echo=FALSE, message=FALSE}
numtrees <- 100
m <- 10
sizeD <- 1000
k = 10

AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)
AppkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m, byrow=TRUE)

#NNorder <- order(NN1, decreasing=TRUE)
# #print(X)
# for(z in 1:sizeD){
# nearK <- order(sqrt(rowSums(sweep(X,2,X[z,])^2)))[1:k]
# AkNN[z,nearK] <- 1
# }

# #NNorder <- order(NN1, decreasing=TRUE)
# forest <- rfrus(X,trees=numtrees, MinParent = k, bagging=0)
# for(z in 1:sizeD){
#   NN1 <- distNN(X[z,], X, forest)
# nearK <- order(NN1,decreasing = TRUE)[1:k]
# AppkNN[z,nearK] <- 1
# }

#NNorder <- order(NN1, decreasing=TRUE)
for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}


#NNorder <- order(NN1, decreasing=TRUE)
forest <- rfrus(X,trees=numtrees, MinParent = k)
for(z in 1:sizeD){
  NN1 <- distNN(X[z,], X, forest)
  AppkNN[z,] = NN1
}

# #NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN != 0)
# plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=10", font.main = 4)
# 
# #NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN[,500] != 0)
# plot(cbind(AppkNN[nnzPts,500],AkNN[nnzPts,500]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (1 sample, #500), n=1000, d=10", font.main = 4)

nnzPts <- which(AppkNN != 0)
ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = AppkNN[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\n10-D Line, n=1000, d=10, k=10, trees=100\nAll Samples (0 Nearness Omitted)")

#plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
#title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=1", font.main = 4)

#NNorder <- order(NN1, decreasing=TRUE)
nnzPts <- which(t(AppkNN[499:501,]) != 0)
ssd <- data.frame(Distance = t(AkNN[499:501,])[nnzPts], Nearness = t(AppkNN[499:501,])[nnzPts])
groupLabels <- c(rep("499",sizeD), rep("500", sizeD), rep("501",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\n10-D Line, n=1000, d=10, k=10, trees=100\nThree Samples (0 Nearness Omitted)")+ geom_jitter()




rm(forest)
```


# noisy 3 dimensional line with outlier

```{r 3DwOutlier, cache = FALSE, echo=FALSE, message=FALSE}
numtrees <- 100
m <- 3
sizeD <- 1000
k = 10

AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)
AppkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

## Remove samples at random for testing.
X <- matrix(sort(runif(m*sizeD)), nrow=sizeD, ncol=m, byrow=TRUE)
X <- (X + runif(m*sizeD)/10)
X[500,1] <- (2*X[500,1]-X[300])/3

# Xprime <- data.frame(X)
# with(Xprime, {
#    scatterplot3d(X1, X2, X3,        # x y and z axis
#                  color="blue", pch=19, # filled blue circles
#                  main="3-D Line w/ Outlier",
#                  xlab="X",
#                  ylab="Y",
#                  zlab="Z")
# })

s3d <- scatterplot3d(X[c(1:499, 501:1000),2], X[c(1:499, 501:1000),3], X[c(1:499, 501:1000),1],        # x y and z axis
                 color="blue", pch=19, # filled blue circles
                 main="Noise 3-D Line w/ Outlier",
                 xlab="X",
                 ylab="Y",
                 zlab="Z")

s3d$points3d(X[500,2], X[500,3], X[500,1], color="green", pch=17)#, cex=8)


#NNorder <- order(NN1, decreasing=TRUE)
# #print(X)
# for(z in 1:sizeD){
# nearK <- order(sqrt(rowSums(sweep(X,2,X[z,])^2)))[1:k]
# AkNN[z,nearK] <- 1
# }

# #NNorder <- order(NN1, decreasing=TRUE)
# forest <- rfrus(X,trees=numtrees, MinParent = k, bagging=0)
# for(z in 1:sizeD){
#   NN1 <- distNN(X[z,], X, forest)
# nearK <- order(NN1,decreasing = TRUE)[1:k]
# AppkNN[z,nearK] <- 1
# }

#NNorder <- order(NN1, decreasing=TRUE)
for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}

#NNorder <- order(NN1, decreasing=TRUE)
forest <- rfrus(X,trees=numtrees, MinParent = k)
for(z in 1:sizeD){
  NN1 <- distNN(X[z,], X, forest)
  AppkNN[z,] = NN1
}

#NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN != 0)
# plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=3", font.main = 4)
# 
# #NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN[,500] != 0)
# plot(cbind(AppkNN[nnzPts,500],AkNN[nnzPts,500]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (1 sample, #500), n=1000, d=3", font.main = 4)

nnzPts <- which(AppkNN != 0)
ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = AppkNN[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\nNoisy Line w/ Outlier, n=1000, d=3, k=10, trees=100\nAll Samples (0 Nearness Omitted)")

#plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
#title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=1", font.main = 4)

#NNorder <- order(NN1, decreasing=TRUE)
nnzPts <- which(t(AppkNN[499:501,]) != 0)
ssd <- data.frame(Distance = t(AkNN[499:501,])[nnzPts], Nearness = t(AppkNN[499:501,])[nnzPts])
groupLabels <- c(rep("499",sizeD), rep("500 (outlier)", sizeD), rep("501",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\nNoisy Line w/ Outlier, n=1000, d=3, k=10, trees=100\nThree Samples (0 Nearness Omitted)")+ geom_jitter()
```

# Swiss Roll

```{r SR, cache = FALSE, echo=FALSE, message=FALSE}
numtrees <- 100
m <- 3
sizeD <- 1000
k = 10

AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)
AppkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

## Remove samples at random for testing.
X <- swissRoll(sizeD/2, size =1, dim3=T)
X <- as.matrix(X)[,2:4]

Xprime <- data.frame(X)
with(Xprime, {
   scatterplot3d(x2, x3, x1,        # x y and z axis
                 color="blue", pch=19, # filled blue circles
                 main="3-D Swiss Roll",
                 xlab="X",
                 ylab="Y",
                 zlab="Z")
})
#X <- as.matrix(X)[,2:4]
#NNorder <- order(NN1, decreasing=TRUE)
# #print(X)
# for(z in 1:sizeD){
# nearK <- order(sqrt(rowSums(sweep(X,2,X[z,])^2)))[1:k]
# AkNN[z,nearK] <- 1
# }

# #NNorder <- order(NN1, decreasing=TRUE)
# forest <- rfrus(X,trees=numtrees, MinParent = k, bagging=0)
# for(z in 1:sizeD){
#   NN1 <- distNN(X[z,], X, forest)
# nearK <- order(NN1,decreasing = TRUE)[1:k]
# AppkNN[z,nearK] <- 1
# }

#NNorder <- order(NN1, decreasing=TRUE)
for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}

#NNorder <- order(NN1, decreasing=TRUE)
forest <- rfrus(X,trees=numtrees, MinParent = k)
for(z in 1:sizeD){
  NN1 <- distNN(X[z,], X, forest)
  AppkNN[z,] = NN1
}

#NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN != 0)
# plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=3", font.main = 4)
# 
# #NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN[,500] != 0)
# plot(cbind(AppkNN[nnzPts,500],AkNN[nnzPts,500]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (1 sample, #500), n=1000, d=3", font.main = 4)

nnzPts <- which(AppkNN != 0)
ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = AppkNN[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point() + labs(title="Distance vs Nearness\nSwiss Roll, n=1000, d=3, k=10, trees=100\nAll Samples (0 Nearness Omitted)")

#plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
#title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=1", font.main = 4)

#NNorder <- order(NN1, decreasing=TRUE)
nnzPts <- which(t(AppkNN[499:501,]) != 0)
ssd <- data.frame(Distance = t(AkNN[499:501,])[nnzPts], Nearness = t(AppkNN[499:501,])[nnzPts])
groupLabels <- c(rep("499",sizeD), rep("500", sizeD), rep("501",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point() + labs(title="Distance vs Nearness\nSwiss Roll, n=1000, d=3, k=10, trees=100\nThree Samples (0 Nearness Omitted)")+ geom_jitter()
```


# Swiss Roll w/ Outlier

```{r SRwOutlier, cache = FALSE, echo=FALSE, message=FALSE}
numtrees <- 100
m <- 3
sizeD <- 1000
k = 10

AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)
AppkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

## Remove samples at random for testing.
X <- swissRoll(sizeD/2, size =1, dim3=T)
X <- as.matrix(X)[,2:4]
X[500,2] <- -1.5
X[500,3] <- 4
X[500,1] <- 1

# Xprime <- data.frame(X)
# with(Xprime, {
#    scatterplot3d(x2, x3, x1,        # x y and z axis
#                  color="blue", pch=19, # filled blue circles
#                  main="3-D Swiss Roll w/ Outlier",
#                  xlab="X",
#                  ylab="Y",
#                  zlab="Z")
# })

s3d <- scatterplot3d(X[c(1:499, 501:1000),2], X[c(1:499, 501:1000),3], X[c(1:499, 501:1000),1],        # x y and z axis
                 color="blue", pch=19, # filled blue circles
                 main="3-D Swiss Roll w/ Outlier",
                 xlab="X",
                 ylab="Y",
                 zlab="Z")

s3d$points3d(X[500,2], X[500,3], X[500,1], color="red", pch=17, cex=1)


#$points3d(seq(10, 20, 2), seq(85, 60, -5), seq(60, 10, -10),
#col = "red", type = "h", pch = 8)

#X <- as.matrix(X)[,2:4]
#NNorder <- order(NN1, decreasing=TRUE)
# #print(X)
# for(z in 1:sizeD){
# nearK <- order(sqrt(rowSums(sweep(X,2,X[z,])^2)))[1:k]
# AkNN[z,nearK] <- 1
# }

# #NNorder <- order(NN1, decreasing=TRUE)
# forest <- rfrus(X,trees=numtrees, MinParent = k, bagging=0)
# for(z in 1:sizeD){
#   NN1 <- distNN(X[z,], X, forest)
# nearK <- order(NN1,decreasing = TRUE)[1:k]
# AppkNN[z,nearK] <- 1
# }

#NNorder <- order(NN1, decreasing=TRUE)
for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}

#NNorder <- order(NN1, decreasing=TRUE)
forest <- rfrus(X,trees=numtrees, MinParent = k)
for(z in 1:sizeD){
  NN1 <- distNN(X[z,], X, forest)
  AppkNN[z,] = NN1
}

#NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN != 0)
# plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=3", font.main = 4)
# 
# #NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN[,500] != 0)
# plot(cbind(AppkNN[nnzPts,500],AkNN[nnzPts,500]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (1 sample, #500), n=1000, d=3", font.main = 4)

nnzPts <- which(AppkNN != 0)
ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = AppkNN[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\nSwiss Roll w/ Outlier, n=1000, d=3, k=10, trees=100\nAll Samples (0 Nearness Omitted)")

#plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
#title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=1", font.main = 4)

#NNorder <- order(NN1, decreasing=TRUE)
nnzPts <- which(t(AppkNN[499:501,]) != 0)
ssd <- data.frame(Distance = t(AkNN[499:501,])[nnzPts], Nearness = t(AppkNN[499:501,])[nnzPts])
groupLabels <- c(rep("499",sizeD), rep("500 (outlier)", sizeD), rep("501",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point() + labs(title="Distance vs Nearness\nSwiss Roll w/ Outlier, n=1000, d=3, k=10, trees=100\nThree Samples (0 Nearness Omitted)")+ geom_jitter()
```


# Sparse Swiss Roll w/ Outlier

```{r SparseSRwOutlier, cache = FALSE, echo=FALSE, message=FALSE}
numtrees <- 100
m <- 3
sizeD <- 900
k = 3

AkNN <- matrix(0, nrow=sizeD, ncol=sizeD)
AppkNN <- matrix(0, nrow=sizeD, ncol=sizeD)

## Remove samples at random for testing.
X <- swissRoll(sizeD/2, size =1, dim3=T)
X <- as.matrix(X)[,2:4]
X[150,2] <- -1.5
X[150,3] <- 4
X[150,1] <- 1

# Xprime <- data.frame(X)
# with(Xprime, {
#    scatterplot3d(x2, x3, x1,        # x y and z axis
#                  color="blue", pch=19, # filled blue circles
#                  main="3-D Swiss Roll w/ Outlier",
#                  xlab="X",
#                  ylab="Y",
#                  zlab="Z")
# })

s3d <- scatterplot3d(X[c(1:149, 151:300),2], X[c(1:149, 151:300),3], X[c(1:149, 151:300),1],        # x y and z axis
                 color="blue", pch=19, # filled blue circles
                 main="Sparse 3-D Swiss Roll w/ Outlier",
                 xlab="X",
                 ylab="Y",
                 zlab="Z")

s3d$points3d(X[150,2], X[150,3], X[150,1], color="red", pch=17, cex=1)


#$points3d(seq(10, 20, 2), seq(85, 60, -5), seq(60, 10, -10),
#col = "red", type = "h", pch = 8)

#X <- as.matrix(X)[,2:4]
#NNorder <- order(NN1, decreasing=TRUE)
# #print(X)
# for(z in 1:sizeD){
# nearK <- order(sqrt(rowSums(sweep(X,2,X[z,])^2)))[1:k]
# AkNN[z,nearK] <- 1
# }

# #NNorder <- order(NN1, decreasing=TRUE)
# forest <- rfrus(X,trees=numtrees, MinParent = k, bagging=0)
# for(z in 1:sizeD){
#   NN1 <- distNN(X[z,], X, forest)
# nearK <- order(NN1,decreasing = TRUE)[1:k]
# AppkNN[z,nearK] <- 1
# }

#NNorder <- order(NN1, decreasing=TRUE)
for(z in 1:sizeD){
AkNN[z,] <- sqrt(rowSums(sweep(X,2,X[z,])^2))
}

#NNorder <- order(NN1, decreasing=TRUE)
forest <- rfrus(X,trees=numtrees, MinParent = k)
for(z in 1:sizeD){
  NN1 <- distNN(X[z,], X, forest)
  AppkNN[z,] = NN1
}

#NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN != 0)
# plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=3", font.main = 4)
# 
# #NNorder <- order(NN1, decreasing=TRUE)
# nnzPts <- which(AppkNN[,500] != 0)
# plot(cbind(AppkNN[nnzPts,500],AkNN[nnzPts,500]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
# title(main = "Nearness vs Distance (1 sample, #500), n=1000, d=3", font.main = 4)

nnzPts <- which(AppkNN != 0)
ssd <- data.frame(Distance = AkNN[nnzPts], Nearness = AppkNN[nnzPts])
ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point()+ labs(title="Distance vs Nearness\nSwiss Roll w/ Outlier, n=300, d=3, k=3, trees=100\nAll Samples (0 Nearness Omitted)")

#plot(cbind(AppkNN[nnzPts],AkNN[nnzPts]), xlab="Nearness", ylab="Distance", pch=19, cex=1)
#title(main = "Nearness vs Distance (all 1000 samples), n=1000, d=1", font.main = 4)

#NNorder <- order(NN1, decreasing=TRUE)
nnzPts <- which(t(AppkNN[149:151,]) != 0)
ssd <- data.frame(Distance = t(AkNN[149:151,])[nnzPts], Nearness = t(AppkNN[149:151,])[nnzPts])
groupLabels <- c(rep("149",sizeD), rep("150 (outlier)", sizeD), rep("151",sizeD ))[nnzPts]
ssd[["Sample"]] <- groupLabels
print(ggplot(aes(x = Nearness, y = Distance, color = Sample), data = ssd) + geom_point() + labs(title="Distance vs Nearness\nSparse Swiss Roll w/ Outlier, n=300, d=3, k=3, trees=100\nThree Samples (0 Nearness Omitted)")+ geom_jitter())
```

