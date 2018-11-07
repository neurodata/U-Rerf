source('../rfr_us.R')
library(ggplot2)
library(reshape2)

# number of trees for forest
numtrees <- 100
# the 'k' of k nearest neighbors


# create a sizeD by m synthetic dataset
X <- as.matrix(iris[,1:4])
k <- sqrt(nrow(X))

# create a similarity matrix using urerf
similarityMatrix <- urerf(X, numtrees, k)$similarityMatrix
#similarityMatrix <- data.frame(similarityMatrix)

meltedSM <- melt(similarityMatrix)

png(file="results/ex11_iris.png")
#heatmap(similarityMatrix$similarityMatrix, Rowv=NA, Colv=NA, col=heat.colors(256), symm=TRUE, main="Heatmap of Iris Dataset Similarity Matrix")
p <- ggplot(data =meltedSM , aes(x=Var1, y=Var2, fill=value)) + geom_tile()
p <- p +  guides(fill=guide_legend(title="Similarity"))
print(p)
dev.off()

