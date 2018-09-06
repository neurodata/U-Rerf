source('../rfr_us.R')
library(ggplot2)


X <- as.matrix(iris[,1:4])

sM <- createSimilarityMatrix(X)

cluss1 <- cluster(sM,3,"average")
table(cluss1, iris$Species)

cluss2 <- cluster(sM,3,"mcquitty")
table(cluss2, iris$Species)

cluss3 <- cluster(sM,3,"kmeans")
table(cluss3, iris$Species)

cluss4 <- cluster(sM,3,"medoids")
table(cluss4, iris$Species)
