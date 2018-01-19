source("rfr_us.R")

context("output")
test_that("bestCutForFeature", {

combinedFeature <- rep(5, 20)
expect_true(is.null(bestCutForFeature(combinedFeature)))
a <- 5
b <- 10
combinedFeature <- c(rep(a,20), rep(b,20))
expect_equal(bestCutForFeature(combinedFeature), c((a+b)/2, 0))
combinedFeature <- c(runif(10), runif(20)+2)
result <- bestCutForFeature(combinedFeature)
expect_equal(result[1], mean(sort(combinedFeature)[10:11]))
expect_true(result[2] > 0)
})


test_that("output is the right size and type", {

X <- as.matrix(iris[,1:4])
similarityMatrix <- createSimilarityMatrix (X, 100, 10)
#similarityMatrix <- X
expect_equal(nrow(similarityMatrix), 150)
expect_equal(ncol(similarityMatrix), 150)
})
