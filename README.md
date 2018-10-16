# U-Rerf
The unsupervised version of Rerf.  Currently U-Rerf creates either an NxN similarity matrix which shows the similarity between each of the N samples to each of the other N-1 samples or U-Rerf can create an NxN distance matrix showing the distance between each of the N samples to each of the other N samples.

## U-Rerf Setup
* Clone the repository
* source the U-Rerf functions -- `source('rfr_us.R')`
* To run the examples
  * install ggplot2 -- `install.packages("ggplot2")`
  * install scatterplot3d -- `install.packages("scatterplot3d")`
  * install mass -- `install.packages("MASS")`
  
## Usage
To train a forest use the urerf or urerfDepth function which is loaded when the rfr_us.R file is sourced.  The function takes three inputs:
1. X - which is an N row by d column matrix of N samples each of d dimensions.  X must be numeric.
1. numTrees - the number of trees used to create the forest.  100 is the default.
1. K or Depth - this is really the minparent used in creating the forest but should be thought of as the K in K nearest neighbors.  If urerfDepth is chosen then the third parameter defines the maximum allowable depth in the forest.

## Examples
Basic usage of U-Rerf.
```
# make sure to put the correct path to rfr_us.R
source('rfr_us.R')

# number of trees for forest
numtrees <- 100
# the 'k' of k nearest neighbors, this parameter is equivalent to minparent found in Random Forests
k <- 10
# set max depth
depth <- 4

# create a sizeD by m synthetic dataset
X <- as.matrix(iris[,1:4])

# create a urerf structure which includes the forest and similarity matrix
urerfStructure <- urerf(X, numtrees, K=k)
urerfStructure <- urerf(X, numtrees, depth=depth)
```

Example use cases can be found in the examples directory.
