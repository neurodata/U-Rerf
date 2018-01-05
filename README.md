# U-Rerf
The unsupervised version of Rerf.  Currently U-Rerf only creates an NxN similarity matrix which shows the similarity between each of the N samples to each of the other N-1 samples.

## U-Rerf Setup
* Clone the repository
* source the U-Rerf functions -- `source('rfr_us.R')`
* To run the examples
  * install ggplot2 -- `install.packages("ggplot2")`
  * install scatterplot3d -- `install.packages("scatterplot3d")`
  
## Usage
To create a similarity matrix use the createSimilarityMatrix function which is loaded when the rfr_us.R file is sourced.  The function takes three inputs:
1. X - which is an N row by d column matrix of N samples each of d dimensions.  X must be numeric.
1. numTrees - the number of trees used to create the forest.  100 is the default.
1. K - this is really the minparent used in creating the forest but should be thought of as the K in K nearest neighbors.
