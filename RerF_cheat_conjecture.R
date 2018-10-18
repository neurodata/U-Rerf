# Minh and Carey's conjecture using RerF
library(stats)
library(rgl)
library(MASS)
library(Matrix)
#library(randomForest)
library(scatterplot3d)
source("rfr_us.R")

num_of_points=1000
num_of_trees = 20
knn=10

set.seed(1)

t=runif(num_of_points, min = 0, max = 1)
t=sort(t)
x1=t^2
x2=2*t*(1-t)
x3=(1-t)^2
data = cbind(x1, x2, x3)

f_signal=urerf(data, num_of_trees, knn)
Rerf_W_signal=f_signal$similarityMatrix

# g_signal=randomForest(data, ntree=num_of_trees, keep.forest=FALSE, proximity=TRUE)
# W_signal=g_signal$proximity

####This is all nice and good, but now we are going to make it into high dim data
# In this case dim is ?
high_dim = 9
matrix_of_0 = matrix(rep(0, num_of_points*(high_dim-3)), nrow = num_of_points, ncol = high_dim -3)
high_dim_data = cbind(data, matrix_of_0)

dim <- high_dim

# change gaus to false to use uniform noise.  change gaus to true to use gaussian noise.
gaus <- FALSE
if(gaus){
	cov_matrix = matrix(rep(0, dim*dim), nrow = dim, ncol = dim)
	diag(cov_matrix) = c(0, 0, 0, rep(70, high_dim-3))
	Sig1 = cov_matrix
	set.seed(2)
	noise = mvrnorm(n = num_of_points, (rep(0, high_dim)), Sig1, tol = 1e-7, empirical = FALSE, EISPACK = FALSE)
} else{
	noise = matrix(0, nrow = num_of_points, ncol = high_dim)
	for(j in 4:high_dim){
		noise[,j] <- runif(num_of_points,0,1)
	}
}

high_dim_noise_data = high_dim_data + noise
rownames(high_dim_noise_data) <- c()
colnames(high_dim_noise_data) <- c()

# g_noise=randomForest(noise[, 4:high_dim], ntree=num_of_trees, keep.forest=FALSE, proximity=TRUE)
# W_noise=g_noise$proximity

f_noise <- urerf(noise[, 4:high_dim], num_of_trees, knn)
Rerf_W_noise=f_noise$similarityMatrix

#######This is the real cheat
# Rerf_W_cheat_sum = matrix(rep(0, num_of_points*num_of_points), nrow = num_of_points, ncol = num_of_points)
# 
# num_of_cheat_trees=num_of_trees
# # cheat = sample(c(1, 2, 3), size = num_of_cheat_trees, replace = TRUE, prob = c(1/choose(9, 3), choose(6, 3)/choose(9, 3), 1- 1/choose(9, 3)-choose(6, 3)/choose(9, 3) ) )
# cheat = sample(c(1, 2), size = num_of_cheat_trees, replace = TRUE, prob = c(1/3, 2/3) )
# i=1
# tree=0
# while(i <= num_of_cheat_trees){
#   if(cheat[i] == 1){
#     #g=randomForest(data, mtry=3, ntree=1, keep.forest=FALSE, proximity=TRUE)
#     f=urerf(data, 1, knn)
#     tree = tree+1
#     Rerf_W_cheat_sum=Rerf_W_cheat_sum + f$similarityMatrix
#   }
#   if(cheat[i] == 2){
#     #g=randomForest(noise[, 4:high_dim], mtry=3, ntree=1, keep.forest=FALSE, proximity=TRUE)
#     f<- urerf(noise[, 4:high_dim], 1, knn)
#     tree = tree+1
#     Rerf_W_cheat_sum=Rerf_W_cheat_sum + f$similarityMatrix
#   }
#   # else do nothing (ignore the case that is the mix of the signal and the noise)
#   i=i+1
# }
# 
# Rerf_W_cheat_ave = Rerf_W_cheat_sum/tree

Rerf_W_mix=(3/dim)*Rerf_W_signal + ((dim-3)/dim)*Rerf_W_noise

Rerf_W_uncheat = urerf(high_dim_noise_data, 20, K=10)$similarityMatrix


print("Rerf uncheat")
print(Rerf_W_uncheat[1,1:100])
#print("Rerf cheat")
#print(Rerf_W_cheat_ave[20,1:100])
print("Rerf mix")
print(Rerf_W_mix[20,1:100])

Rerf_uncheat_result=cmdscale(sqrt(1-Rerf_W_uncheat), k=3)
#Rerf_cheat_result=cmdscale(sqrt(1-Rerf_W_cheat_ave), k=3)
Rerf_mix_result=cmdscale(sqrt(1-Rerf_W_mix), k=3)


# plot3d(uncheat_result[1:100,1], uncheat_result[1:100,2], uncheat_result[1:100,3], col='red', size=3, add = TRUE)
# plot3d(uncheat_result[101:200,1], uncheat_result[101:200,2], uncheat_result[101:200,3], col='green', size=3, add=TRUE)
# plot3d(uncheat_result[201:300,1], uncheat_result[201:300,2], uncheat_result[201:300,3], col='orange', size=3, add=TRUE)
# plot3d(uncheat_result[301:400,1], uncheat_result[301:400,2], uncheat_result[301:400,3], col='blue', size=3, add=TRUE)
# plot3d(uncheat_result[401:500,1], uncheat_result[401:500,2], uncheat_result[401:500,3], size=3, add=TRUE)
# #
# plot3d(cheat_result[1:100,1], cheat_result[1:100,2], cheat_result[1:100,3], col='red', size=3, add = TRUE)
# plot3d(cheat_result[101:200,1], cheat_result[101:200,2], cheat_result[101:200,3], col='green', size=3, add=TRUE)
# plot3d(cheat_result[201:300,1], cheat_result[201:300,2], cheat_result[201:300,3], col='orange', size=3, add=TRUE)
# plot3d(cheat_result[301:400,1], cheat_result[301:400,2], cheat_result[301:400,3], col='blue', size=3, add=TRUE)
# plot3d(cheat_result[401:500,1], cheat_result[401:500,2], cheat_result[401:500,3], size=3, add=TRUE)
# #
# plot3d(mix_result[1:100,1], mix_result[1:100,2], mix_result[1:100,3], col='red', size=3, add = TRUE)
# plot3d(mix_result[101:200,1], mix_result[101:200,2], mix_result[101:200,3], col='green', size=3, add=TRUE)
# plot3d(mix_result[201:300,1], mix_result[201:300,2], mix_result[201:300,3], col='orange', size=3, add=TRUE)
# plot3d(mix_result[301:400,1], mix_result[301:400,2], mix_result[301:400,3], col='blue', size=3, add=TRUE)
# plot3d(mix_result[401:500,1], mix_result[401:500,2], mix_result[401:500,3], size=3, add=TRUE)


plot3d(Rerf_uncheat_result[1:200,1], Rerf_uncheat_result[1:200,2], Rerf_uncheat_result[1:200,3], col='red', size=3, xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3) , zlim = c(-0.3, 0.3))
plot3d(Rerf_uncheat_result[201:400,1], Rerf_uncheat_result[201:400,2], Rerf_uncheat_result[201:400,3], col='green', size=3, add=TRUE)
plot3d(Rerf_uncheat_result[401:600,1], Rerf_uncheat_result[401:600,2], Rerf_uncheat_result[401:600,3], col='orange', size=3, add=TRUE)
plot3d(Rerf_uncheat_result[601:800,1], Rerf_uncheat_result[601:800,2], Rerf_uncheat_result[601:800,3], col='blue', size=3, add=TRUE)
plot3d(Rerf_uncheat_result[801:1000,1], Rerf_uncheat_result[801:1000,2], Rerf_uncheat_result[801:1000,3], size=3, add=TRUE)
#
# plot3d(Rerf_cheat_result[1:200,1], Rerf_cheat_result[1:200,2], Rerf_cheat_result[1:200,3], col='red', size=3, xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3) , zlim = c(-0.3, 0.3))
# plot3d(Rerf_cheat_result[201:400,1], Rerf_cheat_result[201:400,2], Rerf_cheat_result[201:400,3], col='green', size=3, add=TRUE)
# plot3d(Rerf_cheat_result[401:600,1], Rerf_cheat_result[401:600,2], Rerf_cheat_result[401:600,3], col='orange', size=3, add=TRUE)
# plot3d(Rerf_cheat_result[601:800,1], Rerf_cheat_result[601:800,2], Rerf_cheat_result[601:800,3], col='blue', size=3, add=TRUE)
# plot3d(Rerf_cheat_result[801:1000,1], Rerf_cheat_result[801:1000,2], Rerf_cheat_result[801:1000,3], size=3, add=TRUE)
# # #
plot3d(Rerf_mix_result[1:200,1], Rerf_mix_result[1:200,2], Rerf_mix_result[1:200,3], col='red', size=3, xlim = c(-0.3, 0.3), ylim = c(-0.3, 0.3) , zlim = c(-0.3, 0.3))
plot3d(Rerf_mix_result[201:400,1], Rerf_mix_result[201:400,2], Rerf_mix_result[201:400,3], col='green', size=3, add=TRUE)
plot3d(Rerf_mix_result[401:600,1], Rerf_mix_result[401:600,2], Rerf_mix_result[401:600,3], col='orange', size=3, add=TRUE)
plot3d(Rerf_mix_result[601:800,1], Rerf_mix_result[601:800,2], Rerf_mix_result[601:800,3], col='blue', size=3, add=TRUE)
plot3d(Rerf_mix_result[801:1000,1], Rerf_mix_result[801:1000,2], Rerf_mix_result[801:1000,3], size=3, add=TRUE)

plot3d(Rerf_mix_result[,1], Rerf_mix_result[,2], Rerf_mix_result[,3], size=3, add=TRUE)
#
# 
