# Minh and Carey's conjecture using RerF
library(stats)
library(ggplot2)
library(MASS)
library(Matrix)
library(randomForest)
library(scatterplot3d)
library(gridExtra)
source("rfr_us.R")


pdf(file="noiseTest.pdf")


num_of_points=1000
num_of_trees =  200
knn=33
maxNoiseDims <-4

set.seed(9)

normalizeData <- function(X){
		X <- sweep(X, 2, apply(X, 2, min), "-")
		sweep(X, 2, apply(X, 2, max), "/")
	}

fontSize <- 6
leg <- theme(legend.text = element_text(size = fontSize), legend.title=element_text(size = fontSize), plot.title = element_text(size =fontSize,  face="bold"), axis.title.x = element_text(size=fontSize), axis.text.x = element_text(size=fontSize), axis.title.y = element_text(size=fontSize), axis.text.y = element_text(size=fontSize), strip.text.x = element_text(size=fontSize))

NvDPlot <- function(SMfromUrerf, ActualDistance, nnz, numPoints, plotTitle){
	midPoint <- ceiling(numPoints/2)
	ssd <- data.frame(Distance = ActualDistance[,midPoint][nnz], Nearness = SMfromUrerf[,midPoint][nnz])
	p <- ggplot(aes(x = Nearness, y = Distance), data = ssd) + geom_point() + labs(title=plotTitle)
	p <- p + scale_y_log10() + leg
	p + xlab("Similarity") + ylab("Euclidean Distance")
}


nonMetricPlot <- function(SMfromUrerf, ActualDistance, nnz, numPoints, plotTitle){
	midPoint <- ceiling(numPoints/2)
	Rerf_W_signal <- SMfromUrerf / 1.01
	diag(Rerf_W_signal) <- 1
	fit <- isoMDS(as.dist(1-Rerf_W_signal), k=2)
	fit <- data.frame(fit)
	ggplot(aes(x = points.1, y = points.2), data = fit) + geom_point(color = rainbow(numPoints)) + labs(title=plotTitle)+ leg

}


metricPlot <- function(SMfromUrerf, ActualDistance, nnz, numPoints, plotTitle){
	midPoint <- ceiling(numPoints/2)
	Rerf_W_signal <- SMfromUrerf / 1.01
	diag(Rerf_W_signal) <- 1
	fit <- cmdscale(as.dist(1-Rerf_W_signal), eig=FALSE, k=2)
	fit <- data.frame(fit)
	ggplot(aes(x = X1, y = X2), data = fit) + geom_point(color = rainbow(numPoints)) + labs(title=plotTitle)+ leg

}



#Create Data
t=runif(num_of_points, min = 0, max = 1)
t=sort(t)
x1=t^2
x2=2*t*(1-t)
x3=(1-t)^2
data = cbind(x1, x2, x3)

#data <- normalizeData(data)

dMatrixActual <- as.matrix(dist(data))
nnzPts <- which(dMatrixActual[,500] != 0)




#### Run experiments on Urerf minus scaling
Alg <- "Urerf, "
for (gaus in c(TRUE, FALSE)){
	plotNum <- 1
	plots = list()
	for(z in 0:maxNoiseDims){
		print(paste("\n\n\nHW curve with ", z , " UNIFORM noisy dims.  Urerf Nearness vs Distance\n\n\n"))
		high_dim = z+3
		if(z){
			matrix_of_0 = matrix(rep(0, num_of_points*(high_dim-3)), nrow = num_of_points, ncol = high_dim -3)
			high_dim_data = cbind(data, matrix_of_0)
		}else{
			high_dim_data <- data
		}

		dim <- high_dim

		# change gaus to false to use uniform noise.  change gaus to true to use gaussian noise.
		if(z){
			if(gaus){
				cov_matrix = matrix(rep(0, dim*dim), nrow = dim, ncol = dim)
				diag(cov_matrix) = c(0, 0, 0, rep(70, high_dim-3))
				Sig1 = cov_matrix
				noise = mvrnorm(n = num_of_points, (rep(0, high_dim)), Sig1, tol = 1e-7, empirical = FALSE, EISPACK = FALSE)
			} else{
				noise = matrix(0, nrow = num_of_points, ncol = high_dim)
				for(j in 4:high_dim){
					noise[,j] <- runif(num_of_points,0,1)
				}
			}
			high_dim_noise_data <- high_dim_data + noise
		}else{
			high_dim_noise_data <- high_dim_data 
		}

		rownames(high_dim_noise_data) <- c()
		colnames(high_dim_noise_data) <- c()

		Rerf_W_uncheat = urerf(normalizeData(high_dim_noise_data), num_of_trees, K=knn)$similarityMatrix

		if(gaus){
			nType <- "Gaussian"
		}else{
			nType <- "Uniform"
		}

		plots[[plotNum]] <- NvDPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg,z, " ", nType, " noisy dims\nNoisy Nearness Vs 3-D Distance"))
		plotNum <- plotNum+1
		plots[[plotNum]] <- nonMetricPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg, z, " ", nType, " noisy dims\nNonmetric MDS"))
		plotNum <- plotNum+1
		plots[[plotNum]] <- metricPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg, z, " ", nType, " noisy dims\nMetric MDS"))
		plotNum <- plotNum+1
	}
	print(do.call(grid.arrange, c(plots, ncol=3)))
}




#### Run experiments on RF
Alg <- "RF, "
for (gaus in c(TRUE, FALSE)){
	plotNum <- 1
	plots = list()
	for(z in 0:maxNoiseDims){
		print(paste("\n\n\nHW curve with ", z , " UNIFORM noisy dims.  Urerf Nearness vs Distance\n\n\n"))
		high_dim = z+3
		if(z){
			matrix_of_0 = matrix(rep(0, num_of_points*(high_dim-3)), nrow = num_of_points, ncol = high_dim -3)
			high_dim_data = cbind(data, matrix_of_0)
		}else{
			high_dim_data <- data
		}

		dim <- high_dim

		# change gaus to false to use uniform noise.  change gaus to true to use gaussian noise.
		if(z){
			if(gaus){
				cov_matrix = matrix(rep(0, dim*dim), nrow = dim, ncol = dim)
				diag(cov_matrix) = c(0, 0, 0, rep(70, high_dim-3))
				Sig1 = cov_matrix
				noise = mvrnorm(n = num_of_points, (rep(0, high_dim)), Sig1, tol = 1e-7, empirical = FALSE, EISPACK = FALSE)
			} else{
				noise = matrix(0, nrow = num_of_points, ncol = high_dim)
				for(j in 4:high_dim){
					noise[,j] <- runif(num_of_points,0,1)
				}
			}
			high_dim_noise_data <- high_dim_data + noise
		}else{
			high_dim_noise_data <- high_dim_data 
		}

		rownames(high_dim_noise_data) <- c()
		colnames(high_dim_noise_data) <- c()

		Rerf_W_uncheat = randomForest(normalizeData(high_dim_noise_data), ntree=num_of_trees)$proximity

		if(gaus){
			nType <- "Gaussian"
		}else{
			nType <- "Uniform"
		}

		plots[[plotNum]] <- NvDPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg,z, " ", nType, " noisy dims\nNoisy Nearness Vs 3-D Distance"))
		plotNum <- plotNum+1
		plots[[plotNum]] <- nonMetricPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg, z, " ", nType, " noisy dims\nNonmetric MDS"))
		plotNum <- plotNum+1
		plots[[plotNum]] <- metricPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg, z, " ", nType, " noisy dims\nMetric MDS"))
		plotNum <- plotNum+1
	}
	print(do.call(grid.arrange, c(plots, ncol=3)))
}



#### Run experiments on Normalized Urerf
Alg <- "UrerFNorm, "
for (gaus in c(TRUE, FALSE)){
	plotNum <- 1
	plots = list()
	for(z in 0:maxNoiseDims){
		print(paste("\n\n\nHW curve with ", z , " UNIFORM noisy dims.  Urerf Nearness vs Distance\n\n\n"))
		high_dim = z+3
		if(z){
			matrix_of_0 = matrix(rep(0, num_of_points*(high_dim-3)), nrow = num_of_points, ncol = high_dim -3)
			high_dim_data = cbind(data, matrix_of_0)
		}else{
			high_dim_data <- data
		}

		dim <- high_dim

		# change gaus to false to use uniform noise.  change gaus to true to use gaussian noise.
		if(z){
			if(gaus){
				cov_matrix = matrix(rep(0, dim*dim), nrow = dim, ncol = dim)
				diag(cov_matrix) = c(0, 0, 0, rep(70, high_dim-3))
				Sig1 = cov_matrix
				noise = mvrnorm(n = num_of_points, (rep(0, high_dim)), Sig1, tol = 1e-7, empirical = FALSE, EISPACK = FALSE)
			} else{
				noise = matrix(0, nrow = num_of_points, ncol = high_dim)
				for(j in 4:high_dim){
					noise[,j] <- runif(num_of_points,0,1)
				}
			}
			high_dim_noise_data <- high_dim_data + noise
		}else{
			high_dim_noise_data <- high_dim_data 
		}

		rownames(high_dim_noise_data) <- c()
		colnames(high_dim_noise_data) <- c()

		Rerf_W_uncheat = urerfNoNormalize(normalizeData(high_dim_noise_data), num_of_trees, K=knn)$similarityMatrix

		if(gaus){
			nType <- "Gaussian"
		}else{
			nType <- "Uniform"
		}

		plots[[plotNum]] <- NvDPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg,z, " ", nType, " noisy dims\nNoisy Nearness Vs 3-D Distance"))
		plotNum <- plotNum+1
		plots[[plotNum]] <- nonMetricPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg, z, " ", nType, " noisy dims\nNonmetric MDS"))
		plotNum <- plotNum+1
		plots[[plotNum]] <- metricPlot(Rerf_W_uncheat, dMatrixActual,nnzPts, num_of_points, paste(Alg, z, " ", nType, " noisy dims\nMetric MDS"))
		plotNum <- plotNum+1
	}
	print(do.call(grid.arrange, c(plots, ncol=3)))
}

dev.off()
