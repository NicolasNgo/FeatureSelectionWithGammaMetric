###################################################################################
#                       GAMMA-METRIC COMPUTATION FUNCTIONS                        #
###################################################################################

# Set of functions to compute the gamma-metric value of a feature X (or a set of features)
# With regard to a target feature 'class'

## Required packages ----------------------------------------------------------------------
if(!require(corpcor)){
  # Estimate variance-covariance matrix with the shrinkage estimator (v1.6.10)
  devtools::install_version('corpcor', version = '1.6.10')
  library(corpcor)
}
  
if(!require(ggplot2)){
  # Graphics
  devtools::install_version('ggplot2', version = '3.4.3')
  library(ggplot2)
}

## Function to compute the norm of a vector x -----------------------------------------------
norm <- function(x){
  ## Inputs: 
  ##      x: The vector we want to compute the norm 2
  ## Outputs: 
  ##      The value of the norm 2
  return(sqrt(sum(x^2)))
}

## Function to compute the normalisation of a vector x --------------------------------------
normalize <- function(x){
  ## Inputs: 
  ##      x: The vector we want to normalize
  ## Outputs: The normalized x
  return(x/norm(x))
}

## Function returning a vector divided by the sum of its coefficients -----------------------
divsum <- function(x){
  ## Inputs:
  ##      x: The vector we want to divide by its sum
  ## Outputs: 
  ##      The vector divided by the sum of its coefficients
  return(x/sum(x))
}

## Function to compute the coordinates of the contour of the ellipsoids on two dimensions ---
##  This function is mostly useful for plots                                   
ellipse <- function(WkM, mu){
  ## Inputs: 
  ##      WkM: The covariance matrix of class k
  ##       mu: The mean vector of each ellipsoid
  ## Outputs:
  ##        xt: coordinates on axis x of the contour of the ellipsoids on two dimensions
  ##        yt: coordinates on axis y of the contour of the ellipsoids on two dimensions
  ## We used a formula related to the representation of a covariance matrix with an ellipse
  
  # Parameters of the ellipsoid for a 2x2 variance matrix
  a <- WkM[1, 1]
  b <- WkM[1, 2]
  c <- WkM[2, 2]
  
  # eigen values of the covariance matrix. These are also the length of the semi-major and semi-minor axes of the ellipse
  lambda1 <- (a + c)/2 + sqrt(((a - c)/2)^2 + b^2)
  lambda2 <- (a + c)/2 - sqrt(((a - c)/2)^2 + b^2)
  
  # Angle of rotation of the axes of the ellipse
  # b == 0 mean no covariance, features are independents, no rotation
  theta <- ifelse(b != 0, atan2(lambda1 - a, b), ifelse(a < c, pi/2, 0))
  
  # angles t
  t <- seq(0, 2*pi, length.out = 250)
  xt <- sqrt(lambda1)*cos(theta)*cos(t) - sqrt(lambda2)*sin(theta)*sin(t) + mu[1]
  yt <- sqrt(lambda1)*sin(theta)*cos(t) + sqrt(lambda2)*cos(theta)*sin(t) + mu[2]
  
  return(cbind(xt, yt))
}

## Function to compute the distance between k1 and k2 ---------------------------------------
d_k1k2 <- function(mu1, mu2, eigen.values1, eigen.values2, eigen.vectors1, eigen.vectors2){
  ## Inputs: 
  ##      mu1: Coordinates of the mean points of class 1
  ##      mu2: Coordinates of the mean points of class 2
  ##      eigen.values1: vectors of eigen values of the variance-covariance matrix of class 1
  ##      eigen.values2: vectors of eigen values of the variance-covariance matrix of class 2
  ##      eigen.vectors1: matrix of eigen vectors of the variance-covariance matrix of class 1
  ##      eigen.vectors2: matrix of eigen vectors of the variance-covariance matrix of class 2
  
  ## Outputs:
  ##      d_k1k2: distance between ellipses of class 1 and class 2 (on the mean-mean vector with regard to the overlapping)
  
  # Vector mu1mu2
  mu1mu2 <- mu2 - mu1
  
  # Coordinates of the vector mu1mu2 in the ellipse basis
  mu1_coord <- solve(eigen.vectors1) %*% as.matrix(normalize(mu1mu2))
  mu2_coord <- solve(eigen.vectors2) %*% as.matrix(normalize(mu1mu2))
  
  # Compute the normalization factor alpha
  alpha_k1k2 <- sqrt(sum(eigen.values1)) + sqrt(sum(eigen.values2))
  
  # Compute the distances between the center and the edge of the ellipsoid (on mu1mu2's direction)
  dk1 <- 1/sqrt(sum((mu1_coord^2)/(abs(eigen.values1)), na.rm = TRUE))
  dk2 <- 1/sqrt(sum((mu2_coord^2)/(abs(eigen.values2)), na.rm = TRUE))
  
  # Compute d_k1k2
  d_k1k2 <- (1/alpha_k1k2)*(norm(mu1mu2) - (dk1 + dk2))
  
  return(d_k1k2)
}

## Function to compute the gamma-metric ------------------------------------------------------
gammaMetric <- function(X, class, covEstim = c('empiric', 'shrink', 'mixte'), plot = F){
  ## Inputs: 
  ##      X: The matrix of features to compute the gamma-metric value on. Should NOT include class
  ##  class: The target feature
  ##  cov.estim: Which type of variance-covariance matrix estimation should we use ? 
  ##           - 'empiric': operate the empirical estimation of the variance covariance matrix
  ##           - 'shrink': operate the shrinkage estimation with corpcor::cov.shrink function
  ##           - 'mixte': operate the empiric estimation if the number of rows is superiror to the number of features, shrink otherwise
  ##   plot: Should we plot the ellipsoids created (TRUE/FALSE) Only works for 2 features.
  ## Output: 
  ##    The gamma-metric value for feature(s) X
  
  # Matrix X
  X <- as.matrix(X)
  
  # Number of rows
  n <- nrow(X)
  
  # List of the k different classes
  muL <- list()                 # List created to store the coordinates of mu with regards to the class
  eigen.valuesL <- list()       # List created to store the eigen values of the classes
  WkL <- list()                 # List created to store Wk matrix
  eigen.vectorsL <- list()      # List created to store the eigen vectors of the classes
  k <- length(unique(class))    # Number of classes
  length(muL) = length(eigen.valuesL) = length(eigen.vectorsL) = k
  
  # Loop to compute for each class the required elements
  for(i in 1:k){
    XX <- X[class == sort(unique(class))[i],]
    XX <- apply(as.matrix(XX), 2, as.numeric)
    
    # Barycentre
    mu <- colMeans(XX, na.rm = TRUE)
    
    
    # Marix of variance-covariance
    if(covEstim == 'empiric'){
      WkM <- cov(XX, use = 'complete.obs')
    }else if(covEstim == 'shrink'){
      WkM <- corpcor::cov.shrink(XX, verbose = FALSE)
    }else{
      if(nrow(XX) < ncol(XX)){
        WkM <- corpcor::cov.shrink(XX, verbose = FALSE)
      }else{
        WkM <- cov(XX, use = 'complete.obs')
      }
    }
    
    # Eigen values
    eigen.values <- eigen(WkM)$values
    
    # Eigen vectors
    eigen.vectors <- eigen(WkM)$vectors
    
    # Update the lists
    muL[[i]] <- mu
    eigen.valuesL[[i]] <- eigen.values
    eigen.vectorsL[[i]] <- eigen.vectors
    WkL[[i]] <- WkM
    
  }
  
  # Computation of the gamma-metric value
  gamma <- 0
  for(k1 in 1:(k-1)){
    for(k2 in (k1 + 1):k){
      gamma <- gamma + d_k1k2(muL[[k1]], muL[[k2]], eigen.valuesL[[k1]], eigen.valuesL[[k2]], eigen.vectorsL[[k1]], eigen.vectorsL[[k2]])
    }
  }
  
  # Plot
  if(plot == TRUE){
    id <- sample(1:nrow(dat1), size = 200)
    plot(X[id,1], X[id,2], xlab = colnames(X)[1], ylab = colnames(X)[2], col = class, main = paste('gamma =', round(gamma, 2)))
    for (i in 1:k){
      lines(x = c(muL[[i]][1], muL[[i]][1] + (eigen.vectorsL[[i]][1, 1]*sqrt(eigen.valuesL[[i]][1]))),
            y = c(muL[[i]][2], muL[[i]][2] + (eigen.vectorsL[[i]][2, 1]*sqrt(eigen.valuesL[[i]][1]))), col = 'blue')
      lines(x = c(muL[[i]][1], muL[[i]][1] + (eigen.vectorsL[[i]][1, 2]*sqrt(eigen.valuesL[[i]][2]))),
            y = c(muL[[i]][2], muL[[i]][2] + (eigen.vectorsL[[i]][2, 2]*sqrt(eigen.valuesL[[i]][2]))), col = 'blue')
      ellipse_k <- ellipse(WkL[[i]], muL[[i]])
      points(ellipse_k, t = 'l', col = 'green')
      
      # for(k2 in (i+1):k){
      #   if(k2 > k)break
      #   lines(x = c(muL[[i]][1], muL[[k2]][1]), y = c(muL[[i]][2], muL[[k2]][2]))
      # }
    }
  }
  
  return(gamma)
}


