##############################################################################################
#                                 GENERATION OF THE DATA                                     #
##############################################################################################

## Function to generate the data -------------------------------------------------------------
# This is the function used in scenario 1 and 2 to generate the data 
generation <- function(n, n_var_inf, n_noise, beta){
  ## Inputs:
  ##     n: Number of observations to generate
  ##     n_var_inf: Number of informative features
  ##     n_noise: Number of non-informative features
  ##     beta: Effect coefficients associated with the features.
  ## Ouput: The dataframe with the class y, the features named with X and an index for informative
  ##      features and Noise and an index for non-informative features
  
  # Draw normal distribution for the features
  X <- cbind(rep(1, n),                                                               # First column of 1 for the intercept 
             matrix(rnorm(n*(n_var_inf + n_noise), mean = 0, sd = 1), 
                    ncol = (n_var_inf + n_noise)))
  colnames(X) <- c('intercept', paste0('x', 1:n_var_inf), paste0('Noise', 1:n_noise))
  
  # Compute the probability rho_i
  X <- subset(X, select = c(intercept, grep('x', colnames(X)), grep('Noise', colnames(X))))
  rho_i <- exp(X %*% beta)/(1 + exp(X %*% beta))
  
  # Draw observation class
  Y <- rbinom(n = n, size = 1, prob = rho_i)
  
  # Checking if generated non-informative faetures are not too much correlated with the class Y (p.value < 0.1)
  inf_noise <- 0
  i <- 0
  while(length(inf_noise) != 0){
    # With a student test we check if noise features are significantly different from class Y = 0 and Y = 1
    p.value_noise <- sapply(colnames(X)[grep('Noise', colnames(X))], FUN = function(nom_var){
      var.p.value <- var.test(X[which(Y == 1), nom_var], X[which(Y == 0), nom_var])$p.value
      var.equal <- ifelse(var.p.value < 0.05, FALSE, TRUE)
      signif <- t.test(X[which(Y == 1), nom_var], X[which(Y == 0), nom_var], var.equal = var.equal)$p.value
      return(signif)
    })
    
    # We use a risk alpha = 10%
    inf_noise <- names(which(p.value_noise < 0.1))
    
    # We drop and draw the noisy features that have p.value < 0.1
    X[, inf_noise] <- matrix(rnorm(n = n*length(inf_noise), mean = 0, sd = 1), ncol = length(inf_noise), byrow = TRUE)
    X <- subset(X, select = c(intercept, grep('x', colnames(X)), grep('Noise', colnames(X))))
    i <- i+1
  }
  
  data <- data.frame(y = Y, X[, !colnames(X) %in% 'intercept'])
  data$y <- as.factor(data$y)
  
  return(data)
}


## Functions to generate the data in scenario 3 ----------------------------------------------

# This is the function to generate the covariance matrix used as a parameter of the multi-
# -variate Gaussian distribution.
produceRho <- function(s_g, n_g, alpha_max, c = 0.5, n_ind, constant_cor){
  ## Inputs:
  ##   s_g: The size of each group of features (same for every group)
  ##   n_g: The number of groups of features
  ##   alpha_max: The maximum level of correlation
  ##   c: The minimum level of correlation in the case where the correlation is non-constant
  ##   n_ind: Number  of independent group to consider
  ##   constant_cor: Boolean, TRUE for a constant correlation, FALSE for a non-constant correlation
  ## Output: 
  ##   Sigma: A block diagonal matrix of dimension s_g*n_g, symmetric and positive definite. 
  
  # Size of Sigma
  p <- s_g*n_g
  
  # Group of variables independent of each-others and with other group of variables
  independence_group_start <- n_g - n_ind
  
  # Initialize Sigma
  Sigma <- diag(p)
  
  # Define Sigma_g the block matrix (covariance matrix of a group of features)
  Sigma_g <- matrix(alpha_max, s_g, s_g)
  
  # Compute decay weight (in case of non-constant correlation)
  w <- ifelse(constant_cor, 0, -log(c/alpha_max)/(s_g-2)*0.99)
  
  # Compute Sigma_g with decay to correlations
  for(i in 1:s_g){
    Sigma_g[i, ] <- Sigma_g[i, ]*exp(-w*(abs(1:s_g - i) - 1))
  }
  
  # Variance of the features
  diag(Sigma_g) <- 1
  
  for(i in 1:(n_g)){
    ind_x <- 1 + (i - 1)*s_g
    ind_y <- 1 + (i - 1)*s_g
    
    if(i <= independence_group_start){
      Sigma[ind_x:(ind_x + s_g - 1), ind_y:(ind_y + s_g - 1)] <- Sigma_g
    }else if(i > independence_group_start){
      Sigma[ind_x:(ind_x + s_g - 1), ind_y:(ind_y + s_g - 1)] <- diag(s_g)
    }
  }
  
  # Rho is a symmetrical matrix
  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  
  return(Sigma)
}

# This is the function used to generate the beta vector
produceBeta <- function(s_g, n_g, non_zero_coeff, non_zero_group, value){
  ## Inputs:
  ##     s_g: The size of each group
  ##     n_g: The number of groups
  ##     non_zero_coeff: The index of the informative features in each group (the same for all groups)
  ##     non_zero_group: The index of the group that contains at least 1 informative features
  ##     value: The value of the beta coefficients associated with the informative features (the same for all informative features)
  ## Outputs:
  ##     The beta vector, including beta0
  
  # Initialize beta vector with 0
  beta_vec <- rep(0, s_g*n_g)
  
  # Indexes of non-null coefficients
  non_zero_indexes <- as.vector(sapply(s_g*(non_zero_group - 1), FUN = function(index){
    c(index + non_zero_coeff)
  }))
  
  # Replace 0 by value at the indexes of non-null coefficients
  beta_vec <- replace(beta_vec, non_zero_indexes, value)
  
  # We want to have balanced classes 
  # For standard normal distribution centered around zero, it should be 0
  # Balanced => P(Y = 1|X) = 1/2 <=> beta0 = - sum_{j = 1}^{p} E[X_ij]Beta_j
  # In the case where, for all j \in {1,...,p} X_j ~ N(0, 1) then beta0 = 0
  # In the case where we have a different distribution for X, it will be equals to the sum of non-zero beta coefficients times E[X_j].
  
  # For our work we did not considered other distributions
  beta0 <- 0
  
  return(c(beta0, beta_vec))
  
}

# This is the function to generate the data
genData <- function(n, s_g, n_g, beta0, beta, alpha_max, c, n_ind, constant_cor){
  ## Inputs: 
  ##    n: Number of observations
  ##    s_g: Size of each group (same for all groups)
  ##    n_g: Number of groups
  ##    beta0: The value of the intercept (beta0)
  ##    beta: The vector of beta coefficients
  ##    alpha_max: The maximum level of correlation
  ##    c: The minimum level of correlation (in case of non-constant correlation)
  ##    n_ind: The number of independent groups (independent groups are always from the last column of the dataset)
  ##                If we have a single independent group, it will be the last one. If we have 2 it will be the last
  ##                two ones. 
  ##    constant_cor: Boolean value, TRUE for constant correlation within each group, FALSE for non-constant correlation
  ##                within each group.
  ## Output: 
  ##    data: A data.frame with the class y being the first column and drawn observations X after. No intercept column in X.
  
  # Total number of columns
  p <- s_g*n_g
  
  # Produce the correlation matrix
  Sigma <- produceRho(s_g = s_g, n_g = n_g, alpha_max = alpha_max, n_ind = n_ind, constant_cor = constant_cor, c = c)
  
  # Check the definite positive
  if(sum(eigen(Sigma)$values < 0) > 0){
    cat('Specified correlation matrix is not positive definite')
  }else{
    # Initialize X
    X <- matrix(0, ncol = p, nrow = n)
    
    # Generate n instances of p correlated normal variables
    Z <- MASS::mvrnorm(n = n, rep(0, p), Sigma)
    
    # Transform Z to U using cumulative density function
    U <- pnorm(Z)
    
    # Part to modify if we want to include another distribution
    X <- qnorm(U)
  }
  
  # Compute linear predictor
  eta <- beta0 + X %*% beta
  
  # Compute probability of Y = 1
  rho_i <- exp(eta)/(1 + exp(eta))
  
  # Sample, from Bernouilli distribution, the class of i with probability rho_i
  Y <- rbinom(n, 1, p = rho_i)
  
  data <- data.frame(y = Y, X)
  data$y <- as.factor(data$y)
  
  return(data)
}
