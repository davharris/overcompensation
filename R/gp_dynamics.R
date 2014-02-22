# sigma_n: populations can fluctuate this much for non-density reasons
# lengthscale: determines the smoothness of the learned functions
# sigma_f: determines how far curves stray from the grand mean.
gp_dynamics = function(
  N_t, 
  N_t_plus_1, 
  process_noise = sd(c(0, N_t_plus_1)),
  lengthscale = sd(N_t) ,
  sigma_f = sd(c(0, N_t_plus_1)),
  test_values = seq(0, max(N_t, N_t_plus_1) * 1.2, length = 1000)[-1],
  predict = TRUE
){
  
  # Code based on Carl Boettiger's Cholesky method
  # https://github.com/cboettig/nonparametric-bayes/blob/master/R/gp_fit.R
  
  # Notation & algorithm mostly follow Rasmussen and Williams
  # http://www.gaussianprocess.org/gpml/chapters/
  # See especially Chapter 2 and Algorithm 2.1
  
  # Pull out the data from the data frame
  # Add (0,0) as a "data" point: extinction is forever
  x = c(0, N_t)
  y = c(0, N_t_plus_1) - mean(N_t_plus_1) # y is scaled to be mean zero
    
  ## Squared exponential kernel
  SE <- function(Xi,Xj, l=lengthscale){
    sigma_f^2 * exp(-0.5 * (Xi - Xj)^2 / lengthscale^2)
  }
  covar <- function(X, Y) outer(X, Y, SE, lengthscale) 
  K <- covar(x, x)
  
  # Add process noise along kernel's diagonal
  I <-  diag(1, length(x))
  I[1, 1] = 1E-5 # No process noise for the *fact* that curve goes though (0,0)
  K_altered = K + process_noise ^ 2 * I
  
  # Cholesky decomposition method for GP regression
  L <- t(chol(K_altered))
  alpha <- solve(t(L), solve(L, y))
  
  if(predict){
    
    k_star <- covar(x, test_values)
    
    # Expected values along test_values; note that the mean is added back in.
    mu <- t(k_star) %*% alpha + mean(N_t_plus_1) 
    v <- solve(L, k_star)
    
    # Correlated uncertainty along test_values
    var <- covar(test_values,test_values) - t(v) %*% v
  }
  
  
  
  # Can use loglik to optimize process_noise, lengthscale, and sigma_f
  loglik = -.5 * t(y) %*% alpha - sum(log(diag(L))) - length(x) * log(2 * pi) / 2
  
  list(
    K = K_altered,
    L = L,
    mu = if(predict){mu},
    var = if(predict){var},
    loglik = loglik,
    test_values = if(predict){test_values},
    covar = covar,
    v = if(predict){v}
  )
}



