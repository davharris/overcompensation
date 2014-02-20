set.seed(1)

library(MASS) # for mvrnorm


l = 50      # length of time series
KK = 500     # carrying capacity ("K" is for kernel below)
r = 1.8      # population growth rate

# Simulate time series
n = integer(l)
n[1] = KK # Initialize the population at carrying capacity
for(i in 1:(l - 1)){
  expected = n[i] + r * n[i] * (1 - n[i] / KK)
  expected.plus.noise = if(expected == 0){
    0
  }else{
    max(0, rnorm(1, expected, 50))
  }
  n[i + 1] = rpois(1, expected.plus.noise)
}

d = data.frame(n_t = n[-l], n_t_plus_one = n[-1])


# Gaussian process regression ---------------------------------------------

# Code based on Carl Boettiger's Cholesky method
# https://github.com/cboettig/nonparametric-bayes/blob/master/R/gp_fit.R

# Notation & algorithm mostly follow Rasmussen and Williams
# http://www.gaussianprocess.org/gpml/chapters/
# See especially Chapter 2 and Algorithm 2.1

# Pull out the data from the data frame
# Add (0,0) as a "data" point: extinction is forever
x = c(0, d[[1]])
y = c(0, d[[2]]) - mean(d[[2]])

# Tuning parameters for the GP
process_noise = 100  # sigma_n: populations can fluctuate this much for non-density reasons
lengthscale = sd(y) * 10   # smoothness of the learned functions
sigma_f = 1E5        # Function variance is very high: function can go anywhere eventually.

# Range of values to plot
x_seq = seq(0, max(d) * 1.2, length = 1000)

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
k_star <- covar(x, x_seq)


mu <- t(k_star) %*% alpha + mean(d[[2]]) # Expected values along x_seq
v <- solve(L, k_star)
var <- covar(x_seq,x_seq) - t(v) %*% v # Correlated uncertainty along x_seq


# Can use loglik to optimize process_noise, lengthscale, and sigma_f
loglik <- -.5 * t(y) %*% alpha - sum(log(diag(L))) - length(n) * log(2 * pi) / 2


# Examples of possible density dependence relationships from model's posterior
sample_curves = t(mvrnorm(250, mu, covar(x_seq,x_seq) - t(v) %*% v))



# plots -------------------------------------------------------------------
par(mfrow = c(1, 3))

# Plot 1
plot(n, type = "o", xlab = "time", cex = .5, pch = 16, main = "Time series")
abline(h = KK, col = "#00000050", lty = 2, lwd = 2)

# Plot 2
plot(
  d, 
  xlim = c(0, max(x_seq)), 
  ylim = c(0, max(d)) * 1.04, 
  pch = 1, 
  xaxs = "i", 
  yaxs = "i", 
  col = "#0000C0",
  lwd = 3,
  cex = 2,
  main = "Density dependence:\n 250 curves from posterior + \"truth\" in red"
)
# Spaghetti lines showing possible density dependence curves
matplot(
  x_seq, 
  sample_curves, 
  type = "l",
  lty = 1,
  add = TRUE,
  col = "#00000030",
  lwd = 1
)
# "True" relationship in red
curve(
  x + r * x * (1 - x / KK), 
  add = TRUE, 
  col = 2, 
  lwd = 4, 
  from = 0, 
  to = max(x_seq)
)

# Plot 3
slopes = (sample_curves[831, ] - sample_curves[830, ]) / (x_seq[831] - x_seq[830])
plot(
  density(slopes),
  xlab = paste("slope at at N =", round(mean(x_seq[830:831]))),
  main = paste(
    "Posterior for slope at N =", 
    round(mean(x_seq[830:831])),
    "\n (near highest observed value)"
  )
)
abline(v = 0, col = 2, lwd = 2)
