library(MASS) # for mvrnorm


l = 100      # length of time series
KK = 300     # carrying capacity ("K" is for kernel below)
r = 1.6      # population growth rate

# Simulate time series
n = integer(l)
n[1] = KK
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
# plot the time series




# Gaussian process regression ---------------------------------------------

# Code based on Carl Boettiger's Cholesky method
# https://github.com/cboettig/nonparametric-bayes/blob/master/R/gp_fit.R


# Add (0,0) as a "data" point: extinction is forever
x = c(0, d[[1]])
y = c(0, d[[2]]) - mean(d[[2]])

# Range of values to plot
x_seq = seq(0, max(d) * 1.2, length = 1000)

process_noise = 100  # populations can fluctuate this much for non-density reasons
lengthscale = max(d) * 5   # Very smooth function: values should be similar across a broad range
sigma_f = 1E9        # very flat prior: function can go anywhere.

## Squared exponential kernel
SE <- function(Xi,Xj, l=lengthscale) sigma_f * exp(-0.5 * (Xi - Xj)^2 / lengthscale^2)
covar <- function(X, Y) outer(X, Y, SE, lengthscale) 
K <- covar(x, x)
I <-  diag(1, length(x))
I[1, 1] = 1E-5 # No process noise for fact that curve goes though (0,0)
K_altered = K + process_noise ^ 2 * I

# Cholesky decomposition method
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

# Plot 2
plot(
  d, 
  xlim = c(0, max(x_seq)), 
  ylim = c(0, max(d)) * 1.04, 
  pch = 1, 
  xaxs = "i", 
  yaxs = "i", 
  col = "#0000C0",
  lwd = 2,
  cex = 2,
  main = "Density dependence:\n 250 curves from posterior + \"truth\""
)
# Spaghetti lines showing possible density dependence curves
matplot(
  x_seq, 
  sample_curves, 
  type = "l",
  lty = 1,
  add = TRUE,
  col = "#00000040",
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
slopes = (sample_curves[851, ] - sample_curves[850, ]) / (x_seq[851] - x_seq[850])
plot(
  density(slopes),
  xlab = paste("slope at at N =", round(mean(x_seq[850:851]))),
  main = paste("Posterior for slope at N =", round(mean(x_seq[850:851])))
)
abline(v = 0, col = 2, lwd = 2)
