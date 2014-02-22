library(devtools)
load_all()
source("inst/fake_logistic.R") # simulate fake data

pdf("example.pdf", width = 6, height = 2, pointsize = 6)

# Gaussian process regression ---------------------------------------------


# Optimize hyperparameters:
#   * sigma_n: populations can fluctuate this much for non-density reasons
#   * lengthscale: determines the smoothness of the learned functions
#   * sigma_f: determines how far curves stray from the grand mean.
# The real project should sample from the whole posterior distribution of these
#    rather than relying on point estimates
opt = optim(
  par = c(
    process_noise = sd(c(0, d[[2]])),                
    lengthscale = sd(d[[1]]), 
    sigma_f = sd(c(0, d[[2]]))
  ),
  fn = function(par){
    gp = gp_dynamics(
      N_t = d[[1]],
      N_t_plus_1 = d[[2]],
      process_noise = par["process_noise"],
      lengthscale = par["lengthscale"],
      sigma_f = par["sigma_f"],
      predict = FALSE
    )
    -gp$loglik
  }
)


# Fit a GP with the optimal hyperparameters
gp = with(
  as.list(opt$par),
  gp_dynamics(
    N_t = d[[1]],
    N_t_plus_1 = d[[2]],
    process_noise = process_noise,
    lengthscale = lengthscale,
    sigma_f = sigma_f,
    predict = TRUE,
    test_values = seq(0, 1001, length = 1000)
  )
)


# Examples of possible density dependence relationships from model's posterior
n_curves = 500

sample_curves = with(
  gp,
  t(mvrnorm(n_curves, mu, covar(test_values, test_values) - t(v) %*% v))
)



# plots -------------------------------------------------------------------
par(mfrow = c(1, 3))
par(cex.lab = 1.5)
par(cex.axis = 1.5)
par(cex.main = 1.5)

# Plot 1: raw time series
plot(
  n, 
  type = "o", 
  xlab = "time", 
  cex = 1, 
  pch = 16, 
  main = "Time series", 
  xaxs = "i",
  ylab = "N"
)
abline(h = KK, col = "#00000050", lty = 2, lwd = 1)

# Plot 2: Spaghetti lines showing possible density dependence curves
matplot(
  gp$test_values, 
  sample_curves, 
  type = "l",
  lty = 1,
  col = "#00000020",
  lwd = 1/2,
  xlim = c(0, max(gp$test_values)), 
  ylim = c(0, max(d)) * 1.04, 
  xaxs = "i", 
  yaxs = "i", 
  main = paste0(
    "Density dependence:\n",
    n_curves, 
    " curves from posterior + \"truth\" in red"
  ),
  xlab = "N_t",
  ylab = "N_t+1"
)

# Test the slope near the maximum observed value
test_point = which.min(abs(gp$test_values - max(d)))
abline(v = gp$test_values[c(test_point, test_point - 1)], lty = 2)

points(
  d,
  col = "#2020D0",
  bg = "lightgray",
  cex = 1.5,
  pch = 21
)

# "True" relationship in red
curve(
  x + r * x * (1 - x / KK), 
  add = TRUE, 
  col = 2, 
  lwd = 1.5, 
  from = 0, 
  to = max(gp$test_values)
)


# Plot 3: posterior density of slopes
rises = sample_curves[test_point, ] - sample_curves[(test_point - 1), ]
run = gp$test_values[test_point] - gp$test_values[test_point - 1]
slopes = rises / run

plot(
  density(slopes, to = 1, n = 1000),
  xlab = paste("slope at at N =", round(mean(gp$test_values[test_point]))),
  main = paste(
    "Posterior for slope at N =", 
    round(mean(gp$test_values[c(test_point - 1, test_point)])),
    "\n (near highest observed value)"
  ), 
  yaxs = "i"
)
abline(v = 0, col = 2)
rug(slopes, col = "#00000040")

dev.off()
