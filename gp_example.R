library(devtools)
load_all()
source("inst/fake_logistic.R") # simulate fake data

# Gaussian process regression ---------------------------------------------


# Optimize hyperparameters:
#   * sigma_n: populations can fluctuate this much for non-density reasons
#   * lengthscale: determines the smoothness of the learned functions
#   * sigma_f: determines how far curves stray from the grand mean.
# The real project should sample from the whole posterior distribution of these
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
sample_curves = with(
  gp,
  t(mvrnorm(250, mu, covar(test_values, test_values) - t(v) %*% v))
)



# plots -------------------------------------------------------------------
par(mfrow = c(1, 3))

# Plot 1
plot(n, type = "o", xlab = "time", cex = .5, pch = 16, main = "Time series")
abline(h = KK, col = "#00000050", lty = 2, lwd = 2)

# Plot 2
plot(
  d, 
  xlim = c(0, max(gp$test_values)), 
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
  gp$test_values, 
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
  to = max(gp$test_values)
)
test_point = which.min(abs(gp$test_values - max(d)))
abline(v = gp$test_values[c(test_point, test_point - 1)], lty = 2)


# Plot 3
rises = sample_curves[test_point, ] - sample_curves[(test_point - 1), ]
run = gp$test_values[test_point] - gp$test_values[test_point - 1]
slopes = rises / run

plot(
  density(slopes, to = 1),
  xlab = paste("slope at at N =", round(mean(gp$test_values[test_point]))),
  main = paste(
    "Posterior for slope at N =", 
    round(mean(gp$test_values[c(test_point - 1, test_point)])),
    "\n (near highest observed value)"
  ), 
  yaxs = "i"
)
abline(v = 0, col = 2, lwd = 2)
rug(slopes)
