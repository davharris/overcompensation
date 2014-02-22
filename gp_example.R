library(devtools)
load_all()
source("inst/fake_logistic.R") # simulate fake data

# Gaussian process regression ---------------------------------------------

gp = gp_dynamics(
  N_t = d[[1]],
  N_t_plus_1 = d[[2]]
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

# Plot 3
slopes = (sample_curves[831, ] - sample_curves[830, ]) / (gp$test_values[831] - gp$test_values[830])
plot(
  density(slopes),
  xlab = paste("slope at at N =", round(mean(gp$test_values[830:831]))),
  main = paste(
    "Posterior for slope at N =", 
    round(mean(gp$test_values[830:831])),
    "\n (near highest observed value)"
  ), 
  yaxs = "i"
)
abline(v = 0, col = 2, lwd = 2)
rug(slopes)
