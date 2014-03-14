library(corpcor) # for make.positive.definite
library(rethinking)

devtools::load_all()
source("inst/fake_logistic.R")

iter = 10000
logliks = numeric(iter)


proposal = function(state, scale = c(.1, .5, 1)){
  state + rcauchy(length(state), location = 0, scale = scale)
}


# Initialize the chain
state = c(
  log_process_noise = 4, 
  log_lengthscale = 6, 
  log_sigma_f = 6
)
llik = gp_dynamics(
  d[[1]], 
  d[[2]], 
  process_noise = exp(state["log_process_noise"]),
  lengthscale = exp(state["log_lengthscale"]),
  sigma_f = exp(state["log_sigma_f"]),
  predict = FALSE
)$loglik
post = llik

states = matrix(NA, nrow = iter, ncol = length(state))
colnames(states) = names(state)


for(i in 1:iter){
  # Make a proposal
  proposed_state = proposal(state)
  gp = tryCatch(
    gp_dynamics(
      d[[1]], 
      d[[2]], 
      process_noise = exp(proposed_state["log_process_noise"]),
      lengthscale = exp(proposed_state["log_lengthscale"]),
      sigma_f = exp(proposed_state["log_sigma_f"]),
      predict = FALSE
    ),
    error = function(x){
      # message("Metropolis proposal rejected: error in GP")
      gp = gp
      gp$loglik = -Inf
      gp
    }
    
  )
  
  proposed_llik = gp$loglik
  proposed_post = proposed_llik
  
  
  # Metropolis transition rule
  alpha = exp(proposed_post - post)
  if(alpha > runif(1)){
    state = proposed_state
    post = proposed_post
    logliks[i] = gp$loglik
  }else{
    if(i==1){
      logliks[i] = gp$loglik
    }else{
      logliks[i] = logliks[i - 1]
    }
  }
  states[i, ] = state
  
  if(i %% 100 ==0){cat(i, "\n")}
}

matplot(states, type = "l", lty = 1, col = c(1, 2, 4))

library(coda)
effectiveSize(states)

test_values = seq(0, 800, length = 100)

thin = 20
curves = sapply(
  (1:iter)[(1:iter) %% thin == 0],
  function(i){
    gp = gp_dynamics(
      d[[1]], 
      d[[2]], 
      process_noise = exp(states[i, "log_process_noise"]),
      lengthscale = exp(states[i, "log_lengthscale"]),
      sigma_f = exp(states[i, "log_sigma_f"]),
      predict = TRUE,
      test_values = test_values
    )
    with(
      gp,{
        Sigma = make.positive.definite(
          covar(test_values,test_values) - t(v) %*% v
        )
        mvrnorm(1, mu, Sigma)
      }
    )
  }
)

matplot(
  test_values,
  curves, 
  type = "l",
  lty = 1,
  col = "#00000020", 
  ylim = c(0, 700),
  xaxs = "i",
  yaxs = "i",
  lwd = 1
)
with(
  list(
    KK = 500,
    r = 1.8
  ),
  curve(
    x + r * x * (1 - x / KK), 
    add = TRUE, 
    col = 2, 
    lwd = 4, 
    from = 0, 
    to = max(test_values)
  )
)
points(d, col = "lightgray", cex = 1, lwd = 2, pch = 21, bg = "slateblue")
test_point = which.min(abs(test_values - max(d)))



test_point = which.min(abs(test_values - max(d)))
rises = curves[test_point, ] - curves[(test_point - 1), ]
run = test_values[test_point] - test_values[test_point - 1]
slopes = rises / run
abline(v = test_values[c(test_point-1, test_point)], lty = 2)

hist(slopes)
