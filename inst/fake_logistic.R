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

set.seed(NULL)
