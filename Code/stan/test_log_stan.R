#Estimates r for a single variable logistic growth with good data
#With dense but more noisey data, still does well

library('rstan')
library('deSolve')

setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan")

NOISE <- TRUE

r <- 0.1
k <- 100
N <- 1
Time <- 100
deltaT <- 0.1
nstep <- Time/deltaT
t0 <- 0
STD <- 5

simple_model_params <- c(r = r, k = k)
simple_state_val <- c(N = N)
time <- seq(deltaT,Time,deltaT)

simple_log_growth <- function(t, state, params) {
  with(as.list(c(state, params)),{
    dN <- r*N*(1-(N/k))
    list(c(dN))
  })
}

out <- ode(y = simple_state_val, times = time, func = simple_log_growth, parms = simple_model_params)

if (NOISE){
  nums <- out[,'N']
  for (i in 1:length(nums)){
    nums[i] <- nums[i] + rnorm(1, sd=STD)
    if (nums[i] < 0) {
      nums[i] <- 0
    }
  }
} else {
  nums <- out[,'N']
}

y0 <- N

estimates <- stan(file = 'test_log.stan',
                  data = list (
                    Time  = nstep,
                    B0 = array(c(N), dim=1),
                    K = array(c(k), dim=1),
                    z  = nums,
                    t0 = t0,
                    ts = time
                  ),
                  seed = 42,
                  chains = 1,
                  iter = 1000,
                  warmup = 500,
                  refresh = -1
)
#Decreasing warmup to 200 increases sampling time by like 15 minutes
