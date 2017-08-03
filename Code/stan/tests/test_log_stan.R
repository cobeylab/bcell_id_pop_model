#Estimates r for a single variable logistic growth with good data
#With dense but more noisy data, still does well

library('rstan')
library('deSolve')

setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan/tests")

# These lines are for parallelization locally
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

NOISE <- TRUE
MULTI_SERIES <- TRUE

r <- 0.15
k <- 500
N <- 1
Time <- 100
deltaT <- 5
nstep <- Time/deltaT
t0 <- 0
STD <- 1

if (!MULTI_SERIES) {
  NUM_SERIES <- 1
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
    nums <- array(out[,'N'],dim=c(nstep,1))
    for (i in 1:length(nums)){
      nums[i] <- nums[i] + rnorm(1, sd=STD)
      if (nums[i] < 0) {
        nums[i] <- 0
      }
    }
  } else {
    nums <- array(out[,'N'],dim=c(nstep,1))
  }
  
  y0 <- N
} else if (MULTI_SERIES) {
  NUM_SERIES <- 3
  
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
  
  nums <- matrix(rep(out[,'N'], NUM_SERIES),ncol=NUM_SERIES,nrow=nstep)
  
  if (NOISE){
    for (i in 1:length(nums[,1])){
      for (j in 1:length(nums[1,])){
        nums[i,j] <- nums[i,j] + rnorm(1, sd=STD)
        if (nums[i,j] < 0) {
          nums[i,j] <- 0
        }
      }
    }
  }
}

estimates <- stan(file = 'test_log.stan',
                  data = list (
                    N = array(c(NUM_SERIES),dim=1),
                    L  = nstep,
                    B0 = array(c(N), dim=1),
                    K = array(c(k), dim=1),
                    z  = nums,
                    t0 = t0,
                    ts = time
                  ),
                  chains = 4,
                  iter = 2000,
                  warmup = 1000,
                  refresh = 100,
                  control = list(adapt_delta = 0.8)
)

print(estimates)
#Decreasing warmup to 200 increases sampling time by like 15 minutes
#Dramatically increasing data ruins predictive power
#Accurate with less steps also
