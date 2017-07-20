library('rstan')
library('deSolve')

setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan")

NUM_VAR = 5
rs <- c(.1,.15,.2,.25,.3) #r_i's
ks <- c(200,200,200,200,200) #k_i's
as <- as.matrix(read.csv("stan_test.csv", header=FALSE)) #imports matrix of alpha_ij's
Ns <- c(1,1,1,1,1) #initial states
Time <- 100
deltaT <- 1
nstep <- Time/deltaT
time <- seq(deltaT,Time,deltaT) #time vector length and granularity

params <- data.frame(rs,ks,as) #each row contains set of params for ODE
states <- Ns

multivar_log_comp <- function(time,state,params) {
  odes <- rep(NA,NUM_VAR)
  for (i in 1:NUM_VAR) {
    comp_sum <- 0
    for (j in 1:NUM_VAR) { #this is actually faster than an inner product
      comp_sum <- comp_sum + params[i,2+j]*state[j] #because of class forcing
    }
    Ni <- params[i,1]*state[i]*(1-(comp_sum/params[i,2]))
    odes[i] <- Ni
  }
  return(list(odes))
}

out <- ode(y=states, times=time, func=multivar_log_comp, parms=params)
df <- as.data.frame(out)
nums <- as.matrix(df[,-1])

estimates <- stan(file = '5log_test.stan',
                  data = list (
                    N = array(c(NUM_VAR),dim=1),
                    Time  = nstep,
                    B0 = Ns,
                    K = ks,
                    z  = nums,
                    t0 = 0,
                    ts = time
                  ),
                  seed = 42,
                  chains = 1,
                  iter = 1000,
                  warmup = 500,
                  refresh = 100,
                  control = list(adapt_delta = 0.8)
)

print(estimates)