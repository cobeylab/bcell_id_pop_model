library('rstan')
library('deSolve')
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping

setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan/single_beta")

DUMMY_DATA <- 1
NOISE <- 1
PLOT <- 0
STAN <- 1
FIT_PLOT <- 0
STD <- 1

cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7") #colorblind friendly colors

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if (DUMMY_DATA) {
  taus <- c(.7,.7)
  c_50s <- c(2,1.8)
  rho <- 7000
  beta <- 0.25
  omega <- -3
  Bs <- c(1590,1412)
  t_peak <- 7 
  TIME <- 28
  deltaT <- 1
  t0 <- 0
  nstep <- TIME/deltaT
  time <- seq(deltaT,TIME,deltaT)
  
  tv_params <- c(taus=taus,c_50s=c_50s,rho=rho,beta=beta,omega=omega,t_peak=t_peak)
  tv_state <- c(B1=Bs[1],B2=Bs[2],K1=rho/c_50s[1],K2=rho/c_50s[2])
  
  tv_log_compete <- function(time,state,params) {
    with(as.list(c(state, params)),{
      dB1 <- (taus[1]/c_50s[1])*B1*(1-((B1+beta*B2)/K1))
      dB2 <- (taus[2]/c_50s[2])*B2*(1-((B2+beta*B1)/K2))
      if (time <= t_peak) {
        dK1 <- 0
        dK2 <- 0
      } else {
        dK1 <- omega*(time-t_peak)
        dK2 <- omega*(time-t_peak)
      }
      list(c(dB1,dB2,dK1,dK2))
    })
  }
  out <- ode(y=tv_state, times=time, func=tv_log_compete, parms=tv_params)
}

if (PLOT) {
  df <- as.data.frame(out)
  m_df <- melt(df, id=c('time'))
  
  p <- ggplot(m_df, aes(x = time, y = value, color = variable))
  plot <- p + geom_line() + ylab("y") + ggtitle('Two Var. Growth w/ Comp.')
  plot <- plot + scale_color_manual(values=cbbPalette)
  print(plot)
}

data <- out[,-1]

if (STAN) {
  estimates <- stan(file = 'single_beta.stan',
                    data = list (
                      n  = nstep,
                      B0 = Bs,
                      z  = data,
                      t0 = t0,
                      ts = time,
                      AC50 = c_50s,
                      tau = taus,
                      t_peak = t_peak,
                      omega = omega
                    ),
                    chains = 4,
                    iter = 2000,
                    warmup = 1000,
                    refresh = 10,
                    control = list(adapt_delta = 0.8,
                                   max_treedepth = 6)
  )
  
  print(estimates)
}