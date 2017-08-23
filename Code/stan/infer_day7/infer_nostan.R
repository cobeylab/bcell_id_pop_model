library('rstan')
library('deSolve')
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping

setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan/infer_day7")

DUMMY_DATA <- 1
NOISE <- 1
PLOT <- 0
WALK <- 1
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
  TIME <- 35
  deltaT <- 1
  t0 <- 7
  nstep <- ((TIME - 14)/deltaT) + 1
  time <- seq(0,TIME,deltaT)
  
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

if (WALK) {
  taus <- c(.7,.7)
  c_50s <- c(2,1.8)
  rho <- 7000
  beta <- 0.25
  omega <- -3
  stepsize <- 0.1
  times <- seq(14,7,-stepsize)
  valB1 <- rep(NA,length(times))
  valB2 <- rep(NA,length(times))
  valK1 <- rep(NA,length(times))
  valK2 <- rep(NA,length(times))
  valB1[1] = 1590
  valB2[1] = 1412
  valK1[1] = rho/c_50s[1]
  valK2[1] = rho/c_50s[2]
  
  for (i in 1:(length(times)-1)) {
    d1 <- (taus[1]/c_50s[1])*valB1[i]*(1-((valB1[i]+beta*valB2[i])/valK1[i]))
    d2 <- (taus[2]/c_50s[2])*valB2[i]*(1-((valB2[i]+beta*valB1[i])/valK2[i]))
    d3 <- 0
    d4 <- 0
    valB1[i+1] <- valB1[i] - d1*stepsize
    valB2[i+1] <- valB2[i] - d2*stepsize
    valK1[i+1] <- valK1[i] - d3*stepsize
    valK2[i+1] <- valK2[i] - d4*stepsize
  }
  plot(times,valB1)
  lines(times,valB2)
}
