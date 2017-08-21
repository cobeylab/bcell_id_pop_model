library('rstan')
library('deSolve')
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping

setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan/state_var_K")

DUMMY_DATA <- 0
NOISE <- 0
PLOT <- 0
STAN <- 0
FIT_PLOT <- 1
STD <- 10

cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7") #colorblind friendly colors

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if (DUMMY_DATA){
  Ags <- 66 #nM
  AC50 <- 2
  b_i <- 1/AC50
  M <- 6000
  Tp <- 0.9
  t_peak <- 7#days
  omega <- -3
  tau <- 0.9
  
  Time <- 28
  deltaT <- 1
  t0 <- 0
  nstep <- Time/deltaT
  time <- seq(deltaT,Time,deltaT)
  
  state_vals <- c(B=1590,K=M*b_i)
  params <- c(Ags=Ags,
              AC50=AC50,
              b_i=b_i,
              M=M,
              tau=tau,
              t_peak=t_peak)
  
  model_func <- function(t, state, parms){
    with(as.list(c(state, parms)),{
      if (t <= t_peak) {
        dK <- 0
      } else {
        dK <- omega*(t-t_peak)
      }
      dB <- ((tau)/AC50)*B*(1 - B/K)
      list(c(dB,dK))
    })
  }
  dummy <- ode(y = state_vals, times = time, func = model_func, parms = params)
}

if (NOISE) {
  for (i in 1:length(dummy[,1])){
    for (j in 2:length(dummy[1,])){
      dummy[i,j] <- dummy[i,j] + rnorm(1,sd=STD)
      if (dummy[i,j] < 0) {
        dummy[i,j] <- 0
      }
    }
  }
}

if (PLOT) {
  df <- as.data.frame(dummy)
  m_df <- melt(df, id=c('time'))
  
  p <- ggplot(m_df, aes(x = time, y = value, color = variable))
  plot <- p + geom_line() + ylab("y") + ggtitle('Single Epitope dummy data')
  plot <- plot + scale_color_manual(values=cbbPalette)
  print(plot)
}

dummy_data <- dummy[,-1]

if (STAN) {
  estimates <- stan(file = 'state_var_K.stan',
                    data = list (
                      n  = nstep,
                      B0 = state_vals[1],
                      z  = dummy_data,
                      t0 = t0,
                      ts = time,
                      Ags = Ags,
                      AC50 = AC50,
                      bi = b_i,
                      tw = t_peak,
                      omega = omega
                    ),
                    chains = 4,
                    iter = 2000,
                    warmup = 1000,
                    refresh = 100,
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 10)
  )
  
  print(estimates)
}

if(FIT_PLOT){
  Time <- 28
  deltaT <- 1
  t0 <- 0
  nstep <- Time/deltaT
  time <- seq(deltaT,Time,deltaT)
  samples <- extract(estimates, c('theta[1]','theta[2]'))
  fit_M <- mean(samples[[1]])
  fit_tau <- mean(samples[[2]])
  state_vals <- c(fit_B=1590,fit_K=fit_M*b_i)
  params <- c(Ags=Ags,
              AC50=AC50,
              b_i=b_i,
              M=M,
              tau=fit_tau,
              t_peak=t_peak)
  
  fit_func <- function(t, state, parms){
    with(as.list(c(state, parms)),{
      if (t <= t_peak) {
        dK <- 0
      } else {
        dK <- omega*(t-t_peak)
      }
      dB <- ((tau)/AC50)*fit_B*(1 - fit_B/fit_K)
      list(c(dB,dK))
    })
  }
  
  fit <- ode(y = state_vals, times = time, func = fit_func, parms = params)
  
  fit_df <- as.data.frame(fit)
  fit_m_df <- melt(fit_df, id=c('time'))
  
  both <- rbind(m_df,fit_m_df)
  
  p <- ggplot(both, aes(x = time, y = value, color = variable))
  plot <- p + geom_line() + ylab("y") + ggtitle('Sa Fit for rho_m=6000 and tau=0.9')
  plot <- plot + scale_color_manual(values=cbbPalette)
  print(plot)
}
