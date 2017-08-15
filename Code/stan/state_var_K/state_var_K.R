library('rstan')
library('deSolve')
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping

setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan/state_var_K")

DUMMY_DATA <- TRUE
PLOT <- FALSE

cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7") #colorblind friendly colors

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

if (DUMMY_DATA){
  Ags <- 66 #nM
  AC50 <- 2
  b_i <- 0.1
  M <- 13000
  tau <- 0.01
  t_peak <- 21 #days
  omega <- -3
  
  Time <- 42
  deltaT <- 7
  t0 <- 0
  nstep <- Time/deltaT
  time <- seq(deltaT,Time,deltaT)
  
  state_vals <- c(B=1,K=M*b_i)
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
      dB <- ((Ags*tau)/AC50)*B*(1 - B/K)
      list(c(dB,dK))
    })
  }
  dummy <- ode(y = state_vals, times = time, func = model_func, parms = params)
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

estimates <- stan(file = 'state_var_K.stan',
                  data = list (
                    n  = nstep,
                    B0 = state_vals[1],
                    z  = dummy[,-1],
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
                  sample_file = 'state_var_K_samples.csv',
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 10)
)

print(estimates)
