library('rstan')
library('deSolve')
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping


setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan/model1")

DUMMY_DATA <- TRUE
PLOT <- TRUE

cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7") #colorblind friendly colors

if (DUMMY_DATA){
  Ags <- 66 #nM
  AC50 <- 2
  b_i <- 0.05
  b_f <- 0.005
  M <- 13000
  gamma <- 0.1
  tau <- 0.5
  t_peak <- 14 #days
  mu <- 0.1
  
  Time <- 40
  deltaT <- 1
  nstep <- Time/deltaT
  time <- seq(deltaT,Time,deltaT)
  
  state_vals <- c(B=1,A=0)
  params <- c(Ags=Ags,
              AC50=AC50,
              b_i=b_i,
              b_f=b_f,
              M=M,
              gamma=gamma,
              tau=tau,
              t_peak=t_peak,
              mu=mu)
  time <- seq(deltaT,Time,deltaT)
  
  model1dummy <- function(t, state, parms){
    with(as.list(c(state, parms)),{
      if (t < t_peak) {
        dB <- ((Ags*tau)/AC50)*(1 - state[1]/(M*b_i))
      } else {
        dB <- ((Ags*tau)/AC50)*(1 - state[1]/(M*b_f))
      }
      dA <- gamma*state[1] - mu*state[2]
      list(c(dB,dA))
    })
  }
  dummy <- ode(y = state_vals, times = time, func = model1dummy, parms = params)
}

if (PLOT) {
  df <- as.data.frame(dummy)
  m_df <- melt(df, id=c('time'))
  
  p <- ggplot(m_df, aes(x = time, y = value, color = variable))
  plot <- p + geom_line() + ylab("y") + ggtitle('Single Epitope dummy data')
  plot <- plot + scale_color_manual(values=cbbPalette)
  print(plot)
}
