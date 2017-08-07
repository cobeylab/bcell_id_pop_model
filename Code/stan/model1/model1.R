library('rstan')
library('deSolve')
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping


setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/stan/model1")

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
  b_i <- 0.05
  b_f <- 0.005
  M <- 13000
  gamma <- 0.1
  tau <- 0.01
  t_peak <- 21 #days
  mu <- 0.1
  
  Time <- 40
  deltaT <- 4
  t0 <- 0
  nstep <- Time/deltaT
  time <- seq(deltaT,Time,deltaT)
  
  state_vals <- c(B=1,A=0.001)
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
        dB <- ((Ags*tau)/AC50)*B*(1 - B/(M*b_i))
      } else {
        dB <- ((Ags*tau)/AC50)*B*(1 - B/(M*b_f))
      }
      dA <- gamma*B - mu*A
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

dummy_data <- dummy[,-1]

estimates <- stan(file = 'model1.stan',
                  data = list (
                    n  = nstep,
                    B0 = state_vals,
                    z  = dummy[,-1],
                    t0 = t0,
                    ts = time,
                    Ags = Ags,
                    AC50 = AC50,
                    bi = b_i,
                    bf = b_f
                  ),
                  chains = 4,
                  iter = 2000,
                  warmup = 1000,
                  refresh = 1,
                  sample_file = 'model1_samples.csv',
                  control = list(adapt_delta = 0.8)
)

print(estimates)
