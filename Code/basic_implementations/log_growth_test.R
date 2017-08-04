#The point of this file is to get logistic growth implemented
#as well as get a real ODE solver up and running (deSolve)

#install.packages("deSolve")
#install.packages("ggplot2")
#install.packages("reshape2")

ptm <- proc.time() #timing

library("deSolve") #ode solver
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping

ONE_VAR <- TRUE
TWO_VAR <- FALSE
MULTI_VAR <- FALSE

cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7") #colorblind friendly colors

#setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/basic_implementations")

if (ONE_VAR == TRUE) {
  r <- 3.3
  k <- 650
  N <- 1
  
  simple_model_params <- c(r = r, k = k)
  simple_state_val <- c(N = N)
  time <- seq(0,100,0.1)
  
  simple_log_growth <- function(t, state, params) {
    with(as.list(c(state, params)),{
      dN <- r*N*(1-(N/k))
      list(c(dN))
    })
  }
  
  out <- ode(y = simple_state_val, times = time, func = simple_log_growth, parms = simple_model_params)
  
  header <- sprintf("One variable logisitc growth where r=%f & k=%f",r,k)
  plot <- qplot(out[,'time'],out[,'N'], xlab='time',ylab='number of things', geom='line',
        main = header)
  plot <- plot + theme(plot.title = element_text(size=8))
  print(plot)
}

if (TWO_VAR == TRUE) {
  rs <- c(.1,.2)
  ks <- c(100,100)
  as <- c(0.1,0.2)
  Ns <- c(1,1)
  time <- seq(0,100,0.01)
  
  tv_params <- c(r1=rs[1],r2=rs[2],k1=ks[1],k2=ks[2],a12=as[1],a21=as[2])
  tv_state <- c(N1=Ns[1],N2=Ns[2])
  
  tv_log_compete <- function(time,state,params) {
    with(as.list(c(state, params)),{
      dN1 <- r1*N1*(1-((N1+a12*N2)/k1))
      dN2 <- r2*N2*(1-((N2+a21*N1)/k2))
      list(c(dN1,dN2))
    })
  }
  out <- ode(y=tv_state, times=time, func=tv_log_compete, parms=tv_params)
  
  df <- as.data.frame(out)
  m_df <- melt(df, id=c('time'))
  
  p <- ggplot(m_df, aes(x = time, y = value, color = variable))
  plot <- p + geom_line() + ylab("y") + ggtitle('Two Var. Growth w/ Comp.')
  plot <- plot + scale_color_manual(values=cbbPalette)
  print(plot)
}

if (MULTI_VAR == TRUE) {
  #to increase number of populations, only change params and states
  NUM_VAR = 8
  rs <- c(.11,.12,.13,.14,.15,.16,.17,.18) #r_i's
  ks <- c(100,100,100,100,100,100,100,100) #k_i's
  as <- as.matrix(read.csv("comp_matrix.csv", header=FALSE)) #imports matrix of alpha_ij's
  Ns <- c(1,1,1,1,1,1,1,1) #initial states
  time <- seq(0,100,0.01) #time vector length and granularity
  
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
  m_df <- melt(df, id=c('time')) #turns ode output into melted data.frame
  
  p <- ggplot(m_df, aes(x = time, y = value, color = variable)) #basic plot
  plot <- p + geom_line() + ylab("y") + ggtitle('Eight Var. Growth w/ Comp.') #labels
  plot <- plot + scale_color_manual(values=cbbPalette) #better colors
  print(plot)
}

ptm2 <- proc.time()

print(ptm2 - ptm) #timing for efficiency

