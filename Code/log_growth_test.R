#The point of this file is to get logistic growth implemented
#as well as get a real ODE solver up and running (deSolve)

#install.packages("deSolve")
#install.packages("ggplot2")
#install.packages("reshape2")

library("deSolve")
library("ggplot2")
library("reshape2")

ONE_VAR <- FALSE
TWO_VAR <- FALSE
MULTI_VAR <- TRUE

if (ONE_VAR == TRUE) {
  r <- 0.1
  k <- 100
  N <- 1
  
  simple_model_params <- c(r = r, k = k)
  simple_state_val <- c(N = N)
  time <- seq(0,100,0.01)
  
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
  print(plot)
}

if (MULTI_VAR == TRUE) {
  NUM_VAR = 5
  rs <- c(.1,.1,.1,.1,.1)
  ks <- c(100,100,100,100,100)
  #as <- matrix(c(1,0.1,0.1,0.1,0.1,1,0.2,0.1,0.1,0.2,1,0.3,0.1,0.1,0.3,1),nrow=NUM_VAR, ncol=NUM_VAR)
  as <- matrix(c(1,0.1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,0.1,1,0.1,0.1,0.1,0.1,0.1,1)
               ,nrow=NUM_VAR, ncol=NUM_VAR)
  Ns <- c(1,5,10,15,20)
  time <- seq(0,100,0.01)
  
  params <- data.frame(rs,ks,as) #each row contains set of params for ODE
  states <- Ns
  
  multivar_log_comp <- function(time,state,params) {
    odes <- rep(NA,NUM_VAR)
    for (i in 1:NUM_VAR) {
      comp_sum <- 0
      for (j in 1:NUM_VAR) { #this might be inefficient compared to an inner product
        comp_sum <- comp_sum + params[i,2+j]*state[j] #but there are multiple type transformations
      } #that I would have to do
      Ni <- params[i,1]*state[i]*(1-(comp_sum/params[i,2]))
      odes[i] <- Ni
    }
    return(list(odes))
  }
  
  out <- ode(y=states, times=time, func=multivar_log_comp, parms=params)
  
  df <- as.data.frame(out)
  m_df <- melt(df, id=c('time'))
  
  p <- ggplot(m_df, aes(x = time, y = value, color = variable))
  plot <- p + geom_line() + ylab("y") + ggtitle('Five Var. Growth w/ Comp.')
  print(plot)
}





