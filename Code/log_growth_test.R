#The point of this file is to get logistic growth implemented
#as well as get a real ODE solver up and running (deSolve)

#install.packages("deSolve")

library("deSolve")

ONE_VAR <- FALSE
TWO_VAR <- TRUE

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
  plot(out, xlab='time',ylab='number of things',main = header,cex.main=0.7)
}

if (TWO_VAR == TRUE) {
  rs <- c(.1,.2)
  ks <- c(100,100)
  as <- c(0.01,0.01)
  Ns <- c(1,1)
  time <- seq(0,100,0.01)
  
  tv_params <- c(r1=rs[1],r2=rs[2],k1=ks[1],k2=ks[2],a12=as[1],a21=as[2])
  tv_state <- c(N1=Ns[1],N2=Ns[2])
  
  tv_log_compete <- function(time,state,params) {
    with(as.list(c(state, params)),{
      dN1 <- r1*N1*(1-((N1-a12*N2)/k1))
      dN2 <- r2*N2*(1-((N2-a21*N1)/k2))
      list(c(dN1,dN2))
    })
  }
  
  out <- ode(y=tv_state, times=time, func=tv_log_compete, parms=tv_params)
  plot(out[, "time"], out[,"N1"],type='l', col='red',xlab='time',ylab='number',
       main = 'Two Pop. Logisitc Growth w/ Competition', cex.main=0.7)
  lines(out[, "time"], out[,"N2"], col='blue')
}





