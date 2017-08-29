library('deSolve')
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping

ODE <- 0

ac1 <- .001
ac2 <- .0005
inits <- c(B1 = 0.4, B2 = 0.6, t=0)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7") #colorblind friendly colors
TIME <- 28
deltaT <- 0.1
t0 <- 0
nstep <- TIME/deltaT
time <- seq(0,TIME,deltaT)

m_params <- c(ac1=ac1,ac2=ac2)
m_state <- inits

if (ODE){
  markov_pop <- function(t,state,params){
    with(as.list(c(state, params)),{
      P11 <- (1 - (ac1)/(ac1+ac2))*exp(-(ac1+ac2)*t) + (ac1/(ac1+ac2))
      P12 <- (ac2/(ac1+ac2))*(1 - exp(-(ac1+ac2)*t))
      P22 <- (1 - (ac2)/(ac1+ac2))*exp(-(ac1+ac2)*t) + (ac2/(ac1+ac2))
      P21 <- (ac1/(ac1+ac2))*(1 - exp(-(ac1+ac2)*t))
      d1 <- P11 + P21
      d2 <- P22 + P12
      print(d1/2)
      print(d2/2)
      list(c(d1,d2))
    })
  }
  
  out <- ode(y=m_state, times=time, func=markov_pop, parms=m_params)
  
  df <- as.data.frame(out)
  m_df <- melt(df, id=c('time'))
  
  p <- ggplot(m_df, aes(x = time, y = value, color = variable))
  plot <- p + geom_line() + ylab("y") + ggtitle('Niche Frac. Growth')
  plot <- plot + scale_color_manual(values=cbbPalette)
  print(plot)
}

ds <- matrix(ncol = 3, nrow = nstep)
ds[1,] = inits

for (i in 2:length(ds[,1])){
  P11 <- (1 - (ac1)/(ac1+ac2))*exp(-(ac1+ac2)*time[i]) + (ac1/(ac1+ac2))
  P12 <- (ac2/(ac1+ac2))*(1 - exp(-(ac1+ac2)*time[i]))
  P22 <- (1 - (ac2)/(ac1+ac2))*exp(-(ac1+ac2)*time[i]) + (ac2/(ac1+ac2))
  P21 <- (ac1/(ac1+ac2))*(1 - exp(-(ac1+ac2)*time[i]))
  ds[i,1] <- (P11*ds[i-1,1] + P21*ds[i-1,2])
  ds[i,2] <- (P12*ds[i-1,1] + P22*ds[i-1,2])
  ds[i,3] <- time[i]
}

df <- as.data.frame(ds)
colnames(df) <- c('Pop. 1','Pop. 2', 'time')
m_df <- melt(df, id=c('time'))

p <- ggplot(m_df, aes(x = time, y = value, color = variable))
plot <- p + geom_line() + ylab("y") + ggtitle('Niche Frac. Growth')
print(plot)