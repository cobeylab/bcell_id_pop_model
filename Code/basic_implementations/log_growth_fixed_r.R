#Exploring the bredth of behaviors possible with symmetric r terms

#setwd("~/Desktop/CobeyLab/bcell_id_pop_model/Code/basic_implementations")

library("deSolve") #ode solver
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping

cbbPalette <- c("#000000", "#E69F00", "#56B4E9",
                "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7") #colorblind friendly colors

NUM_VAR = 5
r <- 0.3
rs <- rep(r, NUM_VAR) #r_i's
ks <- c(100,110,120,130,140) #k_i's
as <- as.matrix(read.csv("fixed_r.csv", header=FALSE)) #imports matrix of alpha_ij's
Ns <- c(1,1,1,1,1) #initial states
time <- seq(0,100,0.1) #time vector length and granularity

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

p <- ggplot(m_df, aes(x = time, y = value, color = variable))
plot <- p + geom_line() + ylab("y") + ggtitle('Five Var. Growth w/ Fixed r')
plot <- plot + scale_color_manual(values=cbbPalette)
print(plot)