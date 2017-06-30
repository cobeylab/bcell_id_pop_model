#Implementation of basic Z+A 2015 model just to make sure I remember
#exactly what I'm doing here

#Graham Northrup June 30

NUMDAYS <- 60
dt <- 0.05
NUMSTEPS <- NUMDAYS/dt
Hf <- rep(0, NUMSTEPS) #free antigen
Hb <- rep(0, NUMSTEPS) #bound antigen
B <- rep(0, NUMSTEPS) #B cells
A <- rep(0, NUMSTEPS) #Antibodies

k <- 0.01 #1/AU*day
dH <- 0.5 #1/day
s <- 1 #1/day
phi <- #AU
a <- 0.1 #1/day
dA <- 0.1 #1/dat

Hf[1] <- 1000

for (i in 1:NUMSTEPS) {
  
}