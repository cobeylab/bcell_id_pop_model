#UNFINISHED

ac1 <- .2
ac2 <- .1
inits <- c(0.5,0.5)
t <- 0
index <- 2
MAX_INDEX <- 100000
fracs <- matrix(ncol = MAX_INDEX, nrow = 2)
fracs[,1] = inits

while(index <= MAX_INDEX){
  P11 <- (1 - (ac1)/(ac1+ac2))*exp(-(ac1+ac2)*t) + (ac1/(ac1+ac2))
  P12 <- (ac2/(ac1+ac2))*(1 - exp(-(ac1+ac2)*t))
  P22 <- (1 - (ac2)/(ac1+ac2))*exp(-(ac1+ac2)*t) + (ac2/(ac1+ac2))
  P21 <- (ac1/(ac1+ac2))*(1 - exp(-(ac1+ac2)*t))
  as <- c(P11*fracs[1,index-1],P12*fracs[1,index-1],P22*fracs[2,index-1],P21*fracs[2,index-1])
  Atot <- sum(as)
  r1 <- runif(1)
  tau <- -log(r1)/Atot
  r2 <- runif(1)
}