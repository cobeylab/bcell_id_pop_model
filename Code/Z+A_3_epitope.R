#Three epitope Z+A model

library("deSolve") #ode solver
library("ggplot2") #good data vis.
library("reshape2") #data.frame reshaping

k <- 0.01 #1/AU*day
df <- 0.5 #1/day
db <- df #1/day 3 FOR ACM , 0.5 o/w
s <- 1 #1/day
phi <- 10 #AU
a <- 0.1 #1/day
dA <- 0.1 #1/dat
alpha <- 0.0 #0.01 FOR FIM 0.0 o/w
beta <- 0.95 #steric interference
delta <- 0 #ZERO FOR EMM 1 o/w

Hxys <- 1000
Hoys <- 0
Hxos <- 0
Hxyo <- 0
Hoos <- 0
Hoyo <- 0
Hxoo <- 0
Hooo <- 0
Bx <- 1
By <- 1
Bs <- 1
Ax <- 1
Ay <- 1
As <- 512

time <- seq(0,60,0.01)

params <- c(k=k,df=df,db=db,s=s,phi=phi,a=a,dA=dA,alpha=alpha,beta=beta,delta=delta)
states <- c(Hxys=Hxys,Hoys=Hoys,Hxos=Hxos,Hxyo=Hxyo,Hoos=Hoos,Hoyo=Hoyo,
            Hxoo=Hxoo,Hooo=Hooo,Bx=Bx,By=By,Bs=Bs,Ax=Ax,Ay=Ay,As=As)

z_a <- function(t, states, params){
  with(as.list(c(states, params)),{
    dHxys <- -k*Hxys*(Ax+Ay+As) - df*Hxys
    dHoys <- k*((Hxys*Ax) - Hoys*((1-beta)*Ay + As)) - db*Hoys
    dHxos <- k*((Hxys*Ay) - Hoys*((1-beta)*Ax + As)) - db*Hxos
    dHxyo <- k*((Hxys*As) - Hxyo*(Ax+Ay)) - db*Hxyo
    dHoos <- k*((1-beta)*(Hoys*Ay + Hxos*Ax) - Hoos*As) - db*Hoos
    dHoyo <- k*(Hoys*As + Hxyo*Ax - (1-beta)*Hoyo*Ay) - db*Hoyo
    dHxoo <- k*(Hxos*As + Hxyo*Ay - (1-beta)*Hxoo*Ax) - db*Hxoo
    dHooo <- k*(Hoos*As + (1-beta)*(Hoyo*Ay + Hxoo*Ax)) - db*Hooo
    xsum <- Hxyo + (1-beta)*(Hxos + Hxoo) + delta*(Hoys + Hoyo + Hoos + Hooo)
    dBx <- s*Bx*((Hxys + xsum)/(phi + Hxys + xsum))*(1/(1+alpha*xsum))
    ysum <- Hxyo + (1-beta)*(Hoys + Hoyo) + delta*(Hxos + Hxoo + Hoos + Hooo)
    dBy <- s*By*((Hxys + ysum)/(phi + Hxys + ysum))*(1/(1+alpha*ysum))
    ssum <- Hoys + Hxos + Hoos + delta*(Hxyo + Hoyo + Hxoo + Hooo)
    dBs <- s*Bs*((Hxys + ssum)/(phi + Hxys + ssum))*(1/(1+alpha*ssum))
    dAx <- a*Bx - k*(Hxys + Hxyo + (1-beta)*(Hxos + Hxoo))*Ax - dA*Ax
    dAy <- a*By - k*(Hxys + Hxyo + (1-beta)*(Hoys + Hoyo))*Ay - dA*Ay
    dAs <- a*Bs - k*(Hxys + Hoys + Hxos + Hoos)*As - dA*As
    list(c(dHxys,dHoys,dHxos,dHxyo,dHoos,dHoyo,dHxoo,dHooo,dBx,dBy,dBs,dAx,dAy,dAs))
  })
}

out <- ode(y=states, times=time, func=z_a, parms=params)

df <- as.data.frame(out)
m_df <- melt(df, id=c('time'))

p <- ggplot(m_df, aes(x = time, y = value, color = variable))
plot <- p + geom_line() + ylab("y") + ggtitle('3 epitope model')
plot <- plot + scale_color_manual(values=cbbPalette)
plot <- plot %+% subset(m_df, variable %in% c('Bx','By','Bs'))
print(plot)
