#saralamba@gmail.com
#
#set working directory
setwd("D:/Data-Work/GitHub/Rcpp_deSolve")

library(Rcpp)
library(deSolve)

## SIR Rcpp version
sourceCpp("SIR.cpp")

parameters <- c(beta = 0.1, recovery = 0.005, death = 0.001, birth = 0.001)
initial <- c(S = 0.9, I = 0.1, R = 0)
times <- seq(0,300,by = 0.1)

#solve ODE and plot some output parameters  
out <- ode(y = initial, times = times, func = SIRmodel, parms = parameters)
plot(out, select = c("I","extraParm"))

## SIR - R version
SIRmodel.Rversion <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{

    dS  <- birth - beta*I*S - death*S
    dI  <- beta*I*S - recovery*I - death*I
    dR  <- recovery*I - death*R
    
    extraParm = beta*I*S;
    list( c(dS,dI,dR), extraParm = extraParm )
  })
}

out.Rversion <- ode(y = initial, times = times, func = SIRmodel.Rversion, parms = parameters)
plot(out.Rversion)


### Speed comparison
microbenchmark::microbenchmark(
  rcpp = ode(y = initial, times = times, func = SIRmodel, parms = parameters, method = 'rk4'),
  r = ode(y = initial, times = times, func = SIRmodel.Rversion, parms = parameters, method = 'rk4'),
  times = 100
)

#Unit: milliseconds
#expr      min       lq     mean   median       uq      max neval
#rcpp 288.1454 290.7679 298.1657 293.0336 298.7643 375.0034   100
#r 272.2395 276.3817 284.3532 278.3128 282.0885 504.2865   100


### SIRS

# Rcpp version
sourceCpp("SIRS.cpp")
times <- seq(0,10,by = 1/365)
initial <- c(S=999999, I=1, R=0)
parameters <- c(R0=3, durinf=1/52, durimm=4, P=1000000, startvac=3,durvac=2,covvac=20,starttreat =1,covtreat=50,curetime=2/365)

out <- ode(y=initial,times = times, func = SIRSmodel,parms = parameters)
plot(out)

# R version
SIRSmodel.Rversion <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    beta <- R0*durinf
    omega <- 1/durimm
    
    theta <-0
    ##conditions for theta
    if(t < startvac)
      theta <- 0
    if((startvac <= t) && (t < (startvac+durvac)))
      theta <- -(1/durvac)*log(1-(covvac/100))
    if(t > (startvac+durvac))
      theta <- 0
    
    ##conditions for nu
    if(t < starttreat)
      nu <- 1/durinf
    if(t >= starttreat)
      nu <- (1-(covtreat/100))*(1/durinf) + (covtreat/100)*(1/curetime)
    
    dS <- -beta*S*I/P + omega*R - theta*S
    dI <-  beta*S*I/P - nu*I 
    dR <- nu*I - omega*R + theta*S                    
    
    prev <- 100*I/(S+I+R)
    pop <- S+I+R
    
    #case per 1000 per year
    inc <- 1000*beta*S*I/P 
    
    list( c(dS,dI,dR), prev = prev, pop = pop, inc =inc)
  })
}

out.Rversion <- ode(y=initial, times = times, func = SIRSmodel.Rversion, parms = parameters)
plot(out.Rversion)

## speed comparison

microbenchmark::microbenchmark(
  rcpp = ode(y=initial, times = times, func = SIRSmodel, parms = parameters, method = 'rk4'),
  r = ode(y=initial, times = times, func = SIRSmodel.Rversion, parms = parameters, method = 'rk4'),
  times = 100
)

#Unit: milliseconds
#expr      min       lq     mean   median       uq      max neval
#rcpp 619.2207 624.1921 641.6795 626.6914 641.0786 811.7491   100
#r 626.3331 631.5665 668.4641 635.2158 674.3270 967.7877   100


