#set working directory
setwd("D:/Data-Work/MAEMOD/RcppSIR")

library(Rcpp)
library(deSolve)

#compile the model code
sourceCpp("SIR.cpp")

parameters <- c(beta = 0.1, recovery = 0.005, death = 0.001, birth = 0.001)
initial <- c(S = 0.9, I = 0.1, R = 0)
times <- seq(0,300,by = 0.1)

#solve ODE and plot some output parameters  
out <- ode(y = initial, times = times, func = SIRmodel, parms = parameters)
plot(out, select = c("I","extraParm"))
