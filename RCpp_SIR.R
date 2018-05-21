#saralamba@mail.com
#
#set working directory
setwd("D:/Data-Work/MAEMOD/RcppSIR")

library(Rcpp)
library(deSolve)

#compile the model code
#this command will create a function SIRmodel(t, state, parameters) at the  Global env.
   sourceCpp("SIR.cpp")

#but if you would like to compile this Rcpp code ("SIR.cpp") as the dynamic link library (.dll, .o or .so)
#you can do it by using this command  

#system("R CMD SHLIB SIR.cpp")

#if you have the compiling problem about Rcpp.h not found then try to create Makevars file which has these two lines
#   PKG_CPPFLAGS=`Rscript -e 'Rcpp:::CxxFlags()'`
#   PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'` 
#and put this Makevars file in the same path as SIR.cpp and run system(R CMD SHLIB SIR.cpp) again
#this SIR.dll/SIR.o or SIR.so shoud be created.
#to load this dll file type
#   dyn.load("SIR.dll")

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

## Rcpp version
start.time <- Sys.time()

results.Rcpp <- matrix(NA,nrow=100,ncol=length(times))
for(i in 1:100){
  parms <- c(beta = 0.1 + rnorm(1,0,0.1), recovery = 0.005, death = 0.001, birth = 0.001)
  out <- ode(y = initial, times = times, func = SIRmodel, parms = parms)
  results.Rcpp[i,] <- out[,"I"]
} 
cols<-rainbow(100)
plot(results.Rcpp[1,],type='l',ylim = c(0,1))
for(i in 2:100) lines(results.Rcpp[i,],col=cols[i])

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


### R version
start.time <- Sys.time()

results.R <- matrix(NA,nrow=100,ncol=length(times))
for(i in 1:100){
  parms <- c(beta = 0.1 + rnorm(1,0,0.1), recovery = 0.005, death = 0.001, birth = 0.001)
  out <- ode(y = initial, times = times, func = SIRmodel.Rversion, parms = parms)
  results.R[i,] <- out[,"I"]
} 
cols<-rainbow(100)
plot(results.R[1,],type='l',ylim = c(0,1))
for(i in 2:100) lines(results.R[i,],col=cols[i])

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



## Rcpp version
start.time <- Sys.time()
for(i in 1:100){
  out <- ode(y = initial, times = times, func = SIRmodel, parms = parms)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


## R version
start.time <- Sys.time()
for(i in 1:100){
  out <- ode(y = initial, times = times, func = SIRmodel.Rversion, parms = parms)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



