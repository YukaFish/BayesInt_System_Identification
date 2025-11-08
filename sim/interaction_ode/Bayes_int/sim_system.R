# load required packages
library(CollocInfer)
library(gss)
library(fda)
library(deSolve)
library(NLRoot)
library(TMB)
library(tmbstan)
library(grpreg)

# set the seed
set.seed(123)
# define the number of sigma
num_sigma = 10
# simulate the true system
ode.model <- function(t, state, parameters) {
  parameters <- unname(parameters)
  state <- unname(state)
  for (j in 1:5) {
    assign(paste("dstate", 2*j-1, sep = ""), 
           parameters[1,j]*state[2*j-1]-parameters[2,j]*state[2*j-1]*state[2*j])
    assign(paste("dstate", 2*j, sep = ""), 
           parameters[3,j]*state[2*j-1]*state[2*j]-parameters[4,j]*state[2*j])
  }
  return(list(c(dstate1,dstate2,dstate3,dstate4,dstate5,
                dstate6,dstate7,dstate8,dstate9,dstate10)))
}
for (j in 1:5) {
  assign(paste("pair", j, sep = ""), c(1.1+0.2*(j-1), 0.4+0.2*(j-1), 
                                       0.1+0.2*(j-1), 0.4+0.2*(j-1)))
}
pair <- cbind(pair1,pair2,pair3, pair4,pair5)
state <- rep(c(5,4,7.5,5.5,8),each=2) 

# simulate the system
times <- seq(0, 20, length=obs_choice$obs)
td <- seq(0, 20, length = 2001)
x <- deSolve::ode(state, times, ode.model, pair)[,-1]
xdense <- deSolve::ode(state, td, ode.model, pair)[,-1]

# save the simulated system for future use
saveRDS(xdense, 'xdense.rds')
saveRDS(x, paste('x.obs', obs_choice$obs, '.rds',sep=""))

