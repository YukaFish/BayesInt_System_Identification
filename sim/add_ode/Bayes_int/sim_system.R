# load required packages
library(CollocInfer)
library(gss)
library(fda)
library(deSolve)
library(NLRoot)
library(TMB)
library(tmbstan)
library(grpreg)

# set the seed to ensure the x7-x10 be the same every time
set.seed(123)
# define the number of sigma
num_sigma = 10
# simulate the true system: adopt the code from
# Zhang, N., Nanshan, M., and Cao, J. (2022). 
# A joint estimation approach to sparse additive
# ordinary differential equations. Statistics and Computing, 32(5):69.505
poly3 <- function(x) return(c(x,x^2,x^3))
ode.model <- function(t, state, parameters) {
  a0 <- parameters[[1]]
  a1 <- parameters[[2]]
  state1 <- c(sapply(state[1:p1], poly3))
  dstate <- c(a0 + as.vector(state1 %*% a1), a0p)
  return(list(dstate))
}

p0 <- 10; p1 <- 6; p2 <- p0-p1
a10 = 0; a11 = c(1.2, 0.3, -0.6); a12 = c(0.1, 0.2, 0.2)
a20 = 0.4; a21 = c(-2, 0, 0.4); a22 = c(0.5, 0.2, -0.3)
a30 = -0.2; a33 = c(0, 0, 0); a34 = c(0.3, 0.4, 0.1)
a40 = -0.2; a43 = c(0.2,-0.1,-0.2); a44 = c(0, 0, 0)
a50 = 0.05; a55 = c(0, 0, 0); a56 = c(0.1, 0, -0.8)
a60 = -0.05; a65 = c(0, 0, 0.5); a66 = c(0, 0, 0)

a0 <- c(a10,a20,a30,a40,a50,a60)
A1 <- rbind(c(a11,a12),c(a21,a22))
A2 <- rbind(c(a33,a34),c(a43,a44))
A3 <- rbind(c(a55,a56),c(a65,a66))
a1 <- t(as.matrix(Matrix::bdiag(list(A1,A2,A3))))

a0p <- runif(p2, min=-0.2, max=0.2)
state  <- c(c(-2,2,2,-2,-1.5,1.5), rnorm(p2))

# simulate the system
times <- seq(0, 20, length=obs_choice$obs)
td <- seq(0, 20, length = 2001)
x <- deSolve::ode(state, times, ode.model, list(a0,a1))[,-1]
xdense <- deSolve::ode(state, td, ode.model, list(a0,a1))[,-1]

# save the simulated system for future use
saveRDS(xdense, 'xdense.rds')
saveRDS(x, paste('x.obs', obs_choice$obs, '.rds',sep=""))

