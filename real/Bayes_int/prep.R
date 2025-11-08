# load required packages
library(CollocInfer)
library(gss)
library(fda)
library(deSolve)
library(NLRoot)
library(TMB)
library(tmbstan)
library(grpreg)
library(splines)
library(orthogonalsplinebasis)
library(MASS)
library(grplasso)
library(lars)
library(longitudinal)
data(tcell)

# define the number of observations and times
num_sigma = 58; times <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72)

y <- array(NA, dim = c(length(times),num_sigma,34))
for (ii in 1:58) {
  curr.sample <- t(matrix(tcell.34[,ii], nc=length(times)))
  y[,ii,] <- curr.sample
}

# set the value for eta and lambda
eta = c(10,0.6); lambda = 10

# use B-spline to estimate trajectory
knots <- c(0, 2, 4, 6, 8, 18, 24, 32, 48, 72); x_knots = 2
nknots = length(knots)
norder = 4
nbasis = length(knots) + norder - 2
bsbasis = create.bspline.basis(range(knots),nbasis,norder,knots)
# basis values at sampling points
basismat = eval.basis(times, bsbasis)

# number of outer and inner quadrature points in total
out_nquad = 30; in_nquad = 2
# set up outer integral quadrature points and weights
out_quadre = gauss.quad(out_nquad, c(min(times), max(times)))
out_quadpts = out_quadre$pt
out_quadwts = out_quadre$wt
# set up inner integral quadrature points and weights
in_quadpts = matrix(NA, nc = length(out_quadpts),nr = in_nquad)
in_quadwts = matrix(NA, nc = length(out_quadpts),nr = in_nquad)
in_quadre = gauss.quad(in_nquad, c(min(times), out_quadpts[1]))
in_quadpts[,1] = in_quadre$pt
in_quadwts[,1] = in_quadre$wt
for (k in 1:(length(out_quadpts)-1)) {
  in_quadre = gauss.quad(in_nquad, c(out_quadpts[k], out_quadpts[k+1]))
  in_quadpts[,k+1] = in_quadre$pt
  in_quadwts[,k+1] = in_quadre$wt
}
in_quadpts_vec = as.vector(in_quadpts)
in_quadwts_vec = as.vector(in_quadwts)

# values of the integral of basis functions at quadrature points
phi0basis = eval.basis(0, bsbasis, 0)
I0quadbasismat = eval.basis(out_quadpts, bsbasis, 0)
I1quadbasismat = eval.basis(in_quadpts_vec, bsbasis, 0)
phi0basismat = t(matrix(rep(phi0basis, times = length(out_quadwts)),
                        nc = length(out_quadwts)))
