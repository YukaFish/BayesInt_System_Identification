# load required packages
library(CollocInfer)
library(gss)
library(fda)
library(deSolve)
library(NLRoot)
library(TMB)
library(tmbstan)
library(grpreg)
library(foreach)
library(doParallel)

eta = c(20,0.25); lambda = 1000

# use B-spline to estimate trajectory
knots = seq(times[1], max(times),by = obs_choice$basis_by); x_knots = 1
nknots = length(knots)
norder = 4
nbasis = length(knots) + norder - 2
bsbasis = create.bspline.basis(range(knots),nbasis,norder,knots)
# basis values at sampling points
basismat = eval.basis(times, bsbasis)

# number of outer and inner quadrature points in total
out_nquad = obs_choice$out_nquad; in_nquad = 2
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

