# compile required cpp files
compile("spline_int.cpp")
dyn.load(dynlib("spline_int"))

# read the seed list
sed_lst <- read.table('seeds.txt')
sed_lst <- sed_lst[1:100,]

for (j in 1:length(sed_lst)){
  ############## simulate observations #############
  set.seed(sed_lst[j])
  sdy <- apply(x, 2, sd)/snr
  ydata <- replicate(1, x + t(replicate(length(times), 
                     rnorm(num_sigma, mean = 0, sd = sdy))),
                     simplify = FALSE)
  y <- as.matrix(ydata[[1]])
  
  ############## initial values #############
  # initial value for b: smooth the coefficient of basis
  FhNfdPar <- fdPar(bsbasis, int2Lfd(2), 0.01)
  DEfd0 <- smooth.basis(times, as.matrix(y), FhNfdPar)$fd
  b <- DEfd0$coefs
  xx <- I0quadbasismat%*%b
  ix <- I1quadbasismat%*%b
  norm_ix <- normalize(ix)
  ix <- 1/(1+exp(-norm_ix))
  
  # initial value for a, theta
  for (i in 1:num_sigma) {
    assign(paste("int_bs_ix", i, sep = ""), integral_bs(cal_basis(ix[,i],x_knots)))
  }
  pos.len <- seq(0,1,length.out=x_knots+2)
  pos <- pos.len[-c(1,length(pos.len))]
  int_bs_ix <- scale(out_quadpts,scale=F)
  for (ii in 1:num_sigma) {
    int_bs_ix <- cbind(int_bs_ix, get(paste("int_bs_ix", i, sep = "")))
  }
  a0 <- rep(NA, num_sigma)
  a <- array(NA, dim=c(dim(int_bs_ix1)[2], num_sigma, num_sigma))
  graph = matrix(NA, nr=num_sigma,nc=num_sigma)
  for (i in 1:num_sigma) {
    result <- est_coef(scale(xx[,i],scale=F), int_bs_ix)
    a[,i,] <- result$beta
    a0[i] <- result$a0
    graph[i,] <- result$graph
  }
  theta = matrix(rep(0,num_sigma*num_sigma),nc=num_sigma)
  theta[!graph] <- 1e-10; theta[graph] <- 1-1e-10
  
  ############## run the methodology #############
  # define the initial function
  init_logsigma <- rep(log(0.1), num_sigma)
  init_fn <- function() list(logsigma=init_logsigma, b=b, a = a,
                             theta=theta,a0=a0)
  # set the lower bound for parameter
  lower <- c(
    rep(-Inf, length(b)+num_sigma),
    rep(-Inf, length(a)),
    matrix(rep(0,num_sigma*num_sigma),ncol=num_sigma),
    rep(-Inf, num_sigma)
  )
  # set the upper bound for parameter
  upper <- c(
    rep(Inf, length(b)+num_sigma),
    rep(Inf, length(a)),
    matrix(rep(1,num_sigma*num_sigma),ncol=num_sigma),
    rep(Inf, num_sigma)
  )
  # set the parameters
  parameters <- list(logsigma=rep(log(0.1),num_sigma), b=b,
                     a=a, theta=theta,a0=a0
  )
  new_data <- list(y=y, lambda = lambda,
                   I0quadbasismat = as(I0quadbasismat-phi0basismat, "dgCMatrix"),
                   I1quadbasismat = as(I1quadbasismat,'dgCMatrix'), 
                   basismat = as(basismat,'dgCMatrix'),
                   out_nquad = out_nquad,
                   pos = pos, 
                   out_quadpts = out_quadpts, 
                   out_quadwts = out_quadwts, 
                   in_quadwts_vec = in_quadwts_vec,eta=eta)
  options(mc.cores = parallel::detectCores())
  obj <- MakeADFun(new_data, parameters, DLL="spline_int")
  # run the model
  new_fit = tmbstan(obj,
                    lower = lower, upper = upper,
                    init = init_fn,
                    chains = 1, iter = 1000,laplace = FALSE,
                    control = list(adapt_delta = 0.85, max_treedepth = 12))
  saveRDS(new_fit,file = paste('sim.result.', obs_choice$obs, 'obs.', snr, 'snr.', sed_lst[j], '.rds', sep=""))
}

