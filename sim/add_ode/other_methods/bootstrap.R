##### MODIFY THE EXPERIMENT SETTINGS HERE ##########

args = commandArgs(trailingOnly = TRUE)

family = 'g'
N = 40
SNR=25
R = 1
n_tol = 1
NUM_REWEIGHT = 4

source(paste0("experiments/set_main_experiments.R"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#### Automatic arguments ####

# common arguments
r <- 100
# family arguments
family <- c("gaussian", "poisson", "binomial")[match(family, c("g", "p", "b"))]
log_lamth_min <- 0
log_lamth_max <- 0
lam_th_min <- 10^log_lamth_min
lam_th_max <- 10^log_lamth_max

# distribution-specific arguments
if (family == "gaussian") {
  dist_args <- list(SNR = SNR)
} else {
  dist_args <- NULL
}

# ode-generating arguments
p <- 10
t1 <- 0; t2 <- 20
N <- (t2 - t1) / d_t + 1

# main functions
source("R/data_generate_additive_scaled.R")
source("R/bcd_additive.R")
source("R/twostep_additive.R")
source("R/evaluate.R")
# transformer function
TRANSFORMER <- "sigmoid"
transformer <- pracma::sigmoid

# file arguments
header <- paste0(
  "method,neglik,theta.mse,dtheta.mse,",
  "f.mse,f.mse.t,f.mse.f,",
  "ftheta.mse,ftheta.mse.t,ftheta.mse.f,ftheta.mse.sumj,",
  "tpr,fpr"
)
# whether to keep existing records
appendfile <- 0
# file paths
dirname <- "results_boots"
if (!dir.exists(dirname)) dir.create(dirname)
more_suffix <- ""
if (family == "gaussian") more_suffix <- paste0("SNR", SNR)
filename <- paste0(
  dirname, "/", family, "_a_", "N", N - 1, "R", R,
  "NT", n_tol, "NW", NUM_REWEIGHT,
  more_suffix, ".txt"
)
if (!file.exists(filename)) {
  file.create(filename)
  write(header, filename)
} else {
  if (appendfile != 1 && appendfile != 0) appendfile <- 1
  appendfile <- as.logical(appendfile)
  if (!appendfile) {
    close(file(filename, open = "w"))
    write(header, filename)
  }
}

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
times <- seq(0, 20, length=N)
td <- seq(0, 20, length = 2001)
x <- deSolve::ode(state, times, ode.model, list(a0,a1))[,-1]
xdense <- deSolve::ode(state, td, ode.model, list(a0,a1))[,-1]

#### Bootstraps ####
sed_lst <- read.table('seeds.txt')
sed_lst <- sed_lst[1:100,]
for (jth in 1:length(sed_lst)) {
  out <- generate.ydata(p = p, by = d_t, R = R, family = family,
                        dist_args = dist_args,
                        plot.ode = FALSE, plot.legend = FALSE,
                        plot_ydata = FALSE, data.seed = sed_lst[jth])
  
  truth <- out$truth
  ydata <- out$ydata
  times <- out$time
  source("experiments/iterations_multirep.R")
  x_est <- (fda::eval.basis(out$gen$times, res1$inputs$basis_th) %*% res1$Cmat)
  res_est <- ydata[[1]] - x_est
  
  for (ith in 1:300) {
      # Calculate y*= yhat + epsilon_residual
      new_y = ydata[[1]]
      for (y_p in 1:p) {
        new_y[,y_p] = x_est[,y_p] + sample(res_est[,y_p], size = length(res_est[,y_p]), replace = T)
      }
      # Fit a model with new_y in place of the actual response
      ydata[[1]] = new_y
      # save the results and iterate
      ### JADE rw=1 -------------------------
      dirname1 <- paste0(dirname, "/", sed_lst[jth])
      if (!dir.exists(dirname1)) dir.create(dirname1)
      res1_list <- list()
      for (i in 1:num_rep) {
        res1_list[[i]] <- bcd.ode.additive(
          ydata, times,
          family = family,
          dist_args = dist_args,
          lam_th_min = lam_th_min, lam_th_max = lam_th_max,
          delta_lam_th_init = 1,
          delta_update_method = "decrease",
          lam_gam_list = lam_gam_list,
          lam_gam_ratio = "auto",
          delta_lam_gam = 10^0.2,
          knots.theta.dt = 0.1,
          scaling_factor = 1,
          norder_f = 4, nknots_f = 6,
          knots_position = "quantile",
          init_method = init_method,
          lam_gam_init = lam_gam_init,
          lam_gam_init2 = lam_gam_init2,
          penalty = "grpLasso",
          penalty.weighted = TRUE,
          smooth_package = "stats",
          theta_method = "1dH",
          initialize_penalty_weights = initw,
          MAXITER = MAXITER,
          NUM_REWEIGHT = NUM_REWEIGHT,
          final_result_only = FALSE,
          random_order = TRUE,
          tuning_tolerant = n_tol,
          BICobj = "ode",
          eps.conv = 1e-4,
          ada.gam = 1,
          L_lower_bound_factor = 1e-8,
          know_truth = FALSE,
          fit.trajectory = FALSE,
          mode.test = FALSE,
          verbose = FALSE,
          lam.message = FALSE
        )
      }
      
      BICs <- sapply(res1_list, \(x) x$BICvalue)
      minbic.id <- which.min(BICs)
      res1_temp <- res1_list[[minbic.id]]
      
      saveRDS(res1_temp, file=paste0(dirname1, "/", 'res1.',sed_lst[jth],'.', ith, '.rds'))
  }
}

# help function for calculating the confidence intervals
get.bcd.ode.fidelity <- function(res, type = c("all", "data", "ode")) {
  repRbind <- function(v, R) {
    return(do.call(rbind,replicate(R,v,simplify=FALSE)))
  }
  keepDim <- function(a, dims) {
    return(apply(a,dims,as.vector))
  }
  neg_log_likelihood <- function(ydata, theta_est, fam) {
    b <- fam$b
    # negative log-likelihood for exp.family
    theta_est <- as.matrix(theta_est); ydata <- as.matrix(ydata)
    nt <- dim(theta_est)[1]; p <- dim(theta_est)[2]
    if (dim(ydata)[1]!=nt || dim(ydata)[2]!=p) stop("preds, y, df not matched")
    btheta_est <- b(theta_est)
    return(-1/(nt*p) * sum(ydata * theta_est - btheta_est))  # or -1/(nt)
  }
  source("R/utils/transformer.R", local = TRUE)
  type <- match.arg(type)
  ydata <- res$inputs$ydata
  Hmat <- res$inputs$Hmat
  H4int <- res$inputs$H4int
  dH4int <- res$inputs$dH4int
  p  <- dim(ydata)[2]; nt <- dim(ydata)[1]
  R <- res$inputs$R
  rm_id <- res$inputs$rm_id
  
  theta_bspl <- Hmat %*% res$Cmat
  theta_bspl_4int <- H4int %*% res$Cmat
  dtheta_bspl_4int <- dH4int %*% res$Cmat
  Bmat <- form_Bmat(theta_bspl_4int, res$inputs$basis_f, normalized = FALSE,
                    center = res$inputs$center_th, scale = res$inputs$scale_th)
  dtheta_ode_4int <- 
    do.call(cbind,
            lapply(1:p, function(j) {
              res$gam0[j] + (Bmat %*% as.vector(res$gam1[,,j]))
            } )
    )
  times = seq(0,20,length=41)
  t4int = seq(0,20,length=2001)
  eid <- sapply(times, function(x) {
    ind <- which(abs(t4int - x) < 1e-6)
    ind <- ifelse(length(ind) == 0, NA, ind)
  })
  x_est <- matrix(NA,nrow=41,ncol=p)
  x_est[1,] = dtheta_ode_4int[1,]*0.01
  for (jj in 1:p) {
    for (ii in 1:(length(eid)-1)) {
      x_est[ii+1,jj] <- sum(dtheta_ode_4int[1:(eid[ii+1]-1),jj]*0.01)
    }
    x_est[,jj] = scale(x_est[,jj],T,F)
  }
  out <- list()
  if (type %in% c("all", "data")) {
    fidelity.data <-
      neg_log_likelihood(ydata,
                         do.call(rbind, replicate(R, theta_bspl, FALSE)), res$inputs$fam)
    out <- append(out, list(fidelity.data = fidelity.data))
  }
  if (type %in% c("all", "ode")){
    fidelity.ode <- Metrics::mse(dtheta_ode_4int, dtheta_bspl_4int)
    out <- append(out, list(fidelity.ode = fidelity.ode))
    out <- append(out, list(x_est = x_est))
    
  }
  return(out)
}

# calculate the coverage probabilities 
dirname <- "results_boots"
dirname1 <- Sys.glob(paste0(dirname, "/*",sep=""))
cov_prob_lst <- array(NA, dim=c(length(dirname1)-1,N,num_sigma))
for (dir in 1:(length(dirname1)-1)) {
  new_fit_lst <- Sys.glob(paste(dirname1[dir], '/res1.*.rds', sep=""))
  x_est_lst <- list()
  for (new_fit_i in 1:length(new_fit_lst)) {
    x_est_lst[[new_fit_i]] = get.bcd.ode.fidelity(readRDS(new_fit_lst[new_fit_i]),type='all')
  }
  x_est_lst_ma = array(NA,dim=c(length(new_fit_lst),N,10))
  for (new_fit_i in 1:length(new_fit_lst)) {
    x_est_lst_ma[new_fit_i,,] = x_est_lst[[new_fit_i]]$x_est
  }
  x_est_upper = apply(x_est_lst_ma,c(2,3),quantile,0.975)
  x_est_lower = apply(x_est_lst_ma,c(2,3),quantile,0.025)
  cov_prob_lst[dir,,] = (x_est_upper > scale(x,scale=F)) & (x_est_lower < scale(x,scale=F))
}
saveRDS(cov_prob_lst,file='cov_prob_lst.rds')

# plot the coverage probabilities
library(ggplot2)
library(cowplot)
library(latex2exp)
plot_cov_traj_prob <- function(times, cov_traj_lst) {
  p_lst <- list()
  for (i in 1:num_sigma) {
    cov_prob <- data.frame(times = times,freq = apply(cov_traj_lst[,,i], 2, sum)/dim(cov_traj_lst)[1])
    p<-ggplot(cov_prob, aes(x=times, y = freq)) + 
      geom_hline(yintercept=0.95, colour = "red")+
      ylim(c(0,1))+
      geom_line(colour = "black")+
      theme(legend.position="none", axis.title = element_text(size = 15),
            plot.title = element_text(size = 17, face = "bold"),
            axis.text = element_text(size = 13),
            panel.background = element_rect(fill = NA),
            panel.grid.major = element_line(colour = "grey90"))+
      xlab('time') + ylab('probabilities')+ggtitle(TeX(paste('$x_{', i, '}$',sep="")))
    p_lst[[i]] <- p
  }
  p_lst
}
cov_prob_lst <- readRDS('cov_prob_lst.rds')
dev.new(width=10, height=5, noRStudioGD = TRUE)
p_lst <- plot_cov_traj_prob(seq(0,20,length=N), cov_prob_lst)
plot_grid(plotlist=p_lst, ncol=5)

