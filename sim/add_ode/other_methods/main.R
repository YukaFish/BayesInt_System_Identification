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
dirname <- "results"
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
sed_lst <- read.table('seeds.txt')
sed_lst <- sed_lst[1:100,]
#### Experiments ####
for (ith in 1:r) {
  out <- generate.ydata(p = p, by = d_t, R = R, family = family,
                        dist_args = dist_args,
                        plot.ode = FALSE, plot.legend = FALSE,
                        plot_ydata = FALSE, data.seed = sed_lst[ith])

  truth <- out$truth
  ydata <- out$ydata
  times <- out$times
  source("experiments/iterations_multirep.R")
}

# read the final results
source('compare.R')
