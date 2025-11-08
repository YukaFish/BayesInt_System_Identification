obs_41 = list(times_by=0.5, obs=41, basis_by=0.25, out_nquad=200)
# choose the sampling frequency and SNR to get different results
obs_choice <- obs_41
snr <- 25

# simulate the ODE system
source('sim_system.R')
# load preparation file for simulation
source('prep.R')
# load the help function
source('help_func.R')
# load simulation file
source("BayesODE.R")
# read the results
source("compare.R")

