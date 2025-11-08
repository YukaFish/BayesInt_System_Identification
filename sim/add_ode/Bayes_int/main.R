# Note: the default case is 41 observations with SNR = 25; for other cases, 
# please change the variable 'obs_choice' and 'snr'
# simulate three different sampling frequencies for system
obs_41 = list(times_by=0.5, obs=41, basis_by=0.5, out_nquad=100)
obs_21 = list(times_by=1, obs=21, basis_by=1, out_nquad=50)
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

