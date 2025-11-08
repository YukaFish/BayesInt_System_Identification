results <- read.table(paste0("results/gaussian_a_N", N-1, 'R', R, 'NT', n_tol,
                             'NW', NUM_REWEIGHT, 'SNR', SNR, '.txt', sep=''),
                      header=T,sep = ",", dec = ".")

# read the results from JADE
jade <- results[results$method == 'jointest_rw4', ]
summary(jade$ftheta.mse.f)
sd(jade$ftheta.mse.f)
summary(jade$theta.mse)
sd(jade$theta.mse)
summary(jade$tpr)
sd(jade$tpr)
summary(jade$fpr)
sd(jade$fpr)


# read the results for GRADE
grade <- results[results$method == 'GRADE', ]
summary(grade$ftheta.mse.f)
sd(grade$ftheta.mse.f)
summary(grade$theta.mse)
sd(grade$theta.mse)
summary(grade$tpr)
sd(grade$tpr)
summary(grade$fpr)
sd(grade$fpr)
