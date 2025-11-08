# load required packages
library(ggplot2)
library(cowplot)
library(latex2exp)

# read the file list
new_fit_lst <- Sys.glob(paste('sim.result.', 
                              obs_choice$obs, 'obs.', 
                              snr, 'snr.*.rds', sep=""))

# define the dense basis matrix for performance metrics
dense_basismat <- eval.basis(seq(0,20,length=2001), bsbasis)

# define the true graph
true_matrix <- matrix(c(T,T,rep(F,8), T, T, rep(F, 8), rep(F,3), T, rep(F,6),
                        rep(F,2),T,rep(F,7), rep(F, 5), T,rep(F,4), rep(F,4), T,rep(F,5),
                        rep(F,40)),ncol=num_sigma)


# calculate the error of f
cal_f_error <- function(I0quadbasismat, new_fit){
  I0quadbasismat = eval.basis(seq(0,20,length = 2001), bsbasis, 0)
  ix <- I0quadbasismat%*%matrix(get_posterior_mean(new_fit, par = 'b'),ncol=num_sigma)
  for (i in 1:num_sigma) {
    assign(paste("bs_est.x", i, sep = ""), scale_bs(cal_basis(ix[,i],x_knots)))
  }
  a = array(get_posterior_mean(new_fit, par = 'a'),dim=c(dim(bs_est.x1)[2], num_sigma, num_sigma))
  est_f_result <- array(NA, dim=c(num_sigma, num_sigma, dim(dense_basismat)[1]))
  for (i in 1:num_sigma) {
    for (j in 1:num_sigma) {
      est_f_result[i,j,] <- get(paste("bs_est.x", i, sep = ""))%*%(a[,j,i])
    }
  }
  f_error_matrix <- matrix(NA, nc=num_sigma, nr=num_sigma)
  true_f_xdense <- true_f_matrix(cbind(seq(0,20,length=2001), xdense))
  for (i in 1:num_sigma) {
    for (j in 1:num_sigma) {
      f_error_matrix[i,j] <- mean((true_f_xdense[i,j,]-scale(est_f_result[i,j,],scale=F))^2)
    }
  }
  f_error_matrix
}

# calculate tpr and fpr
cal_pr <- function(new_fit) {
  graph <- matrix(NA, nc=num_sigma, nr=num_sigma)
  theta <- matrix(get_posterior_mean(new_fit, par = 'theta'),ncol=num_sigma)
  graph[theta>0.5] <-T; graph[theta<0.5] <- F
  list(tpr = sum(graph[true_matrix])/sum(true_matrix), 
       fpr = sum(graph[!true_matrix])/(num_sigma*num_sigma-sum(true_matrix)))
}

# calculate the error of x
cal_x_error <- function(new_fit){
  dense_basismat = eval.basis(seq(0,20,length=2001), bsbasis, 0)
  est_x <- dense_basismat%*%matrix(get_posterior_mean(new_fit, par = 'b'),ncol=num_sigma)
  mean((est_x-xdense)^2)
}

# calculate the coverage probabilities for trajectories
cov_prob_traj_f <- function(b_iter, a_iter, a0_iter) {
  est_f_lst <- array(NA, dim = c(dim(b_iter)[1], out_nquad, num_sigma))
  for (m in 1:(dim(b_iter)[1])) {
    ix <- I1quadbasismat%*%matrix(b_iter[m,],ncol=num_sigma)
    norm_ix <- normalize(ix)
    ix <- 1/(1+exp(-norm_ix))
    
    for (x_index in 1:num_sigma) {
      for (i in 1:num_sigma) {
        assign(paste("bs_ix", i, sep = ""), integral_bs(cal_basis(ix[,i],x_knots)))
      }
      est_f <- a0_iter[m,x_index]*scale(out_quadpts,scale=F)
      for (i in 1:num_sigma) {
        est_f <- est_f + get(paste("bs_ix", i, sep = ""))%*%(a_iter[m,,x_index,i])
      }
      est_f_lst[m,,x_index] <- est_f
    }
  }
  result_0.975 <- (apply(est_f_lst, c(2,3), quantile, probs = 0.975))
  result_0.025 <- (apply(est_f_lst, c(2,3), quantile, probs = 0.025))
  cov_prob <- matrix(NA, nr = dim(result_0.025)[1], nc = num_sigma)
  xsim <- data.frame(out_quadpts)
  xtrue <- cbind(seq(0,20,length=2001), xdense)
  xtrue <- data.frame(xtrue)
  names(xtrue)[names(xtrue) == 'V1'] <- 'time'
  xtrueFunc <- lapply(2:(ncol(xdense)+1), function(j)
    approxfun(xtrue[, "time"], xtrue[, j]))
  xsim <- cbind(xsim, sapply(xtrueFunc, function(f) f(xsim$out_quadpts)))
  for (mm in 1:num_sigma) {
    cov_prob[,mm] <- ((result_0.975[,mm] > scale(xsim[,mm+1],scale=F)) & (result_0.025[,mm] < scale(xsim[,mm+1],scale=F)))
  }
  cov_prob
}

# calculate the performance metric
tpr = fpr =  c()
cov_prob_traj_lst = array(NA, dim=c(length(new_fit_lst),out_nquad,num_sigma))
f_true_lst = f_false_lst = f_sum_lst = traj_error_lst = c()
for (new_fit_i in 1:length(new_fit_lst)) {
  # read the result from fit
  new_fit <-readRDS(new_fit_lst[new_fit_i]) 
  new_fit_sum <- data.frame(summary(new_fit)$summary)
  a0_iter <- rstan::extract(new_fit, pars = 'a0')[[1]]
  a_iter <- array(rstan::extract(new_fit, pars = 'a')[[1]],dim=c(dim(a0_iter)[1], x_knots+3, num_sigma, num_sigma))
  b_iter <- rstan::extract(new_fit, pars = "b")[[1]]
  
  # calculate tpr and fpr
  pr <- cal_pr(new_fit)
  tpr <- c(tpr, pr$tpr)
  fpr <- c(fpr, pr$fpr)
  
  # calculate the coverage probabilities for trajectories
  cov_prob_traj_lst[new_fit_i,,] <- cov_prob_traj_f(b_iter, a_iter, a0_iter)
  
  # calculate the MSE for f
  f_error_matrix <- cal_f_error(I0quadbasismat, new_fit)
  f_true_lst <- c(f_true_lst, mean(f_error_matrix[true_matrix]))
  f_false_lst <- c(f_false_lst, mean(f_error_matrix[!true_matrix]))
  f_sum_lst <- c(f_sum_lst, mean(f_error_matrix))
  
  # calculate the MSE for trajectories
  traj_error_lst <- c(traj_error_lst, cal_x_error(new_fit))
  
}

# save the results for future use
result <- list(cov_prob_traj_lst = cov_prob_traj_lst, 
               f_true_lst = f_true_lst,f_false_lst = f_false_lst, 
               f_sum_lst = f_sum_lst,
               traj_error_lst = traj_error_lst,
               tpr=tpr, fpr=fpr)
saveRDS(result, file = paste('results.obs.', obs_choice$obs, '.snr.', snr, '.rds', sep=''))

# read the results
perf <- readRDS(paste('results.obs.', obs_choice$obs, '.snr.', snr, '.rds', sep=''))
summary(perf$tpr)
sd(perf$tpr)
summary(perf$fpr)
sd(perf$fpr)
summary(perf$f_true_lst)
sd(perf$f_true_lst)
summary(perf$f_false_lst)
sd(perf$f_false_lst)
summary(perf$f_sum_lst)
sd(perf$f_sum_lst)
summary(perf$traj_error_lst)
sd(perf$traj_error_lst)

cov_traj_lst <- perf$cov_prob_traj_lst
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
dev.new(width=10, height=5, noRStudioGD = TRUE)
p_lst <- plot_cov_traj_prob(out_quadpts, cov_traj_lst)
plot_grid(plotlist=p_lst, ncol=5)

