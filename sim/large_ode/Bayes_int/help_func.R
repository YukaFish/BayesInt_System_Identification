# normalize the data matrix 
normalize <- function(theta_mat) {
  center_th <- 0.5 * (apply(theta_mat, 2, max) + apply(theta_mat, 2, min))
  scale_th <- ((apply(theta_mat, 2, max) - apply(theta_mat, 2, min)))
  theta_mat <- as.matrix(theta_mat)
  theta_mat <- sweep(theta_mat, 2, center_th, FUN = "-")
  theta_mat <- sweep(theta_mat, 2, scale_th, FUN = "/")
  return(theta_mat)
}

# get the percentile of x according to the position
get_percentile<-function(x,nspecs){
  nknots<-nspecs[2]
  norder<-nspecs[1]
  if(nknots>0){
    out<-seq(from=0, to=1, length.out=(nknots+2) )
    out<-quantile(x,probs=out )
  } else{out=NULL }
  return(list(percentile=out,norder=norder))
}

# calculate the basis matrix for x
cal_basis <- function(ix, x_knots){
  specs <- get_percentile(ix, c(3,x_knots))
  knts<-specs$percentile[-c(1,length(specs$percentile))]
  bs_ix<-bs(ix,degree=3,knots=knts)
  bs_ix
}

# get the integral of basis matrix
integral_bs <- function(bs_x) {
  int_bs_x =  temp1 = c()
  for (i in 1:dim(bs_x)[2]) {
    temp1 <- cbind(temp1, in_quadwts_vec*bs_x[,i])
  }
  temp2 <- apply(temp1, 2, cumsum)
  for (i in 1:out_nquad) {
    int_bs_x <- rbind(int_bs_x, temp2[i*in_nquad,])
  }
  scale(int_bs_x,center=T,scale=F)
}

# scale the basis functions
scale_bs <- function(bs) {
  bs <- scale(bs,center=T, scale=F)
  bs
}

# help function: get the BIC of regression
get.BIC <- function(nlambda, N, response, preds, dfs, fam = NULL) {
  BICs <- sapply(
    seq_len(nlambda),
    function(i) {
      MSE <- mean((response - preds[, i])^2)
      bic <- N * log(MSE) + dfs[i] * log(N)
      return(bic)
    }
  )
}

# give the initial estimates for coefficients of additive components
est_coef <- function(y, x) {
  index <- c(0,rep(c(1:(num_sigma)),each = (ncol(x)-1)/num_sigma))
  fit <- grpreg::grpreg(x, y, group=index, lambda=10^seq(-2, 2, 0.1), 
                        penalty="grLasso", family="gaussian")
  lambdas <- fit$lambda
  nlambda <- length(lambdas)
  preds <-predict(fit, x)
  dfs <- fit$df
  BICs <- get.BIC(nlambda, times, y, preds, dfs)
  minbic.id <- which.min(BICs)
  beta <- matrix(as.vector(fit$beta[,minbic.id])[-c(1:2)],nc=num_sigma)
  graph <- apply(beta != 0, 2, any)
  beta[beta==0] <- rnorm(length(beta[beta==0]),0,1e-10)
  a0 <- as.vector(fit$beta[,minbic.id])[2]
  return(list(beta = beta,a0 = a0, graph = graph))
}

# true functions
true_f_matrix <- function(xsim) {
  true_f <- array(0, dim=c(num_sigma, num_sigma, dim(xsim)[1]))
  true_f[1,1,] <- scale(1.2*xsim[,2]+0.3*xsim[,2]^2-0.6*xsim[,2]^3,scale=F)
  true_f[1,2,] <- scale(-2*xsim[,2]+0.4*xsim[,2]^3,scale=F)
  true_f[2,1,] <- scale(0.1*xsim[,3]+ 0.2*xsim[,3]^2 + 0.2*xsim[,3]^3,scale=F)
  true_f[2,2,] <- scale(0.5*xsim[,3]+0.2*xsim[,3]^2-0.3*xsim[,3]^3,scale=F)
  true_f[4,3,] <- scale(0.3*xsim[,5]+ 0.4*xsim[,5]^2+ 0.1*xsim[,5]^3,scale=F)
  true_f[3,4,] <- scale(0.2*xsim[,4]-0.1*xsim[,4]^2 -0.2*xsim[,4]^3,scale=F)
  true_f[6,5,] <- scale(0.1*xsim[,7]-0.8*xsim[,7]^3,scale=F)
  true_f[5,6,] <- scale(0.5*xsim[,6]^3,scale=F)
  true_f
}