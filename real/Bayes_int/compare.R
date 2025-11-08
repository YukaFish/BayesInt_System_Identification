# load required packages
library(ggplot2)
library(cowplot)
library(latex2exp)
library(xtable)
library(igraph)
data(tcell)

# read the result
new_fit <- readRDS(Sys.glob(paste('real.*rds',sep="")))
runtime <- readRDS(Sys.glob(paste('time.*rds',sep="")))
new_fit_sum <- data.frame(summary(new_fit)$summary)
graph <- matrix(NA, nc=num_sigma, nr=num_sigma)
theta <- matrix(get_posterior_mean(new_fit, par = 'theta'),ncol=num_sigma)
graph[theta>0.5] <-T; graph[theta<0.5] <- F
sum(graph)
order(apply(graph,2,sum))
out.gene.45 <- which(graph[,45]==TRUE)

# plot the trajectories with credible interval
dense_basismat = eval.basis(seq(0,72,length=2001), bsbasis, 0)
est.x <- dense_basismat%*%matrix(get_posterior_mean(new_fit, par = 'b'),ncol=num_sigma)
est.x.97.5 <- dense_basismat %*% matrix(apply(rstan::extract(new_fit, pars = "b")[[1]], 2, quantile, probs=0.975), ncol = num_sigma)
est.x.2.5 <- dense_basismat %*% matrix(apply(rstan::extract(new_fit, pars = "b")[[1]], 2, quantile, probs=0.025), ncol = num_sigma)
dense_times <- seq(0,72,length=2001)
colors <- c("observations" = "black", "estimated trajectories" = "red")
plot_traj_lst <- c(45,46,49,55)
p_lst <- list()
for (index in 1:length(plot_traj_lst)) {
  obs <- c()
  for (m in 1:length(times)) {
    obs <- rbind(obs, expand.grid(times[m],as.numeric(y[m,plot_traj_lst[index],])))
  }
  p_df <- data.frame(times = dense_times,freq = est.x[, plot_traj_lst[index]],
                     low=est.x.2.5[, plot_traj_lst[index]],
                     upp=est.x.97.5[, plot_traj_lst[index]])
  p <- ggplot() +
    geom_point(data=obs,aes(y=Var2,x=Var1,color='observations'),size=2,shape=17) +
    geom_ribbon(data=p_df,aes(ymin = low, 
                              ymax = upp, x=dense_times),
                colour = NA, fill = "red", alpha = 0.2) +
    geom_point(data=p_df,aes(y=freq,x=times,color="estimated trajectories"),size=0.6) +
    geom_line(data=p_df,aes(y=freq,x=times,color="estimated trajectories"),size=0.6) +
    scale_color_manual(values = colors)+
    theme(legend.position="none", axis.title = element_text(size = 15),
          axis.text = element_text(size = 13),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = "grey90"))+
    xlab('time (hours)') + ylab('gene expression') + 
    ggtitle(paste('Gene ',plot_traj_lst[index],sep=""))
  p_lst[[index]] = p
}
dev.new(width=8, height=8, noRStudioGD = TRUE)
plot_grid(plotlist=p_lst, ncol=2)


# plot the estimated regulation functions of genes influenced by gene 45 (FYB)
for (i in 1:num_sigma) {
  assign(paste("bs_est.x", i, sep = ""), scale_bs(cal_basis(est.x[,i],x_knots)))
}
a = array(get_posterior_mean(new_fit, par = 'a'),dim=c(num_sigma,num_sigma,dim(bs_est.x1)[2]))

est_f_result <- array(NA, dim=c(num_sigma, num_sigma, dim(dense_basismat)[1]))
for (i in 1:num_sigma) {
  for (j in 1:num_sigma) {
    est_f_result[i,j,] <- get(paste("bs_est.x", i, sep = ""))%*%(a[i,j,])
  }
}

p_lst <- list()
for (index in 1:length(out.gene.45)) {
  p_df <- data.frame(gene = est.x[,45],func = est_f_result[45,out.gene.45[index],])
  p <- ggplot() +
    geom_line(data=p_df,aes(y=func,x=gene,color="black"),size=0.6) +
    theme(legend.position="none", axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = "grey90"),
          plot.title = element_text(hjust = 1))+
    xlab('') + ylab('') + 
    ggtitle(paste('Gene',out.gene.45[index], 'By FYB', sep=" "))
  p_lst[[index]] = p
}
dev.new(width=10, height=5, noRStudioGD = TRUE)
plot_grid(plotlist=p_lst, ncol=5)+
  draw_label("Gene Expression", x=0.5, y=  0, vjust=-0.5, angle= 0) +
  draw_label("Function", x=  0, y=0.5, vjust= 1.5, angle=90)

# generate the variable selection result table
out_lst = in_lst = c()
for (i in 1:num_sigma) {
  out_lst <- c(out_lst, paste(which(graph[,i]==TRUE), collapse=", "))
  in_lst <- c(in_lst, paste(which(graph[i,]==TRUE), collapse=", "))
}
graph_df <- data.frame('Gene' = 1:58, 
                       'Inward influence genes' = in_lst,
                       'Outward influence genes' = out_lst)
row.names(graph_df) = NULL
graph_latex <- xtable(graph_df, label ='realresult',caption ='')
print(graph_latex, floating = FALSE, include.rownames=FALSE)

adjm <- matrix(sample(0:1, 100, replace = TRUE, prob = c(0.9, 0.1)), ncol = 10)
g1 <- graph_from_adjacency_matrix(t(graph))
colrs <- c("red", "gray50")
V(g1)$color <- colrs[2]
V(g1)$color[45] <- colrs[1]
edge.start <- ends(g1, es=E(g1), names=F)[,1]
edge.col <- V(g1)$color[edge.start]
edge.size <- ifelse(V(g1)$color[edge.start]=="red", 2,1)
colrs <- c("pink", "white")
V(g1)$color <- colrs[2]
V(g1)$color[45] <- colrs[1]
dev.new()
plot(g1,edge.arrow.size=.5,arrow.width=.3,vertex.size=9,
     vertex.label.color	= "black",
     edge.width = edge.size, margin=0,
     edge.color=edge.col,vertex.label.cex = 0.8,layout=layout.circle)
