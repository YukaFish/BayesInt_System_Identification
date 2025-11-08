# load required packages
library(ggplot2)
library(cowplot)
library(latex2exp)
library(xtable)
data(tcell)

# read the result
new_fit_lst <-Sys.glob(paste('real_eta*rds',sep=""))
time_fit_lst <- Sys.glob(paste('time_eta*rds',sep=""))

out.gene.45_lst <- list()
freq.gene.lst<- c()
sum.graph.lst <- c()
for (i in 1:length(new_fit_lst)) {
  new_fit<-readRDS(new_fit_lst[i])
  new_fit_sum <- data.frame(summary(new_fit)$summary)
  graph <- matrix(NA, nc=num_sigma, nr=num_sigma)
  theta <- matrix(get_posterior_mean(new_fit, par = 'theta'),ncol=num_sigma)
  graph[theta>0.5] <-T; graph[theta<0.5] <- F
  sum.graph.lst[i] <- sum(graph)
  freq.gene.lst <- c(freq.gene.lst, order(apply(graph,2,sum))[44:58])
  out.gene.45 <- which(graph[,45]==TRUE)
  out.gene.45_lst[[i]] = out.gene.45
}
saveRDS(out.gene.45_lst, file='out45.lst.rds')
saveRDS(freq.gene.lst, file='freq.lst.rds')
saveRDS(sum.graph.lst, file='sum.graph.rds')

sum.graph.lst <- readRDS('sum.graph.rds')
out.gene.45_lst <- readRDS('out45.lst.rds')
freq.gene.lst <- readRDS('freq.lst.rds')

# print the summary table for computing time
time_lst <- c()
for (i in 1:length(time_fit_lst)) {
  time_lst <- c(time_lst, readRDS(time_fit_lst[i]))
}
summary(time_lst)

# print the summary for selection of non-zero functions
summary(sum.graph.lst)

# print the frequency of top 15 genes with highest number of outward influence genes
freq.gene.lst_vec=c()
for (i in 1:length(freq.gene.lst)) {
  freq.gene.lst_vec <- c(freq.gene.lst_vec, freq.gene.lst[[i]])
}
freq.gene.lst_data <- as.data.frame(table(freq.gene.lst_vec))
freq.gene.lst_data <- freq.gene.lst_data[order(freq.gene.lst_data$Freq),]
print(xtable(freq.gene.lst_data),include.rownames=FALSE)

# print the frequency of each gene being an outward influence gene of FYB (gene 45)
out.gene.45_vec <- c()
for (i in 1:length(out.gene.45_lst)) {
    out.gene.45_vec <- c(out.gene.45_vec, out.gene.45_lst[[i]])
}
out.gene.45_data <- as.data.frame(table(out.gene.45_vec))
out.gene.45_data <- out.gene.45_data[order(out.gene.45_data$Freq),]
print(xtable(out.gene.45_data),include.rownames=FALSE)


