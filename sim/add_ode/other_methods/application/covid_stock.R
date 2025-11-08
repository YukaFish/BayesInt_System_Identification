library(tidyverse)
library(lubridate)

source("R/bcd_additive_mixed.r")

TRANSFORMER <- "sigmoid"
transformer <- pracma::sigmoid

topic_stock_data <- readRDS(
  "./application/data/topic_stock_data.rds"
)


for (i in 1:10) {
  start_date <- ymd(paste0("2020-", 1, "-01"), tz = "GMT") +
    months(2 * (i - 1))
  end_date <- start_date + months(2) - days(1)
  
  idx <- which(
    topic_stock_data$dates < end_date &
      topic_stock_data$dates >= start_date
  )
  
  lam_gam_list <- c(
    10^seq(0, 2, 0.25),
    10^seq(-3, -0.5, 0.5),
    10^(-4)
  )
  
  for (i in 1:3) {
    res <- tryCatch(
      bcd.ode.additive.mixed(
        list(topic_stock_data$ydata[idx, ]), seq_along(idx),
        family = topic_stock_data$family,
        dist_args = NULL,
        lam_th_min = 1e0, lam_th_max = 1e0,
        delta_update_method = "decrease",
        lam_gam_list = lam_gam_list,
        lam_gam_ratio = "auto",
        delta_lam_gam = NULL,
        knots.theta.dt = 1,
        dt.int = 0.1,
        scaling_factor = 1,
        norder_f = 4, nknots_f = 6,
        init_method = "baseline",
        lam_gam_init = 10^seq(-4, 1, 0.1),
        penalty = "grpLasso",
        penalty.weighted = TRUE,
        smooth_package = "stats",
        theta_method = "1dH",
        initialize_penalty_weights = FALSE,
        MAXITER = 4,
        NUM_REWEIGHT = 2,
        random_order = TRUE,
        BICobj = "ode",
        eps.conv = 1e-4,
        ada.gam = 1,
        L_lower_bound_factor = 1e-8,
        know_truth = FALSE,
        fit.trajectory = FALSE,
        mode.test = FALSE
      ),
      error = function(e) {
        return(NULL)
      }
    )
    
    if (is.null(res)) {
      # lam_gam_list <- lam_gam_list * 10
      next
    }
    
    # num_nonzero <- sum(apply(res$graph.all, 3, any), na.rm = TRUE)
    # if (num_nonzero < 0.5 * length(lam_gam_list)) {
    #   lam_gam_list <- lam_gam_list * 0.5
    # }
    
    break
  }
  
  saveRDS(
    res,
    paste0(
      "./application/result/res_",
      start_date, "_", end_date, ".rds"
    )
  )
}
