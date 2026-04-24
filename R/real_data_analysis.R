# Real data analysis using all estimators

library(readxl)
library(car)

source("R/estimators_two_parameter.R")
source("R/mse_functions.R")

load_bodyfat <- function() {
  readxl::read_excel("data/raw/bodyfat.xlsx")
}

run_real_data_estimators <- function() {
  
  data <- load_bodyfat()
  
  # CHANGE these names if your Excel column names are different
  y <- data$HDL
  X <- as.matrix(data[, c("Weight", "Height", "Waist", "Arm")])
  
  X <- scale(X)
  y <- scale(y)
  
  xtx <- t(X) %*% X
  ols <- ols_estimator(X, y)
  
  eig <- eigen(xtx)
  lambda <- eig$values
  alpha_hat <- t(eig$vectors) %*% ols
  sigma2 <- as.numeric(t(y - X %*% ols) %*% (y - X %*% ols) / (nrow(X) - ncol(X) - 1))
  
  # Simple starting parameters
  k <- sigma2 / sum(alpha_hat^2)
  d <- min(sigma2 / alpha_hat^2)
  d_vec <- rep(d, ncol(X))
  q <- 1
  
  estimates <- list(
    OLS  = ols,
    LTE  = lte_estimator(xtx, ols, k, d),
    TP   = tp_estimator(xtx, ols, k, d),
    NBE  = nbe_estimator(xtx, ols, k, d),
    YC   = yc_estimator(xtx, ols, k, d),
    MRT  = mrt_estimator(xtx, ols, k, d),
    DK   = dk_estimator(xtx, ols, k, d),
    GTP  = gtp_estimator(xtx, ols, k, d_vec),
    LKL  = lkl_estimator(xtx, ols, k, d),
    NTP  = ntp_estimator(xtx, ols, k, d),
    NRT  = nrt_estimator(xtx, ols, k, d),
    TPR  = tpr_estimator(xtx, ols, k, q),
    NTPR = ntpr_estimator(xtx, ols, k)
  )
  
  beta_ols <- as.numeric(ols)
  
  result_table <- do.call(rbind, lapply(names(estimates), function(name) {
    beta_hat <- as.numeric(estimates[[name]])
    mse <- sum((beta_hat - beta_ols)^2)
    
    c(
      Estimator = name,
      X1 = beta_hat[1],
      X2 = beta_hat[2],
      X3 = beta_hat[3],
      X4 = beta_hat[4],
      MSE = mse
    )
  }))
  
  result_table <- as.data.frame(result_table)
  
  write.csv(result_table, "results/tables/real_data_estimators.csv", row.names = FALSE)
  
  return(result_table)
}
