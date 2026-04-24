# Real data analysis using all estimators

library(readxl)
library(car)

source("R/estimators_two_parameter.R")
source("R/mse_functions.R")

load_bodyfat <- function() {
  readxl::read_excel("data/raw/bodyfat.xlsx")
}

data <- load_bodyfat()
  
# CHANGE these names if your Excel column names are different
  y <- data$HDL
  X <- as.matrix(data[, c("Weight", "Height", "Waist", "Arm")])
# Correlation
cor(X)
  
xtx <- t(X) %*% X

# VIF
library(car)
model <- lm(y ~ x1 + x2 + x3 + x4)
vif(model)

# Condition Number
n <- nrow(x)
p <- ncol(x)
# Eigen values
l <- eigen(t(x) %*% x)$values

CN = sqrt(max(l)/min(l))

# Eigen Vector
w <- eigen(t(x) %*% x)$vectors
alphaht <- t(w) %*% (solve(t(x) %*% x) %*% t(x) %*% y)
# Identity matrix
I <- diag(p)

MSE_ols <- sigma2ols[1,1] * sum(1/l)
MSE_lte <- sigma2ols[1,1] * sum((d_lte - l)^2 / (l * (l + k_lte)^2)) + sum((d_lte + k_lte)^2 * alphaht^2 / (l + k_lte)^2)
MSE_tp1 <- sigma2ols[1,1] * sum((l + k_tp_am * d_tp1)^2 / (l * (l + k_tp_am)^2)) + sum(k_tp_am^2 * (d_tp1 - 1)^2 * alpha^2 / (l + k_tp_am)^2)
MSE_nbe1 <- sigma2ols[1,1] * sum(l * (l + (d_nbe1 + k_hk))^2 / ((l + 1)^2 * (l + k_hk)^2)) + sum(((1 - d_nbe1) * l + k_hk)^2 * alphaht^2 / ((l + k_hk)^2 * (l + 1)^2))
MSE_yc1 <- sigma2ols[1,1] * sum(l * (l + d_yc1)^2 / ((l + 1)^2 * (l + k_yc_am)^2)) + sum(((k_yc_am + 1 - d_yc1) * l + k_yc_am)^2 * alpha^2 / ((l + k_yc_am)^2 * (l + 1)^2))
MSE_mautp <- sigma2ols[1,1] * sum(l * (l + 2 * k_mautp - k_mautp * d_mautp)^2 * (l + k_mautp * d_mautp)^2 / (l + k_mautp)^3) +
  sum((l * (l + 2 * k_mautp - k_mautp * d_mautp) * (l + k_mautp * d_mautp) - (l + k_mautp)^3)^2 * alpha^2 / (l + k_mautp)^2)
MSE_mrt <- sigma2ols[1,1] *  sum(l/(l + k_mrt * (1 + d_mrt))^2) +
  sum(k_mrt^2 * (1 + d_mrt)^2 * alpha^2 / (l + k_mrt * (1 + d_mrt))^2)
MSE_dk <- sigma2ols[1,1] * sum((l - k_dk * (1 + d_dk))^2 / (l * (l + k_dk * (1 + d_dk))^2)) +
  sum(4 * k_dk^2 * (1 + d_dk)^2 * alpha^2 / (l * (l + k_dk + (1 + d_dk))^2))
MSE_gtp <- sigma2ols[1,1] * sum((l + k_gtp * d_gtp)^2 / (l * (l + k_gtp)^2)) +
  sum(k_gtp^2 * (d_gtp - 1)^2 * alpha^2 / (l + k_gtp)^2)
MSE_mlrt <- sigma2ols[1,1] * sum(l * (l + d_mlrt)^2 / ((l + k_mlrt * (1 + d_mlrt))^2 * (l + 1)^2)) +
  sum((l * (1 - d_mlrt) + k_mlrt * (1 + d_mlrt) * (l + 1))^2 * alpha^2 / ((l + k_mlrt * (1 + d_mlrt))^2 * (l + 1)^2))
MSE_nrt <- sigma2ols[1,1] * sum(l / (l + k_nrt + d_nrt)^2) + sum((k_nrt + d_nrt)^2 * alpha^2 / (l + k_nrt + d_nrt)^2)
MSE_lkl <- sigma2ols[1,1] * sum((l - k_lkl)^2 * (l + d_lkl)^2 / (l * (l + k_lkl)^2 * (l + 1)^2)) +
  sum(((l + d_lkl) * (l - k_lkl) - (l + k_lkl) * (l + 1))^2 * alpha^2 / ((l + k_lkl)^2 * (l + 1)^2))
MSE_tpr <- sigma2ols[1,1] * sum(q1[1,1]^2 * l / (l + k_tpr1)^2) +
  sum((q1[1,1] * l - l - k_tpr1)^2 * alpha^2 / (l + k_tpr1)^2)
MSE_new <- sigma2ols[1,1] * sum(q_new1^2 * l / (l + k_tpr1)^2) +
  sum((q_new1 * l - l - k_tpr1)^2 * alpha^2 / (l + k_tpr1)^2)













