# MSE calculation

mse_estimator <- function(beta_hat, beta_true) {
  sum((beta_hat - beta_true)^2)
}
