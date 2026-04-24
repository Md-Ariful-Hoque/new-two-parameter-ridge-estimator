# Estimator functions

ols_estimator <- function(x, y) {
  solve(t(x) %*% x) %*% t(x) %*% y
}

ridge_estimator <- function(x, y, k) {
  p <- ncol(x)
  solve(t(x) %*% x + k * diag(p)) %*% t(x) %*% y
}

ntpr_estimator <- function(x, y, k) {
  xtx <- t(x) %*% x
  p <- ncol(x)

  ols <- ols_estimator(x, y)

  q_min <- min(eigen((xtx + k * diag(p)) %*% solve(xtx),
                     symmetric = TRUE,
                     only.values = TRUE)$values)

  q_min * solve(xtx + k * diag(p)) %*% xtx %*% ols
}
