# Two-parameter estimator functions

ols_estimator <- function(x, y) {
  solve(t(x) %*% x) %*% t(x) %*% y
}

lte_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  solve(xtx + k * diag(p)) %*% (xtx - d * diag(p)) %*% ols
}

tp_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  solve(xtx + k * diag(p)) %*% (xtx + k * d * diag(p)) %*% ols
}

nbe_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  solve(xtx + diag(p)) %*%
    (xtx + (d + k) * diag(p)) %*%
    solve(xtx + k * diag(p)) %*%
    xtx %*% ols
}

yc_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  solve(xtx + diag(p)) %*%
    (xtx + d * diag(p)) %*%
    solve(xtx + k * diag(p)) %*%
    xtx %*% ols
}

mrt_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  solve(xtx + k * (1 + d) * diag(p)) %*% xtx %*% ols
}

dk_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  solve(xtx + k * (1 + d) * diag(p)) %*%
    (xtx - k * (1 + d) * diag(p)) %*% ols
}

gtp_estimator <- function(xtx, ols, k, d_vec) {
  p <- ncol(xtx)
  D <- diag(d_vec)
  solve(xtx + k * diag(p)) %*% (xtx + k * D) %*% ols
}

lkl_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  C <- solve(xtx + diag(p)) %*% (xtx + d * diag(p))
  A <- solve(xtx + k * diag(p)) %*% (xtx - k * diag(p))
  C %*% A %*% ols
}

ntp_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  solve(xtx + k * d * diag(p)) %*%
    (xtx - k * d * diag(p)) %*% ols
}

nrt_estimator <- function(xtx, ols, k, d) {
  p <- ncol(xtx)
  solve(xtx + (k + d) * diag(p)) %*% xtx %*% ols
}

tpr_estimator <- function(xtx, ols, k, q) {
  p <- ncol(xtx)
  q * solve(xtx + k * diag(p)) %*% xtx %*% ols
}

ntpr_estimator <- function(xtx, ols, k) {
  p <- ncol(xtx)

  q_min <- min(eigen((xtx + k * diag(p)) %*% solve(xtx),
                     symmetric = TRUE,
                     only.values = TRUE)$values)

  q_min * solve(xtx + k * diag(p)) %*% xtx %*% ols
}
