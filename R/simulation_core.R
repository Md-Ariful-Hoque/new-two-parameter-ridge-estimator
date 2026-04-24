# Core simulation functions

source("R/data_generation.R")
source("R/estimators.R")
source("R/mse_functions.R")

one_replication <- function(n, p, rho, sigma) {
  x <- generate_x(n, p, rho)

  eig <- eigen(t(x) %*% x)
  beta <- eig$vectors[, which.max(eig$values)]
  beta <- beta / sqrt(sum(beta^2))

  y <- generate_y(x, beta, sigma)

  ols <- ols_estimator(x, y)

  k <- sum((eig$values / abs(t(eig$vectors) %*% ols))^2) *
    (sigma^2 / max((t(eig$vectors) %*% ols)^2))

  ntpr <- ntpr_estimator(x, y, k)

  data.frame(
    OLS = mse_estimator(ols, beta),
    NTPR = mse_estimator(ntpr, beta)
  )
}

run_simulation <- function(n, p, rho, sigma, size = 5000) {
  res <- replicate(
    size,
    one_replication(n, p, rho, sigma),
    simplify = FALSE
  )

  res <- do.call(rbind, res)
  colMeans(res)
}
