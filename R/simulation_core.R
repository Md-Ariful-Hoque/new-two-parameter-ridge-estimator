# Core simulation functions

source("R/data_generation.R")
source("R/estimators.R")
source("R/mse_functions.R")

one_replication <- function(n, p, rho, sigma) {
  x 

  eig <- eigen(t(x) %*% x)
  w <- eigen(X)$vectors
  # True beta
  b <- matrix(w[,1]) 

  Z <- x %*% w
  alpha <- t(w) %*% b
  
  y 

  ols 
  mse_ols[i] <- t(ols - alpha) %*% (ols - alpha)

 
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
