# Data generation functions

generate_x <- function(n, p, rho) {
  z <- matrix(rnorm(n * (p + 1)), nrow = n, ncol = p + 1)
  x <- matrix(0, nrow = n, ncol = p)

  for (j in 1:p) {
    x[, j] <- sqrt(1 - rho^2) * z[, j] + rho * z[, p + 1]
  }

  scale(x)
}

generate_y <- function(x, beta, sigma) {
  as.vector(x %*% beta + rnorm(nrow(x), mean = 0, sd = sigma))
}
