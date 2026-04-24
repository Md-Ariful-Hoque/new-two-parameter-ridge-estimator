# Data generation functions

z <- matrix(rnorm(n * p), nrow = n, ncol = p)
# Generate x
x <- sqrt(1 - rho^2) * z + rho * z[,p]

X <- t(x) %*% x
l <- eigen(X)$values
del <- diag(l)
w <- eigen(X)$vectors
b <- matrix(w[,1])

# Identity matrix
I <- diag(p)

# error term
e <- matrix(rnorm(n, mean = 0, sd = sqrt(sigma)))

# canonical form
Z <- x %*% w

alpha <- t(w) %*% b

# Generate y
y <- Z %*% alpha + e
