# Simulation with outliers

library(MASS)

source("../R/simulation_core.R")
source("../R/outlier_functions.R")

run_outlier_sim <- function(n, p, rho, sigma, size = 5000) {

  res <- replicate(size, {

    x <- generate_x(n, p, rho)

    eig <- eigen(t(x) %*% x)
    beta <- eig$vectors[, which.max(eig$values)]
    beta <- beta / sqrt(sum(beta^2))

    y <- generate_y(x, beta, sigma)
    y <- add_outliers(y)

    ols <- ols_estimator(x, y)

    k <- 0.1
    ntpr <- ntpr_estimator(x, y, k)

    c(
      OLS = mse_estimator(ols, beta),
      NTPR = mse_estimator(ntpr, beta)
    )

  })

  colMeans(t(res))
}

# Example run
res_out <- run_outlier_sim(100, 5, 0.90, 5)

write.csv(res_out, "../results/tables/outlier_results.csv")
print(res_out)
