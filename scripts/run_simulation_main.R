# Main simulation script (clean data)

library(MASS)

source("../R/simulation_core.R")

# Settings
n_vals   <- c(30, 50, 100, 200)
p_vals   <- c(3, 5, 10)
rho_vals <- c(0.80, 0.90, 0.95, 0.99)
sigma_vals <- c(1, 5, 10)

results_list <- list()
counter <- 1

for (n in n_vals) {
  for (p in p_vals) {
    for (rho in rho_vals) {
      for (sigma in sigma_vals) {

        cat("Running:", n, p, rho, sigma, "\n")

        res <- run_simulation(n, p, rho, sigma)

        results_list[[counter]] <- data.frame(
          n = n,
          p = p,
          rho = rho,
          sigma = sigma,
          t(res)
        )

        counter <- counter + 1
      }
    }
  }
}

final_results <- do.call(rbind, results_list)

# Save
write.csv(final_results, "../results/tables/main_simulation.csv", row.names = FALSE)
saveRDS(final_results, "../results/saved_objects/main_simulation.rds")

print(final_results)
