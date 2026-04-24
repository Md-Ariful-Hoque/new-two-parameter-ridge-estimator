# Create tables and plots

library(ggplot2)

# Load results
res <- read.csv("../results/tables/main_simulation.csv")

# Plot example
ggplot(res, aes(x = n, y = NTPR, color = as.factor(rho))) +
  geom_line() +
  facet_grid(p ~ sigma) +
  theme_minimal()

# Save plot
ggsave("../results/figures/mse_plot.png")
