# Real data application

library(car)

source("../R/estimators.R")

# Load your dataset (put in data/raw/)
data <- read.csv("../data/raw/bodyfat.csv")

# Example model
y <- data$y   # change to your response
X <- as.matrix(data[, -1])  # predictors

# OLS
ols <- solve(t(X) %*% X) %*% t(X) %*% y

# Ridge parameter
k <- 0.1

ntpr <- ntpr_estimator(X, y, k)

# Diagnostics
vif_values <- vif(lm(y ~ X))
cond_number <- kappa(X)

results <- list(
  OLS = ols,
  NTPR = ntpr,
  VIF = vif_values,
  ConditionNumber = cond_number
)

saveRDS(results, "../results/saved_objects/real_data.rds")
print(results)
