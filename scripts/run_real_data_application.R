library(readxl)

source("../R/data_loading.R")
source("../R/estimators_two_parameter.R")

data <- load_bodyfat()

# Adjust based on your dataset
y <- data[[1]]
X <- as.matrix(data[, -1])

ols <- ols_estimator(X, y)
ntpr <- ntpr_estimator(t(X) %*% X, ols, k = 0.1)

print(ols)
print(ntpr)
