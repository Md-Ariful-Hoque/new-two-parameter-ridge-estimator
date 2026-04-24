# Outlier contamination functions

add_outliers <- function(y, proportion = 0.25) {
  n <- length(y)
  n_out <- floor(proportion * n)

  out_index <- sample(1:n, n_out)
  out_value <- quantile(y, 0.75) + 5 * IQR(y)

  y[out_index] <- out_value
  y
}
