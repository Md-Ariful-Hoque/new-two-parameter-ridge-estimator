# Outliers
    num_out <- round(0.25 * n) # 25% outliers
    out_ind <- sample(1:n, num_out, replace = F)
    e[out_ind] <- quantile(e, 0.75) + 5 * IQR(e)
