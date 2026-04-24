# Two-parameter estimator functions

ols <- solve(del) %*% t(Z) %*% y

# Estimated alpha
alphaht <- t(w) %*% ols

sigma2ols <- t(y - x %*% ols) %*% (y - x %*% ols)/(n - p - 1)

# 1 Liu Type Two Parameter Estimator
    k_lte <- (l[1] - 100 * l[p]) / 99
    if (k_lte <= 0) {
      k_lte <- abs(k_lte)
    }
    d_lte <- sum((sigma2ols[1,1] - k_lte * alphaht^2) / (l + k_lte)^2) /
      sum((sigma2ols[1,1] + l * alphaht^2) / (l * (l + k_lte)^2))
    if (d_lte <= 0) {
      d_lte <- abs(d_lte)
    }
    A_lte <- solve(del + k_lte * I) %*% (del - d_lte * I)
    lte <- A_lte %*% ols

# 1 Liu Type Two Parameter Estimator
    k_lte <- (l[1] - 100 * l[p]) / 99
    if (k_lte <= 0) {
      k_lte <- abs(k_lte)
    }
    d_lte <- sum((sigma2ols[1,1] - k_lte * alphaht^2) / (l + k_lte)^2) /
      sum((sigma2ols[1,1] + l * alphaht^2) / (l * (l + k_lte)^2))
    if (d_lte <= 0) {
      d_lte <- abs(d_lte)
    }
    A_lte <- solve(del + k_lte * I) %*% (del - d_lte * I)
    lte <- A_lte %*% ols

# 2. Two Parameter Estimator
    d_hat_tp <- min(alphaht^2 / ((sigma2ols[1,1] / l) + alphaht^2)) - 0.01
    
    k_tp_am <- (1 / p) * sigma2ols[1,1] / sum(alphaht^2 - d_hat_tp * ((sigma2ols[1,1] / l) + alphaht^2))
    if (k_tp_am <= 0) {
      k_tp_am <- abs(k_tp_am)
    }
    
    d_tp1 <- sum((k_tp_am * alphaht^2 - sigma2ols[1,1]) / (l + k_tp_am)^2) /
      sum(k_tp_am * (sigma2ols[1,1] + alphaht^2 * l) / (l * (l + k_tp_am)^2))
    if (d_tp1 <= 0 || d_tp1 >= 1) {
      d_tp1 <- d_hat_tp
    }
    
    Atp1 <- solve(del + k_tp_am * I) %*% (del + k_tp_am * d_tp1 * I)
    tp1 <- Atp1 %*% ols

# 3. New Biased Estimator Based on Ridge
    
    k_hk <- sigma2ols[1,1] / sum(alphaht^2)
    if (k_hk <= 0) {
      k_hk <- abs(k_hk)
    }
    
    d_nbe1 <- sum(l * (alphaht^2 - sigma2ols[1,1]) / ((l + 1)^2 * (l + k_hk))) /
      sum(l * (l * alphaht^2 + sigma2ols[1,1]) / ((l + 1)^2 * (l + k_hk)^2))
    
    Anbe1 <- solve(del + I) %*% (del + (d_nbe1 + k_hk) * I) %*% solve(del + k_hk * I) %*% (del)
    nbe1 <- Anbe1 %*% ols

# 4. Yang and Chang
    d_min_yc <- min(sigma2ols[1,1] / alphaht^2)
    k_yc <- (sigma2ols[1,1] * (l + d_min_yc) - (1 - d_min_yc) * l * alphaht^2) /
      ((l + 1) * alphaht^2)
    
    # mean
    k_yc_am <- mean(k_yc) 
    if (k_yc_am <= 0) {
      k_yc_am <- abs(k_yc_am)
    }
    
    d_yc1 <-  sum((((k_yc_am + 1) * l + k_yc_am) * l * alphaht^2 - l^2 * sigma2ols[1,1]) / ( (l + 1)^2 * (l + k_yc_am)^2 )) /
      sum((sigma2ols[1,1] + l * alphaht^2) * l / ( (l + 1)^2 * (l + k_yc_am)^2 ))
    if (d_yc1 <= 0 || d_yc1 >= 1) {
      d_yc1 <- d_min_yc
    }
    Ayc1 <- solve(del + I) %*% (del + d_yc1 * I) %*% solve(del + k_yc_am * I) %*% (del)
    yc1 <- Ayc1 %*% ols

# 5. Modified Almost Unbiased Liu Estimator
    d_hat_maule <- min(sigma2ols[1,1]/alphaht^2)
    k_hat_maule <- (sigma2ols[1,1] * (l + 1)^2 * l - (1 - d_hat_maule)^2 * l * (sigma2ols[1,1] + alphaht^2)) /
      (alphaht^2 * (l + 1)^2)
    k_maule <- mean(k_hat_maule)
    if (k_maule <= 0) {
      k_maule <- abs(k_maule)
    }
    d_maule <- 1 - sqrt(abs(sum(l * (sigma2ols[1,1] - k_maule * alphaht^2) / ( (l + k_maule)^2 * (l + 1)^2 ))) /
                          sum(l * (sigma2ols[1,1] + l * alphaht^2) / ( (l + 1)^4 * (l + k_maule)^2 )))
    
    if (d_maule <= 0 || d_maule >= 1) {
      d_maule <- d_hat_maule
    }
    Amaule <- (I - (1 - d_maule)^2 * solve(del + I) %*% solve(del + I)) %*% solve(del + k_maule * I) %*% (del)
    maule <- Amaule %*% ols

# 6. Modified Almost Unbiased Two Parameter Estimator
    k_mautp <- p * sigma2ols[1,1] / sum(alphaht^2) 
    if (k_mautp <= 0) {
      k_mautp <- abs(k_mautp)
    }
    d_mautp <- min(alphaht^2 / ( (sigma2ols[1,1] / alphaht^2) + alphaht^2 ))
    if (d_mautp <= 0) {
      d_mautp <- abs(d_mautp)
    }
    Amautp <- (I - k_mautp^2 * (1 - d_mautp)^2 * solve(del + k_mautp * I) %*% solve(del + k_mautp * I)) %*% 
      solve(del + k_mautp * I) %*% (del)
    mautp <- Amautp %*% ols

# 7. Modified Ridge Type 
    d_hat_mrt <- min(sigma2ols[1,1] / alphaht^2)
    k_mrt <- p * sigma2ols[1,1] / sum((1 + d_hat_mrt) * alphaht^2)
    if (k_mrt <= 0) {
      k_mrt <- abs(k_mrt)
    }
    d_mrt1 <- (sigma2ols[1,1] / (k_mrt * alphaht^2)) - 1 
    d_mrt <- p / sum(1/d_mrt1)
    if (d_mrt <= 0 || d_mrt >= 1) {
      d_mrt <- d_hat_mrt
    }
    Amrt <- solve(del + k_mrt * (1 + d_mrt) * I) %*% del 
    mrt <- Amrt %*% ols

# 8. A New Biased Estimator 
    d_hat_dk <- min(sigma2ols[1,1] / alphaht^2)
    k_dk <- min(sigma2ols[1,1] / ( (1 + d_hat_dk) * ( (sigma2ols[1,1] / l) + 2 * alphaht^2 ) ))
    if (k_dk <= 0) {
      k_dk <- abs(k_dk)
    }
    d_dk <- min(sigma2ols[1,1] * l / (k_dk * (sigma2ols[1,1] + 2 * l * alphaht^2)))
    if (d_dk <= 0 || d_dk >= 1) {
      d_dk <- d_hat_dk
    }
    
    Adk <- solve(del + k_dk * (1 + d_dk) * I) %*% (del - k_dk * (1 + d_dk) * I)
    dk <- Adk %*% ols

# 9. Generalized Two Parameter 
    d_hat_gtp <- sigma2ols[1,1] / alphaht^2
    k_hat_gtp <- sigma2ols[1,1] / (alphaht^2 - d_hat_gtp * ( (sigma2ols[1,1] / l) + alphaht^2))
    k_gtp <- mean(k_hat_gtp)
    if (k_gtp <= 0) {
      k_gtp <- abs(k_gtp)
    }
    d_gtp <- (k_gtp * alphaht^2 - sigma2ols[1,1]) * l / (k_gtp * (sigma2ols[1,1] + alphaht^2 * l))
    d_gtp <- as.numeric(d_gtp)
    D_gtp <- diag(d_gtp)   
    
    Agtp <- solve(del + k_gtp * I) %*% (del + k_gtp * D_gtp)
    gtp <- Agtp %*% ols

# 10. Modified Liu Ridge Type
    d_hat_mlrt <-  max((alphaht^2 - sigma2ols[1,1]) * l / (sigma2ols[1,1] + l * alphaht^2))
    k_mlrt <- max((sigma2ols[1,1] * (l + d_hat_mlrt) - l * (1 - d_hat_mlrt) * alphaht^2) / ( (1 + d_hat_mlrt) * (l + 1) * alphaht^2))
    if (k_mlrt <= 0) {
      k_mlrt <- abs(k_mlrt)
    }
    d_mlrt <- sum((((l + k_mlrt * (l + 1)) * (l^2 - k_mlrt * l^2 - k_mlrt * l) *  alphaht^2 - sigma2ols[1,1] * l * (l^2 - k_mlrt * l^2 - k_mlrt * l))) /
                    (l + 1)^2) / 
      sum((sigma2ols[1,1] * (l^2 - k_mlrt * l^2 - k_mlrt * l) + (l - k_mlrt * (l + 1)) * (l^2 - k_mlrt * l^2 - k_mlrt * l) * alphaht^2) / (l + 1)^2)
    if (d_mlrt <= 0 || d_mlrt >= 1) {
      d_mlrt <- d_hat_mlrt
    }
    Amlrt <- solve(del + I) %*% (del + d_mlrt * I) %*% solve(del + k_mlrt * (1 + d_mlrt) * I) %*% del
    mlrt <- Amlrt %*% ols

# 11. New Biased Regression Two Parameter
    d_hat_nbr <- min(alphaht^2 / sigma2ols[1,1])
    k_nbr <- max((-(l^2 * alphaht^2 * (3 - d_hat_nbr) + sigma2ols[1,1] * l * (1 - d_hat_nbr)) / (2 * (sigma2ols[1,1] * d_hat_nbr + l * alphaht^2 * (1 + d_hat_nbr)))) +
                   (l * sqrt(l * alphaht^4 * (d_hat_nbr - 3)^2 + 2 * l * sigma2ols[1,1] * alphaht^2 * (5 - 2 * d_hat_nbr + d_hat_nbr^2) * sigma2ols[1,1]^2 * (1 + d_hat_nbr)^2) /
                      (2 * (sigma2ols[1,1] * d_hat_nbr + l * alphaht^2 * (1 + d_hat_nbr)))))
    if (k_nbr <= 0) {
      k_nbr <- abs(k_nbr)
    }
    
    d_nbr <- min((l^2 * (sigma2ols[1,1] - 3 * alphaht^2 * k_nbr) - l * k_nbr * (sigma2ols[1,1] + alphaht^2 * k_nbr)) /
                   (k_nbr * (k_nbr - l) * (sigma2ols[1,1] + l * alphaht^2)))
    if (d_nbr <= 0 || d_nbr >= 1) {
      d_nbr <- d_hat_nbr
    }
    R <- solve(del + k_nbr * I) %*% (del + k_nbr * d_nbr * I)
    W <- solve(del + k_nbr * I)
    M <- del - k_nbr * I
    nbr <- R %*% W %*% M %*% ols

# 12. On the Biased Two Parameter
    d_hat_btp <- min(alphaht^2 / sigma2ols[1,1])
    
    k_btp <- (sigma2ols[1,1] * l * (d_hat_btp - l) + alphaht^2 * l^2 * (d_hat_btp + 1)) / 
      (sigma2ols[1,1] * (1 + d_hat_btp) * (d_hat_btp - l) - alphaht^2 * l * (1 + d_hat_btp) * (2 * l - d_hat_btp + 1))
    k_btp <- mean(k_btp)
    if (k_btp <= 0) {
      k_btp <- abs(k_btp)
    }
    
    d_btp <- (-(sigma2ols[1,1] * l - sigma2ols[1,1] * k_btp + alphaht^2 * l^2 + sigma2ols[1,1] * l * k_btp + 2 * alphaht^2 * l^2 * k_btp) /
                (2 * sigma2ols[1,1] * k_btp + alphaht^2 * l * k_btp)) + 
      (sqrt(abs(alphaht^4 * l^2 * k_btp * (2 * l * k_btp + k_btp + l) + sigma2ols[1,1] * alphaht^2 * l^2 * k_btp * (k_btp - l) + sigma2ols[1,1] * alphaht^2 * (2 * sigma2ols[1,1] * l * k_btp + k_btp + l) + (sigma2ols[1,1]^2 * l * k_btp * (k_btp - l)))) /
         (sigma2ols[1,1] * k_btp + alphaht^2 * l * k_btp))
    d_btp <- mean(d_btp)
    if (d_btp <= 0 || d_btp >= 1) {
      d_btp <- d_hat_btp
    }
    
    E <- solve(del + I) %*% (del - d_btp * I)
    H <- solve(del + k_btp * (1 + d_btp) * I) %*% (del - k_btp * (1 + d_btp) * I)
    btp <- E %*% H %*% ols

# 13. New Two Parameter Estimator
    d_hat_ntp <- min(sigma2ols[1,1] / alphaht^2)
    k_ntp <- sigma2ols[1,1] / sum(d_hat_ntp * (2 * alphaht^2 + (sigma2ols[1,1] / l)))
    if (k_ntp <= 0) {
      k_ntp <- abs(k_ntp)
    }
    d_ntp <- sigma2ols[1,1] / sum(k_ntp * (2 * alphaht^2 + (sigma2ols[1,1] / l)))
    if (d_ntp <= 0) {
      d_ntp <- abs(d_ntp)
    }
    Antp <- solve(del + k_ntp * d_ntp * I) %*% (del - k_ntp * d_ntp * I)
    ntp <- Antp %*% ols

# 14. A new ridge type
    d_hat_nrt <- min(sigma2ols[1,1] / alphaht^2)
    k_nrt <- min((sigma2ols[1,1] / alphaht^2) - d_hat_nrt)
    if (k_nrt <= 0) {
      k_nrt <- abs(k_nrt)
    }
    d_nrt <- min((sigma2ols[1,1] / alphaht^2) - k_nrt)
    if (d_nrt <= 0) {
      d_nrt <- abs(d_nrt)
    }
    Anrt <- solve(del + (k_nrt + d_nrt) * I) %*% del
    nrt <- Anrt %*% ols

# 15. Liu-Kibria-Lukman
    k_lkl <- min(sigma2ols[1,1] / (2 * alphaht^2 + (sigma2ols[1,1] / l)))
    if (k_lkl <= 0) {
      k_lkl <- abs(k_lkl)
    }
    d_lkl <- min((l * (alphaht^2 - sigma2ols[1,1]) + l * k_lkl * (2 * alphaht^2 * l + alphaht^2 - sigma2ols[1,1])) / 
                   (sigma2ols[1,1] * (l - k_lkl) + alphaht^2 * l * (l - k_lkl)))
    if (d_lkl <= 0 || d_lkl >= 1) {
      d_lkl <- min(alphaht^2 / sigma2ols[1,1])
    }
    C <- solve(del + I) %*% (del + d_lkl * I)
    A <- solve(del + k_lkl * I) %*% (del -k_lkl * I)
    lkl <- C %*% A %*% ols

# 16. Two parameter ridge
    k_hat_tpr1 <- 0.5
    q_tpr1 <- sum(alphaht^2 * l / (l + k_hat_tpr1)) / sum((sigma2ols[1,1] * l + alphaht^2 * l^2) / (l + k_hat_tpr1)^2)
    k_tpr1 <- mean((q_tpr1 * sum(sigma2ols[1,1] * l) + (q_tpr1 - 1) * sum(alphaht^2 * l^2)) / (sum(alphaht^2 * l^2)))

    q1 <- t(t(Z) %*% y) %*% solve(del + k_tpr1 * I) %*% t(Z) %*% y / 
      t(t(Z) %*% y) %*% solve(del + k_tpr1 * I) %*% del %*% solve(del + k_tpr1 * I) %*% t(Z) %*% y

    Atpr1 <- q1[1,1] * solve(del + k_tpr1 * I) %*% del 
    tpr1 <- Atpr1 %*% ols
# 17. New two parameter
    q_new1 <- min(eigen((del + k_hat_tpr1 * I) %*% solve(del))$values)
    Anew1 <- q_new1 * solve(del + k_hat_tpr1 * I) %*% del
    new1 <- Anew1 %*% ols
