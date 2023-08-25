# Test 

# TODO: try variance reduction; control variates or antithetic sampling

normal_mean <- 0
normal_sd <- seq(.1, 4, length.out = 50)

log_LN_means <- vector(mode = "numeric", length = length(normal_sd))
log_LN_means_MC <- vector(mode = "numeric", length = length(normal_sd))

K <- 100000

for(l in seq_along(normal_sd)) {

  # Closed-form. 
  log_LN_means[l] <- exp(normal_mean + normal_sd[l]^2/2)
  # log_LN_means[l] <- normal_mean + normal_sd[l]^2/2
  
  # Monte Carlo. 
  x_samp <- rnorm(n = K, mean = normal_mean, sd = normal_sd[l])
  log_LN_means_MC[l] <- mean(exp(x_samp))
  # log_LN_means_MC[l] <- matrixStats::logSumExp(x_samp) - log(K)
  
}


errs <- log_LN_means - log_LN_means_MC

plot(log_LN_means_MC, log_LN_means)
lines(c(min(log_LN_means_MC), max(log_LN_means_MC)), c(min(log_LN_means_MC), max(log_LN_means_MC)), col = "red")

plot(normal_sd, errs)
