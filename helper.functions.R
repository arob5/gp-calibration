

# Log unnormalized isotropic multivariate normal density
dmvnorm.log.unnorm <- function(y, u, tau, f) {
  n <- length(y)
  mu.vec <- rep(f(u), n)
  (n/2)*log(tau) - (tau/2)*sum((y - mu.vec)^2)
}

# Log unnormalized isotropic multivariate normal density, as function of 
# sufficient statistic
dmvnorm.log.unnorm.SS <- function(SS, tau, n) {
  (n/2)*log(tau) - (tau/2)*SS
}

save.gaussian.llik.plot <- function(y.obs, X.pred, out.dir, file.name, f) {
  
  N.pred <- nrow(X.pred)
  llik.pred <- matrix(NA, nrow = N.pred, ncol = 1)
  for(i in seq(1, N.pred)) {
    llik.pred[i, 1] <- dmvnorm.log.unnorm(y.obs, X.pred[i,1], tau, f)
  }
  
  png(file.path(out.dir, file.name), width=600, height=350)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  matplot(X.pred, llik.pred, type = 'l', 
          lty=1, xlab="u", ylab="Unnormalized Log-Likelihood", 
          main = "Exact Log-Likelihood", col="red")
  dev.off()
  
}

# GP predictive mean and interval estimates at test points
save.gp.pred.mean.plot <- function(interval.pct, tau, y.obs, X, X.pred, SS, SS.pred, 
                                   gp.means, gp.se, llik.gp.uq, out.dir, base.file.name, f) {
  p <- 1 - (1 - interval.pct)/2
  N <- nrow(X)
  N.pred <- nrow(X.pred)
  
  gp.pred.upper <- qnorm(p, gp.means, gp.se)
  gp.pred.lower <- qnorm(p, gp.means, gp.se, lower.tail = FALSE)
  
  #
  # Plot 1: Sufficient statistic plot (GP means)
  #
  png(file.path(out.dir, paste0(base.file.name, "_SS.png")), width=600, height=350)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(X.pred, gp.means, xlab = 'u', type = 'l', lty = 1, 
       ylab = paste0('GP Predictive Mean and ', 100*interval.pct, '% CI'),
       main = 'GP Predictive Mean of SS', 
       ylim = c(min(gp.pred.lower), max(gp.pred.upper)),
       col = 'blue')
  points(X, SS, pch = 16, col="black")
  lines(X.pred, gp.pred.upper, lty=1, col="gray")
  lines(X.pred, gp.pred.lower, lty=1, col="gray")
  lines(X.pred, SS.pred, lty=1, col="red")
  legend("right", inset=c(-0.2,-0.3), 
         legend = c("GP pred mean", "Design points", paste0(100*interval.pct, "% CI"), "True SS"), 
         col = c("blue", "black", "gray", "red"), lty = c(1, NA, 1, 1), pch = c(NA, 16, NA, NA))
  dev.off()
  
  #
  # Plot 2: Unnormalized log-likelihood evaluated at GP means
  #
  
  # True log-likelihood at design points
  llik.design <- matrix(NA, nrow = N, ncol = 1)
  for(i in seq(1, N)) {
    llik.design[i,1] <- dmvnorm.log.unnorm(y.obs, X[i,1], tau, f)
  }
  
  # True log-likelihood at test points
  llik.pred <- matrix(NA, nrow = N.pred, ncol = 1)
  for(i in seq(1, N.pred)) {
    llik.pred[i, 1] <- dmvnorm.log.unnorm(y.obs, X.pred[i,1], tau, f)
  }
  
  llik.gp.mean <- dmvnorm.log.unnorm.SS(gp.means, tau, n)
  llik.gp.upper <- dmvnorm.log.unnorm.SS(gp.pred.upper, tau, n)
  llik.gp.lower <- dmvnorm.log.unnorm.SS(gp.pred.lower, tau, n)
  
  png(file.path(out.dir, paste0(base.file.name, "_llik.png")), width=600, height=350)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(X.pred, llik.gp.mean, type="l", lty=1, xlab="u", 
       ylab=paste0("Unnormalized Log-Likelihood and ", 100*interval.pct, "% CI"),
       main = "GP Mean Approx Log-Likelihood", col="blue")
  points(X, llik.design, pch=16, col="black")
  lines(X.pred, llik.gp.upper, lty=1, col="gray")
  lines(X.pred, llik.gp.lower, lty=1, col="gray")
  lines(X.pred, llik.pred, lty=1, col="red")
  legend("right", inset=c(-0.2,-0.3), 
         legend = c("GP pred mean", "Design points", paste0(100*interval.pct, "% CI"), "True llik"), 
         col = c("blue", "black", "gray", "red"), lty = c(1, NA, 1, 1), pch = c(NA, 16, NA, NA))
  dev.off()
  
  # Plot 3: "Integrated out" GP approx log-likelihood
  png(file.path(out.dir, "gp_integrated_out_llik.png"), width=600, height=350)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(X.pred, llik.gp.uq, type = "l", lty = 1, xlab = "u", ylab="Unnormalized Log-Likelihood",
       main = "GP Integrated Out Approx Log-Likelihood", col = "blue")
  points(X, llik.design, pch=16, col="black")
  lines(X.pred, llik.pred, lty=1, col="red")
  legend("right", inset=c(-0.2,-0.3), 
         legend = c("GP integrated out llik", "Design points", "True llik"), 
         col = c("blue", "black", "red"), lty = c(1, NA, 1), pch = c(NA, 16, NA))
  dev.off()
  
  
}



