

# Log unnormalized isotropic multivariate normal density
dmvnorm.log.unnorm <- function(y, u, tau, f) {
  n <- length(y)
  mu.vec <- rep(f(u), n)
  (n/2)*log(tau) - (tau/2)*sum((y - mu.vec)^2)
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


# True log-likelihood at test points


