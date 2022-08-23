
preprocess.settings <- function(settings) {
  
  # Expand relative directories into full paths
  settings$output.dir <- file.path(settings$base.dir, settings$output.dir)
  settings$mcmc.brute.force.stan.path <- file.path(settings$base.dir, settings$mcmc.brute.force.stan.path)
  settings$mcmc.gp.stan.path <- file.path(settings$base.dir, settings$mcmc.gp.stan.path)
  settings$mcmc.gp.mean.path <- file.path(settings$base.dir, settings$mcmc.gp.mean.path)
  
  if(is.null(settings$tau.true)) {
    settings$tau.true <- 1 / settings$sigma.true^2
  }
  
  if(is.null(settings$sigma.true)) {
    settings$sigma.true <- 1 / sqrt(settings$tau.true)
  }
  
  # Init tau.gamma.shape so that the hyperprior mean of tau equals the true tau value
  if(is.null(settings$tau.gamma.shape) & (settings$tau.gamma.rate == 1.0)) {
    settings$tau.gamma.shape <- settings$tau.true
  }
  
  # Init u.gaussian.mean so that hyperprior mean of u equals the true u value
  if(is.null(settings$u.gaussian.mean)) {
    settings$u.gaussian.mean <- settings$u.true
  }
  
  if(any(as.logical(settings[c("mcmc.gp.stan", "mcmc.gp.mean.stan", "mcmc.pecan")])) && is.null(settings$X)) {
    settings$X <- matrix(seq(qnorm(.01, settings$u.true, settings$sigma.true), 
                             qnorm(.99, settings$u.true, settings$sigma.true), length = settings$N), ncol = settings$k)
  }
  
  return(settings)
  
}


# Log unnormalized isotropic multivariate normal density
log.dmvnorm <- function(y, u, tau, f, normalize = TRUE) {
  n <- length(y)
  mu.vec <- rep(f(u), n)
  log.dens <- 0.5*n*log(tau) - 0.5*tau*sum((y - mu.vec)^2)
  
  if(normalize) return(log.dens - 0.5*n*log(2*pi))
  return(log.dens)
}

# Log unnormalized isotropic multivariate normal density, as function of 
# sufficient statistic
dmvnorm.log.unnorm.SS <- function(SS, tau, n) {
  (n/2)*log(tau) - (tau/2)*SS
}

save.gaussian.llik.plot <- function(y.obs, X.pred, out.dir, file.name, f, tau, normalize = TRUE) {
  
  N.pred <- nrow(X.pred)
  llik.pred <- matrix(NA, nrow = N.pred, ncol = 1)
  for(i in seq(1, N.pred)) {
    llik.pred[i, 1] <- log.dmvnorm(y.obs, X.pred[i,1], tau, f, normalize)
  }
  
  if(normalize) {
    y.lab <- "Log-Likelihood"
  } else {
    y.lab <- "Unnormalized Log-Likelihood"
  }
  
  png(file.path(out.dir, file.name), width=600, height=350)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  matplot(X.pred, llik.pred, type = 'l', 
          lty=1, xlab="u", ylab=y.lab, 
          main = "Exact Log-Likelihood", col="red")
  dev.off()
  
}


save.gp.pred.mean.plot <- function(gp.obj, X, X.pred, SS, SS.pred, f, interval.pct, out.path, exp.pred = FALSE) {
  # Determine size of confidence interval to plot
  p <- 1 - (1 - interval.pct)/2
  
  # Number of design and test points
  N <- nrow(X)
  N.pred <- nrow(X.pred)
  
  # Calculate predictive means, standard errors, and confidence intervals
  gp.pred <- predict_gp(X.pred, gp.obj)
  gp.pred.mean <- gp.pred$mean
  gp.pred.se <- sqrt(gp.pred$var)
  gp.pred.upper <- qnorm(p, gp.pred.mean, gp.pred.se)
  gp.pred.lower <- qnorm(p, gp.pred.mean, gp.pred.se, lower.tail = FALSE)
  
  # If specified, exponentiate the predictions (used if the response variable was the log of the sufficient 
  # statistic but want to save plot on the original scale). The sufficient statistic now follows
  # a log-normal process. 
  if(exp.pred) {
    gp.pred.mean.lnorm <- exp(gp.pred.mean + 0.5 * gp.pred.se^2)
    gp.pred.se.lnorm <- sqrt((exp(gp.pred.se^2) - 1) * exp(2*gp.pred.mean + gp.pred.se^2))
    gp.pred.upper <- qlnorm(p, gp.pred.mean, gp.pred.se)
    gp.pred.lower <- qlnorm(p, gp.pred.mean, gp.pred.se, lower.tail = FALSE)
    
    SS <- exp(SS)
    SS.pred <- exp(SS.pred)
    gp.pred.mean <- gp.pred.mean.lnorm
    gp.pred.se <- gp.pred.se.lnorm
  }
  
  # Sufficient statistic plot
  png(out.path, width=600, height=350)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(X.pred, gp.pred.mean, xlab = 'u', type = 'l', lty = 1, 
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
  
}


# GP predictive mean and interval estimates at test points
save.gp.pred.mean.plot.old <- function(interval.pct, tau, y.obs, X, X.pred, SS, SS.pred, 
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
  lines(X.pred, llik.gp.mean, lty=1, col="green")
  points(X, llik.design, pch=16, col="black")
  lines(X.pred, llik.pred, lty=1, col="red")
  legend("right", inset=c(-0.2,-0.3), 
         legend = c("GP integrated out llik", "GP pred mean llik","Design points", "True llik"), 
         col = c("blue", "green", "black", "red"), lty = c(1, 1, NA, 1), pch = c(NA, NA, 16, NA))
  dev.off()
  
  
}


save.posterior.intervals.plot <- function(posterior.samples, pars, out.path, int.prob = 0.5, 
                                          int.outer.prob = 0.9, point.est = "median", append.title = NULL) {
  color_scheme_set("red")
  if(is.null(append.title)) {
    append.title <- ""
  } else {
    append.title <- paste0(": ", append.title)
  }
  
  mcmc_intervals(posterior.samples, 
                 pars = pars, 
                 prob = int.prob, 
                 prob_outer = int.outer.prob, 
                 point_est = point.est) + 
    ggtitle(paste0("Posterior Intervals", append.title)) + 
    ylab("Parameter") + 
    xlab(paste0("Inner ", 100*int.prob, "%; Outer ", 100*int.outer.prob, "%"))
  ggsave(out.path, bg = "white")
}


save.posterior.hist.plot <- function(posterior.samples, pars, out.path, append.title = NULL) {
  color_scheme_set("red")
  if(is.null(append.title)) {
    append.title <- ""
  } else {
    append.title <- paste0(": ", append.title)
  }
  
  mcmc_hist(posterior.samples, pars = pars) + 
    ggtitle(paste0("Posterior Histogram", append.title))
  ggsave(out.path, bg = "white")
}


save.posterior.kernel.dens.plot <- function(posterior.samples, pars, out.path, append.title = NULL) {
  color_scheme_set("red")
  if(is.null(append.title)) {
    append.title <- ""
  } else {
    append.title <- paste0(": ", append.title)
  }
  
  mcmc_dens(posterior.samples, pars = pars) + 
    ggtitle(paste0("Posterior Kernel Density Estimates", append.title))
  ggsave(out.path, bg = "white")
}

save.trace.plot <- function(posterior.samples, pars, out.path, append.title = NULL) {
  color_scheme_set("mix-blue-red")
  if(is.null(append.title)) {
    append.title <- ""
  } else {
    append.title <- paste0(": ", append.title)
  }
  
  mcmc_trace(posterior.samples, pars = pars,
             facet_args = list(ncol = 1, strip.position = "left")) + 
    ggtitle(paste0("Trace Plots", append.title))
  ggsave(out.path, bg = "white")
  
}


save.posterior.plots <- function(posterior.samples, pars, out.dir, int.prob = 0.5, int.outer.prob = 0.9,
                                 point.est = "median", append.filename = NULL, append.title = NULL) {
  color_scheme_set("red")                

  # Intervals plot
  plot.path <- file.path(out.dir, paste0("post.int", append.filename, ".png"))
  save.posterior.intervals.plot(posterior.samples, pars, plot.path, int.prob, 
                                int.outer.prob, point.est, append.title)
  
  # Histogram
  plot.path <- file.path(out.dir, paste0("post.hist", append.filename, ".png"))
  save.posterior.hist.plot(posterior.samples, pars, plot.path, append.title)
  
  # Kernel Density
  plot.path <- file.path(out.dir, paste0("post.kern.dens", append.filename, ".png"))
  save.posterior.kernel.dens.plot(posterior.samples, pars, plot.path, append.title)
  
  # Trace Plot
  plot.path <- file.path(out.dir, paste0("trace.plot", append.filename, ".png"))
  save.trace.plot(posterior.samples, pars, plot.path, append.title)
  
}


par.mcmc.results.to.arr <- function(mcmc.results, pars, n.mcmc, warmup.frac) {
  n.warmup <- ceiling(warmup.frac * n.mcmc)
  warmup.sel <- seq(1, n.warmup)
  n.samples <- n.mcmc - n.warmup
  samples <- array(NA, c(n.samples, length(mcmc.results), length(pars)))
  
  for(chain in seq_along(mcmc.results)) {
    for(p in seq_along(pars)) {
      samples[,chain,p] <- mcmc.results[[chain]][[pars[p]]][-warmup.sel]
    }
  }
  
  dimnames(samples)[[2]] <- paste0("chain:", seq_along(mcmc.results))
  dimnames(samples)[[3]] <- pars
  
  return(samples)
  
}


mcmc.summary <- function(samples, pars) {
  
  summary.list <- list(summary = create.mcmc.summary.table(samples, pars), 
                       c_summary = lapply(seq_len(dim(samples)[2]), 
                                          function(chain) create.mcmc.summary.table(samples, pars, chain)))
  
  return(summary.list)  

}

create.mcmc.summary.table <- function(samples, pars, chain = NULL) {
  
  # If chain is NULL, calculate statistics using samples from all chains
  if(is.null(chain)) {
    chain <- seq_len(dim(samples)[2])
  }
  
  summary.table <- matrix(NA, nrow = length(pars), ncol = 4, 
                          dimnames = list(pars, c("mean", "sd", "n_eff", "Rhat")))
  summary.table[,"mean"] <- apply(samples[,chain, pars, drop=FALSE], 3, mean)
  summary.table[,"sd"] <- apply(samples[,chain, pars, drop=FALSE], 3, sd)
  summary.table[,"n_eff"] <- apply(samples[,chain, pars, drop=FALSE], 3, ess_bulk)
  summary.table[,"Rhat"] <- apply(samples[,chain, pars, drop=FALSE], 3, Rhat)
  
  return(summary.table)
  
}












