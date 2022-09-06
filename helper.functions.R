
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
  if(is.null(settings$tau.gamma.shape)) {
    settings$tau.gamma.shape <- settings$tau.true * settings$tau.gamma.rate
  }
  
  for(j in seq_len(settings$k)) {
    
    # Init u.gaussian.mean so that hyperprior mean of u equals the true u value
    if((length(settings$u.gaussian.mean) < j) || is.null(settings$u.gaussian.mean[[j]])) {
      settings$u.gaussian.mean[[j]] <- settings$u.true[[j]]
    }
    
    # If Gaussian std devs are NULL, initialize so that the priors have the specified
    # coefficients of variation. 
    if((length(settings$u.gaussian.sd) < j) || is.null(settings$u.gaussian.sd[[j]])) {
      if(is.null(settings$u.prior.coef.var[[j]])) {
        stop("Either <u.gaussian.sd> or <u.prior.coef.var> must be non-NULL for calibration parameter ", j)
      }
    
      settings$u.gaussian.sd[[j]] <- settings$u.prior.coef.var[[j]] * settings$u.gaussian.mean[[j]]
    }
    
    # If calibration parameter range is NULL, set to (-Inf, Inf).
    if((length(settings$u.rng) < j) || is.null(settings$u.rng[[j]])) {
      settings$u.rng[[j]] <- c(-Inf, Inf)
    }
  }
  
  # If not explicitly passed, generate design via Latin Hypercube Sampling
  run.algs <- settings[c("mcmc.gp.stan", "mcmc.gp.mean.stan", "mcmc.pecan")]
  if(!any(sapply(run.algs, is.null)) && any(as.logical(run.algs)) && is.null(settings$X)) {
    LHS.X.pred <- is.null(settings$X.pred)
    if(LHS.X.pred) {
      X.LHS <- LHS.train.test(settings$N, settings$N.pred, settings$k, settings$u.gaussian.mean, settings$u.gaussian.sd)
      settings$X <- X.LHS$X
      settings$X.pred <- X.LHS$X.test
    } else {
      settings$X <- randomLHS(settings$N, settings$k)
    }
  }
  
  # Convert lists to more convenient types
  settings$u.true <- as.vector(settings$u.true, mode = "numeric")
  settings$u.gaussian.mean <- as.vector(settings$u.gaussian.mean, mode = "numeric")
  settings$u.gaussian.sd <- as.vector(settings$u.gaussian.sd, mode = "numeric")
  settings$u.prior.coef.var <- as.vector(settings$u.prior.coef.var, mode = "numeric")
  settings$u.rng <- matrix(c(sapply(settings$u.rng, function(x) x)), nrow = settings$k, ncol = 2, byrow = TRUE)
  
  return(settings)
  
}


# Log isotropic multivariate normal density
log.dmvnorm <- function(y, u, tau, f, normalize = TRUE) {
  n <- length(y)
  mu.vec <- rep(f(u), n)
  log.dens <- 0.5*n*log(tau) - 0.5*tau*sum((y - mu.vec)^2)
  
  if(normalize) return(log.dens - 0.5*n*log(2*pi))
  return(log.dens)
}

# Log unnormalized isotropic multivariate normal density, as function of 
# sufficient statistic
dmvnorm.log.SS <- function(SS, tau, n, normalize = TRUE) {
  log.dens <- (n/2)*log(tau) - (tau/2)*SS
  if(normalize) return(log.dens - 0.5*n*log(2*pi))
  return(log.dens)
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


save.gp.pred.llik.plot <- function(gp.log.obj, tau, n, X, X.pred, log.SS, log.SS.pred, out.path) {

  # Calculate predictive means and standard errors
  gp.pred <- predict_gp(X.pred, gp.log.obj)
  gp.pred.mean <- gp.pred$mean
  gp.pred.se <- sqrt(gp.pred$var)
  
  # Calculate likelihood approximations
  # M <- lnorm.mgf.estimate.analytic(0.5 * tau, gp.pred.mean, gp.pred.se)
  M <- rep(NA, length(gp.pred.mean))
  scale.factors <- rep(NA, length(gp.pred.mean))
  for(i in seq_along(M)) {
    mgf.results <- lnorm.mgf.estimate(0.5*tau, gp.pred.mean[i], gp.pred.se[i], scale = TRUE)
    # mgf.results <- lnorm.mgf.estimate.analytic(0.5 * tau, gp.pred.mean[i], gp.pred.se[i], scale = TRUE)
    # mgf.results <- lognormal_mgf_numerical_approx(0.5 * tau, gp.pred.mean[i], gp.pred.se[i], 
    #                                               1000, 1e-10, 1, scale = TRUE)
    M[i] <- mgf.results$mgf
    scale.factors[i] <- mgf.results$scale.factor
  }

  logM <- log(M)

  # Likelihood evaluations to plot
  llik.gp.approx <- 0.5 * n * log(tau) + logM - 0.5 * n * log(2*pi) - scale.factors
  llik.exact.design <- dmvnorm.log.SS(exp(log.SS), tau, n)
  llik.exact.pred <- dmvnorm.log.SS(exp(log.SS.pred), tau, n)
  
  # Log-likelihood plot
  range.values <- c(llik.gp.approx, llik.exact.design, llik.exact.pred)
  png(out.path, width=600, height=350)
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(X.pred, llik.gp.approx, xlab = 'u', type = 'l', lty = 1, 
      ylim = range(range.values[!is.infinite(range.values)]),
      ylab = "Log-Likelihood",
      main = 'GP Approx Log-Likelihood', 
      col = 'blue')
  points(X, llik.exact.design, pch = 16, col="black")
  lines(X.pred, llik.exact.pred, lty=1, col="red")
  legend("right", inset=c(-0.2,-0.3),
         legend = c("GP Approx", "Design points", "True llik"), col = c("blue", "black", "red"),
         lty = c(1, NA, 1), pch = c(NA, 16, NA))
  dev.off()
  
  return(invisible(M))
  
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
  
  llik.gp.mean <- dmvnorm.log.SS(gp.means, tau, n, FALSE)
  llik.gp.upper <- dmvnorm.log.SS(gp.pred.upper, tau, n, FALSE)
  llik.gp.lower <- dmvnorm.log.SS(gp.pred.lower, tau, n, FALSE)
  
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


save.SS.tau.samples.plot <- function(posterior.samples, out.path) {
  
  # Combine samples from all chains
  posterior.samples <- posterior.samples[,,c("tau", "SS")]
  SS <- as.vector(sapply(seq(1, dim(posterior.samples)[2]), function(chain) posterior.samples[, chain, "SS"]))
  tau <- as.vector(sapply(seq(1, dim(posterior.samples)[2]), function(chain) posterior.samples[, chain, "tau"]))
  
  # Add regression line
  reg <- lm(tau ~ SS, data = data.frame(SS = SS, tau = tau))

  # Save scatter plot
  png(out.path, width=600, height=350)
  plot(SS, tau, main = "Sufficient statistic GP samples vs. tau samples",
       xlab = "SS", ylab = "tau")
  abline(reg, col = "red")
  
  mtext(paste0("tau = ", format(round(reg$coefficients[1], 2), nsmall = 2), 
               " + ", format(round(reg$coefficients[2], 2), nsmall = 2), " * SS"), side = 4, col = "red")
  dev.off()

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

lnorm.mgf.estimate.analytic <- function(s, mu, sigma, scale = FALSE, scale.factor = NULL) {
  
  W <- lambertWp(s * exp(mu) * sigma^2)
  exp.arg <- 0.5 * (1 / sigma^2) * (W^2 + 2*W)
  
  if(scale) {
    if(is.null(scale.factor)) {
      scale.factor <- exp.arg
    }
  } else {
    scale.factor <- 0
  }
  
  
  list.out <- list(mgf = exp(scale.factor - exp.arg) / sqrt(1 + W), 
                   scale.factor = scale.factor)
  
  return(list.out)
  
}


# To check approximation, also estimate MGF via Monte Carlo and analytic approximation
lnorm.mgf.estimate <- function(s, mu, sigma, N = 10000, scale = FALSE, scale.factor = NULL) {
  lnorm.samples <- rlnorm(N, meanlog = mu, sdlog = sigma)
  if(scale) {
    if(is.null(scale.factor)) {
      scale.factor <- min(s * lnorm.samples)
    }
  } else {
    scale.factor <- 0
  }
  
  list.out <- list(mgf = mean(exp(scale.factor - s * lnorm.samples)), 
                   scale.factor = scale.factor)

  return(list.out)
  
}


lognormal_mgf_numerical_approx <- function(s, mu, sigma, num_eval, tol, M, scale = FALSE, scale.factor = NULL) {
  eps <- 0.5 * tol  
  cut_lower <-  mu - sigma * sqrt(2) * sqrt(log(M/eps))
  cut_upper <- log(1/s) + log(log(M/eps))
  dx <- (cut_upper - cut_lower) / num_eval
  x <- seq(cut_lower, cut_upper, length.out = num_eval) 
  fx_arg <- s * exp(x) + 0.5 * (x - mu)^2 / sigma^2
  
  if(scale) {
    if(is.null(scale.factor)) {
      scale.factor <- median(fx_arg)
    }
  } else {
    scale.factor <- 0
  }

  fx <- exp(scale.factor - fx_arg) 
  
  list.out <- list(mgf = 0.5 * dx * (1 / sqrt(2.0 * pi)) * (1 / sigma) * (fx[1] + fx[num_eval] + 2.0*sum(fx[2:(num_eval-1)])), 
                   scale.factor = scale.factor)
  
  return(list.out)
  
}


calc.SS <- function(X, f, y.obs, log.SS = FALSE) {
  # Run full model at given locations
  y <- apply(X, 1, f)
  
  # Calculate sufficient statistic at training and test locations
  N <- nrow(X)
  SS <- rep(0, N)

  for(i in seq(1, N)) {
    SS[i] <- sum((y[i] - y.obs)^2)
  }
  
  if(log.SS) return(log(SS))
  return(SS)

}


LHS.train.test <- function(N, N.test, k, mu.vals, sd.vals) {
  # Generate train and test sets
  X.combined <- randomLHS(N + N.test, k)
  X <- X.combined[1:N,,drop=FALSE]
  X.test <- X.combined[-(1:N),,drop=FALSE]
  
  # Apply inverse CDS transform using prior distributions.
  for(j in seq_len(k)) {
    X[,j] <- qnorm(X[,j], mu.vals[j], sd.vals[j])
    X.test[,j] <- qnorm(X.test[,j], mu.vals[j], sd.vals[j])
  }
  
  return(list(X = X, X.test = X.test))
  
}


fit.GPs <- function(gp.libs, X, SS, log.SS = FALSE) {
  
  # If modeling log(SS) instead of SS
  if(log.SS) {
    SS <- log(SS)
  }
  
  # Fit GP regressions
  gp.fits <- vector(mode = "list", length = length(gp.libs))
  gp.fits.index <- 1
  if("mlegp" %in% gp.libs) {
    gp.fits[[gp.fits.index]] <- mlegp(X, SS, nugget.known = 0, constantMean = 1)
    names(gp.fits)[[gp.fits.index]] <- "mlegp"
    gp.fits.index <- gp.fits.index + 1
  }
  
  if("hetGP" %in% settings$gp.library) {
    gp.fits[[gp.fits.index]] <- mleHomGP(X, SS, covtype = "Gaussian", known = list(g = .Machine$double.eps))
    names(gp.fits)[[gp.fits.index]] <- "hetGP"
    gp.fits.index <- gp.fits.index + 1
  }
  
  if(gp.fits.index != length(gp.fits) + 1) {
    stop("Invalid GP library detected.")
  }
  
  # Order to match order of gp.libs
  return(gp.fits[gp.libs])
  
}


save.cv.SS.plot <- function(cv.obj, out.path, log.SS = FALSE) {
  
  num.itr.cv <- cv.obj$num.itr.cv
  num.pred <- length(cv.obj$SS.test[[1]])
  
  SS.true <- c(sapply(cv.obj$SS.test, as.numeric))
  SS.pred <- c(sapply(cv.obj$hetGP$pred.mean, as.numeric))
  cv.num <- rep(seq(1, num.itr.cv), rep(num.pred, num.itr.cv))
  
  cv.colors <- colorspace::rainbow_hcl(num.itr.cv)
  dt <- data.table(SS_pred = SS.pred, SS_true = SS.true, cv_itr = paste0("cv", cv.num))
  
  png(out.path, width=600, height=350)
  if(log.SS) {
    plot(log(SS_true) ~ log(SS_pred), col = colorspace::rainbow_hcl(num.itr.cv)[cv.num], data = dt)
  } else {
    plot(SS_true ~ SS_pred, col = colorspace::rainbow_hcl(num.itr.cv)[cv.num], data = dt)
  }
  dev.off()

}




