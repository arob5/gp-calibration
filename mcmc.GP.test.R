# mcmc.GP.test.R
# Contains function mcmc.GP.test(), which is a simplified version of the PEcAn
# function mcmc.GP() to be used for testing purposes and comparison to 
# alternative algorithms. 
#
# Andrew Roberts
# Dependencies: helper.functions.R

library(TruncatedNormal)
library(mvtnorm)

# TODO: 
#   - generalize so that each dimension can be updated independently 
#   - Pass in priors for u instead of hard-coding Gaussian prior
#   - Generalize u prior function to allow for multivariate u

mcmc.GP.test <- function(gp.obj, n, n.itr, u.rng, SS.joint.sample, resample.tau, tau.gamma.shape, tau.gamma.rate, 
                         u.prior.mean, u.prior.sd, u0, proposal.vars, adapt.frequency, adapt.min.scale, accept.rate.target, 
                         log.normal.process) {

  accept.count <- 0
  tau.samples <- matrix(NA, nrow = n.itr, ncol = 1)
  u.samples <- matrix(NA, nrow = n.itr, ncol = ncol(gp.obj$X))
  SS.samples <- matrix(NA, nrow = n.itr, ncol = 1)
  
  # Initial parameter values
  SS0 <- sample.SS(gp.obj, u0, log.normal.process)
  tau.curr <- sample.tau(n, SS0, tau.gamma.shape, tau.gamma.rate)
  u.curr <- u0
  
  # Initial proposal covariance matrix
  cov.proposal <- diag(proposal.vars, nrow = length(proposal.vars))
  
  for(itr in seq(1, n.itr)) {
    
    # Adapt proposal covariance matrix
    if((itr > 2) && ((itr - 1) %% adapt.frequency) == 0) {
      cov.proposal <- adapt.cov.proposal(cov.proposal, u.samples[(itr - adapt.frequency):(itr - 1),,drop=FALSE], 
                                         adapt.min.scale, accept.count / adapt.frequency, accept.rate.target)
      accept.count <- 0
    }
    
    # Propose new calibration parameters
    u.proposal <- TruncatedNormal::rtmvnorm(1, mu = c(u.curr), sigma = cov.proposal, lb = u.rng[1,], ub = u.rng[2,])
    u.proposal <- matrix(u.proposal, nrow = 1)
    
    # Sample sufficient statistics corresponding to current and proposed values of the calibration parameter
    if(SS.joint.sample) {
      SS.joint.samples <- sample.SS(gp.obj, rbind(u.curr, u.proposal), log.normal.process)
      SS.curr <- SS.joint.samples[1]
      SS.proposal <- SS.joint.samples[2]
    } else {
      SS.curr <- sample.SS(gp.obj, as.matrix(u.curr, nrow = 1), log.normal.process)
      SS.proposal <- sample.SS(gp.obj, u.proposal, log.normal.process)
    }
    
    # Re-sample precision parameter tau
    if(resample.tau) {
      tau.proposal <- sample.tau(n, SS.proposal, tau.gamma.shape, tau.gamma.rate)
    } else {
      tau.proposal <- tau.curr
    }
    
    # Calculate log-conditional density u given tau
    u.curr.cond.dens <- calc.u.log.conditional.density(u.curr, n, tau.curr, SS.curr, u.prior.mean, u.prior.sd)
    u.proposal.cond.dens <- calc.u.log.conditional.density(u.proposal, n, tau.proposal, SS.proposal, u.prior.mean, u.prior.sd)
                                                           
    # Accept or reject proposed u value
    if(accept.u.proposal(u.curr.cond.dens, u.proposal.cond.dens, u.curr, u.proposal, cov.proposal, u.rng)) {
      u.curr <- u.proposal
      SS.curr <- SS.proposal
      accept.count <- accept.count + 1
    }
    SS.samples[itr,] <- SS.curr
    u.samples[itr,] <- u.curr
   
    # Gibbs step: sample tau conditional on current value of u
    tau.curr <- sample.tau(n, SS.curr, tau.gamma.shape, tau.gamma.rate)
    tau.samples[itr, 1] <- tau.curr
  }
  
  samples.list <- list(u = u.samples, tau = tau.samples, SS = SS.samples)
  
  return(samples.list)
  
}


sample.SS <- function(gp.obj, X.pred, log.normal.process) {

  gp.pred <- predict_gp(X.pred, gp.obj, pred.cov = TRUE)

  # Sampled value is log(SS), so exponentiate.
  if(log.normal.process) {
    return(exp(as.vector(rmvnorm(1, mean = c(gp.pred$mean), sigma = gp.pred$cov))))
  }
  
  repeat {
    SS <- rmvnorm(1, mean = c(gp.pred$mean), sigma = gp.pred$cov)
    if(all(SS >= 0)) break
  }
  
  return(as.vector(SS))
  
}


sample.tau <- function(n, SS, gamma.shape, gamma.rate) {
  return(rgamma(1, shape = gamma.shape + 0.5*n, rate = gamma.rate + 0.5*SS))
}


calc.u.log.conditional.density <- function(u, n, tau, SS, u.prior.mean, u.prior.sd, rng) {
  
  log.prior <- sum(sapply(seq_along(u.prior.mean), function(j) dnorm(u[,j], u.prior.mean[j], u.prior.sd[j], log = TRUE)))
  0.5*n*log(tau) - 0.5*tau*SS + log.prior
  
}


accept.u.proposal <- function(u.curr.cond.dens, u.proposal.cond.dens, u.curr, u.proposal, cov.proposal, rng) {
  # Proposal adjustment for detailed balance
  log.dens.proposal <- u.proposal.cond.dens + TruncatedNormal::dtmvnorm(u.curr, u.proposal, cov.proposal, lb = rng[1], 
                                                                        ub = rng[2], log = TRUE, B = 100)
  log.dens.curr <- u.curr.cond.dens + TruncatedNormal::dtmvnorm(u.proposal, u.curr, cov.proposal, lb = rng[1], 
                                                                ub = rng[2], log = TRUE, B = 100)
  
  # Determine acceptance
  return(as.logical(exp(log.dens.proposal - log.dens.curr) > runif(1)))
  
}


adapt.cov.proposal <- function(cov.proposal, sample.history, min.scale, accept.rate, accept.rate.target) {
  if(accept.rate == 0) {
    return(min.scale * cov.proposal)
  } else {
    cor.estimate <- stats::cor(sample.history)
    scale.factor <- max(accept.rate / accept.rate.target, min.scale)
    stdev <- apply(sample.history, 2, stats::sd)
    scale.mat <- scale.factor * diag(stdev, nrow = length(stdev))
    return(scale.mat %*% cor.estimate %*% scale.mat)
  }
}




