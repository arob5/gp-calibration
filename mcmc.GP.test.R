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
                         u.prior.mean, u.prior.sd) {
  accept.count <- 0
  tau.samples <- matrix(NA, nrow = n.itr, ncol = 1)
  u.samples <- matrix(NA, nrow = n.itr, ncol = ncol(gp.obj$X.obs))
  
  for(itr in seq(1, n.itr)) {
    # TODO: adapt proposal variance
    
    # Propose new calibration parameters
    u.proposal <- TruncatedNormal::rtmvnorm(1, mu = c(u.curr), sigma = cov.proposal, lb = u.rng[,1], ub = u.rng[,2])
    
    # Sample sufficient statistics corresponding to current and proposed values of the calibration parameter
    if(SS.joint.sample) {
      SS.samples <- sample.SS(gp.obj, rbind(u.curr, u.proposal))
      SS.curr <- SS.samples[1]
      SS.proposal <- SS.samples[2]
    } else {
      SS.curr <- sample.SS(gp.obj, u.curr)
      SS.proposal <- sample.SS(gp.obj, u.proposal)
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
    if(accept.u.proposal(u.curr.cond.dens, u.proposal.cond.dens, u.curr, u.proposal, cov.proposal)) {
      u.curr <- u.proposal
      SS.curr <- SS.proposal
      accept.count <- accept.count + 1
    }
    u.samples[itr,] <- u.curr
   
    # Gibbs step: sample tau conditional on current value of u
    tau.curr <- sample.tau(n, SS.curr, tau.gamma.shape, tau.gamma.rate)
    tau.samples[itr, 1] <- tau.curr
     
  }
  
}


sample.SS <- function(gp.obj, X.pred) {

  gp.pred <- predict_gp(X.pred, gp.obj, pred.cov = TRUE)

  repeat {
    SS <- rmvnorm(1, mean = gp.pred$mean, sigma = gp.pred$cov)
    if(all(SS >= 0)) break
  }
  
  return(as.vector(SS))
  
}


sample.tau <- function(n, SS, gamma.shape, gamma.rate) {
  return(rgamma(1, shape = gamma.shape + 0.5*n, rate = gamma.rate + 0.5*SS))
}


calc.u.log.conditional.density <- function(u, n, tau, SS, u.prior.mean, u.prior.sd, rng) {
  0.5*log(tau) - 0.5*tau*SS - log(dnorm(u, u.prior.mean, u.prior.sd))
}


accept.u.proposal <- function(u.curr.cond.dens, u.proposal.cond.dens, u.curr, u.proposal, cov.proposal) {
  # Proposal adjustment for detailed balance
  log.dens.proposal <- u.proposal.cond.dens + TruncatedNormal::dtmvnorm(u.curr, u.proposal, cov.proposal, lb = rng[,1], 
                                                                        ub = rng[,2], log = TRUE, B = 100)
  log.dens.curr <- u.curr.cond.dens + TruncatedNormal::dtmvnorm(u.proposal, u.curr, cov.proposal, lb = rng[,1], 
                                                                ub = rng[,2], log = TRUE, B = 100)
  
  # Determine acceptance
  return(exp(log.dens.proposal - log.dens.curr) > runif(1))
  
}










