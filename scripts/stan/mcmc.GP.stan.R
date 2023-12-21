##' Samples from an approximate posterior distribution, in which the true likelihood
##' is replaced by a Gaussian Process (GP) approximation. The GP is already assumed 
##' to be fit prior to calling this function and must be passed as an argument. To 
##' account for the uncertainty introduced by the GP approximation, the likelihood
##' evaluation at a specific parameter value is averaged over the predictive distribution
##' given by the GP; i.e. the GP is analytically integrated out. This function automatically
##' generated a .stan file, which is then used to run MCMC. 
##'
##' @name mcmc.GP.stan
##' @title mcmc.GP.stan
##' @export
##'
##' @param gp.param.list Named list of tuned Gaussian Process hyperparameters. See 
##'                      'create.gp.params.list()'
##' @param x0 initial values
##' @param nmcmc number of iterations
##' @param rng range of knots
##' @param ar.target acceptance rate target
##' @param settings PEcAn settings list
##' @param priors prior list
##' @param run.block is this a new run or making the previous chain longer
##' @param n.of.obs number of observations
##' @param llik.fn list that contains likelihood functions
##' @param hyper.pars hyper parameters
##' 
##' @author Andrew Roberts
mcmc.GP.stan <- function(gp.param.list) {
  
  
  
}


