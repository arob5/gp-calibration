
generate.stan.script <- function(prior.table, write.to.file = TRUE) {
  
  data.block <- define.stan.data.block()
  params.block <- define.stan.params.block()
  model.block <- define.stan.model.block(prior.table)
  stan.script <- paste(data.block, params.block, model.block, sep = '\n\n')
  return(stan.script)
}


define.stan.data.block <- function() {
  'data {
    int<lower=1> n;  // Number of observations for output variable
    int<lower=1> k;  // Number of calibration parameters
    int<lower=1> N;  // Number of knots
    matrix[N, k] X;  // Design matrix
    vector[N] y;    // Process model evaluations at design points
    real<lower=0> a; // Shape parameter for Gamma hyperprior on tau
    real<lower=0> b; // Rate parameter for Gamma hyperprior on tau
    vector<lower=0>[k] gp_rho; // Vector of lengthscale parameters
    real<lower=0> gp_alpha; // Marginal standard deviation
    real<lower=0> gp_sigma; // Nugget (standard deviation)
    vector[k] gp_beta; // Vector of regression coefficients for linear mean function
  }'
}


define.stan.params.block <- function() {
  'parameters {
    real<lower=0> tau; // Precision parameter in Gaussian model 
    vector[k] u;       // Calibration parameters
   }'
}


define.stan.model.block <- function(prior.table) {
  
  model.block <- 'model {'
  model.block <- append.priors(model.block, prior.table)
  model.block <- paste(model.block, '\n}')
}


append.priors <- function(model.block, prior.table) {
  # Parameters in statistical noise model
  tau.prior <- 'tau ~ gamma(a, b);'
  
  # Calibration paramters
  calibration.param.priors <- apply(prior.table, 1, convert.prior.to.stan)
  calibration.param.priors <- paste0('u[', 1:nrow(prior.table), '] ~ ', calibration.param.priors, ';')
  calibration.param.priors <- paste0(calibration.param.priors, collapse = '\n')
  
  # Combine and append to existing model block
  model.block <- paste(model.block, tau.prior, calibration.param.priors, sep = '\n')
  
  return(model.block)
}


convert.prior.to.stan <- function(prior.info) {
  # Supports weibull, normal, uniform, gamma
  dist <- prior.info[['distn']]
  valid.dists <- c('unif', 'gamma', 'weibull', 'normal')
  if(!(dist %in% valid.dists)) {
    stop('Distribution ', dist, ' is not valid.')
  }
  
  # Convert between PEcAn and stan specification, if necessary. 
  if(dist == 'unif') dist <- 'uniform'
  
  stan.string <- paste0(dist, '(', prior.info[['parama']])
  if(!is.na(prior.info[['paramb']])) stan.string <- paste0(stan.string, ',', prior.info[['paramb']])
  stan.string <- paste0(stan.string, ')')

  return(stan.string)
  
}


##' Extracts tuned hyperparameters from a Gaussian Process (GP) regression   
##' and stores them in a list intended to be passed to generate.stan.script().
##' These parameters are then used to define the kernel/covariance function
##' in stan. Currently only Gaussian/squared exponential kernels are supported.
##' This function is intended to be independent of the specific package used to
##' fit the GP; helper functions can be added to accommodate new packages as 
##' needed. 
##' 
##' @name create.gp.params.list
##' @title create.gp.params.list
##' @export
##' 
##' @param gp Fit GP object from one of the supported libraries. Currently this is
##'           limited to an object of class gp from package mlegp.
##' @param library String specifying the library used to fit the GP. 
##' 
##' @author Andrew Roberts
create.gp.params.list <- function(gp, library) {

  if(library == 'mlegp') {
    gp.params <- convert.mlegp.params(gp)
  }
  
  return(gp.params)
  
}


convert.mlegp.params <- function(gp) {
  
  if(length(gp$nugget) > 1) {
    stop('convert.mlegp.params() does not support non-constant nugget.')
  }
  
  if(length(gp$Bhat) > 1) {
    stop('convert.mlegp.params() does not support non-constant mean function.')
  }
  
  gp.params <- list(gp_rho = gp$beta,
                    gp_alpha = sqrt(gp$sig2), 
                    gp_sigma = gp$nugget, 
                    gp_beta = gp$Bhat)
                    
  
  return(gp.params)
  
}






