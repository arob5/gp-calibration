#
# gp_mcmc_functions.r
# MCMC algorithms utilizing a Gaussian process (GP) emulator to accelerate 
# inference. 
#
# Andrew Roberts
#
# Depends:
#    gpWrapper, llikEmulator, general_helper_functions, mcmc_helper_functions

library(parallel)

# -----------------------------------------------------------------------------
# Code Standards and Conventions:
#
# This file contains MCMC algorithms designed to function with objects from the   
# `llikEmulator` class. The emphasis is on log-likelihood emulators constructed
# from Gaussian process (GP) approximations. The philosophy is for these 
# algorithms to be able to work with an exact (non-emulated) likelihood as well.
# That is, for llikEmulator objects with attribute `exact_llik = TRUE`, the 
# approximate MCMC algorithms here should sample from the exact posterior.
#
# Note that there are multiple types of parameters that may be samples in an 
# MCMC algorithm. The main parameters of interest here are those that 
# correspond to the inputs of the log-likelihood emulator, and are by 
# convention referred to as "par" in this file. The remaining set of parameters
# are typically those defining the likelihood; e.g., a variance parameter. 
# These are referred to as "likelihood parameters".
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Core MCMC algorithms. 
#
# At present, there are three primary approximate MCMC schemes, as well as one 
# wrapper function to allow access to the samplers in the `BayesianTools` R 
# package for exact MCMC.
#    1.) mcmc_noisy_llik: An MCMC algorithm defined with respect to a stochastic
#        approximation of the log-likelihood. This function implements methods 
#        that only require sampling from the log-likelihood surrogate, thus 
#        representing "noisy" versions of a traditional Metropolis-Hastings
#        scheme. The function arguments allow for a user to run a Monte Carlo 
#        within Metropolis or a pseudo marginal algorithm.
#    2.) mcmc_lik_approx: A Metropolis-Hastings algorithm defined with respect 
#        to a deterministic approximation of the true likelihood (alternatively,
#        can view this as an approximation to the unnormalized posterior 
#        density). This function relies of the `calc_lik_approx()` method of 
#        the llikEmulator object.
#    3.) mcmc_acc_prob_approx: An approximate Metropolis-Hastings like algorithm 
#        defined by a deterministic approximation to the true acceptance 
#        probability. This function relies on the helper function 
#        `get_gp_mh_acc_prob_approx()`.
#    4.) mcmc_bt_wrapper: A wrapper that provides access to the MCMC algorithms
#        implemented in the `BayesianTools` package. Hence, this function is 
#        only for inference with exact log-likelihood objects.
#
# At this point, these algorithms are only applicable when the true likelihood
# is of a multiplicative Gaussian form. This is due to the way that the 
# likelihood parameters are treated by the functions, and should be 
# generalized to more general likelihoods in the future.
# -----------------------------------------------------------------------------

# Need to check that the llik_em is correct for this function. 
# Infer which sig2 need to be learned by the `lik_par_fixed` attribute.
# How to map `sig2_prior+_params` onto the proper parameters? 
# TODO: need to have checking in llikSumEmulator to make sure the same input 
# parameters are used for each term. 

mcmc_noisy_llik <- function(llik_em, par_prior, par_init=NULL, sig2_init=NULL,
                            sig2_prior=NULL, mode="mcwmh", use_joint=TRUE,  
                            n_itr=50000L, cov_prop=NULL, log_scale_prop=NULL, 
                            adapt_cov_prop=TRUE, adapt_scale_prop=TRUE, 
                            adapt=adapt_cov_prop||adapt_scale_prop, 
                            return_prop_sd=FALSE, accept_rate_target=0.24, 
                            adapt_factor_exponent=0.8, adapt_factor_numerator=10, 
                            adapt_interval=200, ...) {
  # Args:
  #    llik_em: an object for which `is_llik_em(llik_em)` is TRUE.
  #    par_prior: data.frame defining the prior distribution for `par`. 
  #    par_init: numeric, d-dimensional vector, the initial condition. If NULL,
  #              samples from the prior to determine the initial condition.
  #    sig2_init: numeric, the initial condition for the variance parameters. 
  #    sig2_prior: list defining the prior distribution on the variance 
  #                parameters.
  #    mode: character, either "mcwmh" or "pseudo-marginal" to run either the Monte
  #          Carlo within Metropolis-Hastings or Pseudo-Marginal algorithm.
  #    use_joint: logical, whether to utilize the log-likelihood emulator
  #               joint distribution across different parameter values, or treat
  #               the emulator as independent on a parameter-by-parameter basis.
  #    n_itr: integer, the number of MCMC iterations.
  #    return_prop_sd: logical, if TRUE returns the trajectory of the standard 
  #                    deviations of the proposal distribution; that is, the 
  #                    square root of the diagonal of the proposal covariance.
  #    Remaining arguments define the proposal covariance adaptation, and are 
  #    passed to `adapt_MH_proposal_cov()`.
  #
  # Returns:
  # list, with the main element named "samp". This element is itself a list with
  # elements "par" and "sig2". The elements of the outer list other than "samp"  
  # store other information regarding the MCMC run (initial conditions, proposal 
  # adaptation settings, etc.). 
  # TODO: ensure that `sig2_init` is named vector, if non-NULL. 
  # And that the names include the names of the outputs where sig2 must be 
  # learned. Ensure `par_init` has names as well. 

  # TODO: for some reason, `inherits(llik_em, "llikEmulator")` always 
  # evaluates to FALSE when this is called from parLapply(). Need to figure 
  # out why. 
  # assert_that(is_llik_em(llik_em))

  # Objects to store samples and other outputs.
  d <- llik_em$dim_input
  par_samp <- matrix(nrow=n_itr, ncol=d)
  chain_info <- matrix(nrow=n_itr, ncol=d+2)
  colnames(par_samp) <- llik_em$input_names
  colnames(chain_info) <- c("llik", "lprior", llik_em$input_names)
  
  # Set initial conditions. 
  if(is.null(par_init)) par_init <- sample_prior(par_prior, n=1L)[1,]
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  lprior <- get_lprior_dens(par_prior)
  lprior_curr <- lprior(par_curr)
  
  # Setup for `sig2` (observation variances). Safe to assume that all of the 
  # non-fixed likelihood parameters are `sig2` since this is verified by 
  # `validate_args_mcmc_gp_noisy()` above.
  learn_sig2 <- !is.null(sig2_prior)
  term_labels_learn_sig2 <- names(learn_sig2)[learn_sig2]
  include_sig2_Gibbs_step <- (length(term_labels_learn_sig2) > 0)
  sig2_curr <- sig2_init[term_labels_learn_sig2] # Only includes non-fixed variance params.
  
  if(include_sig2_Gibbs_step) {
    N_obs <- unlist(llik_em$get_llik_attr("N_obs", labels=term_labels_learn_sig2))
    sig2_samp <- matrix(nrow=n_itr, ncol=length(sig2_curr_learn))
    sig2_samp[1,] <- sig2_curr
  } else {
    sig2_curr <- NULL
    sig2_samp <- NULL
  }
  
  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0
  times_adapted <- 0
  
  # Store initial chain info.
  llik_curr <- drop(llik_em$sample(par_curr, ...))
  chain_info[1,] <- c(llik_curr, lprior_curr, 
                      exp(log_scale_prop) * sqrt(diag(cov_prop)))
  
  # Variable to store error condition, if it occurs.
  err <- NULL
  
  # Main MCMC loop. 
  tryCatch(
    {
      for(itr in 2L:n_itr) {
        #
        # Metropolis step for calibration parameters.
        #
        
        # Random walk proposal. 
        par_prop <- par_curr + (exp(log_scale_prop) * L_cov_prop %*% matrix(rnorm(d), ncol=1))[,1]
        
        # Compute prior. 
        lprior_prop <- lprior(par_prop)
 
        # Immediately reject if proposal has prior density zero (which will 
        # often happen when the prior has been truncated to stay within the 
        # design bounds).
        if(is.infinite(lprior_prop)) {
          par_samp[itr,] <- par_samp[itr-1,]
        } else {
          # Sample log-likelihood emulator.
          llik_samp <- sample_mh_llik(llik_em, par_curr, par_prop, 
                                      llik_curr=llik_curr, 
                                      mode=mode, use_joint=use_joint, ...)
          llik_curr <- llik_samp[1]
          llik_prop <- llik_samp[2]

          # Accept-Reject step.
          lpost_curr <- lprior_curr + llik_curr
          lpost_prop <- lprior_prop + llik_prop
          alpha <- min(1.0, exp(lpost_prop - lpost_curr))
          
          if(runif(1) <= alpha) {
            par_samp[itr,] <- par_prop
            par_curr <- par_prop
            lprior_curr <- lprior_prop
            llik_curr <- llik_prop
            accept_count <- accept_count + 1L
          } else {
            par_samp[itr,] <- par_curr
          }
        }
        
        # Adapt proposal covariance matrix and scaling term.
        if(adapt && (((itr-1) %% adapt_interval) == 0)) {
          times_adapted <- times_adapted + 1L
          adapt_list <- adapt_MH_proposal_cov(cov_prop=cov_prop, 
                                              log_scale_prop=log_scale_prop, 
                                              times_adapted=times_adapted, 
                                              adapt_cov=adapt_cov_prop, 
                                              adapt_scale=adapt_scale_prop,
                                              samp_interval=par_samp[(itr-adapt_interval+1):itr,,drop=FALSE], 
                                              accept_rate=accept_count/adapt_interval, 
                                              accept_rate_target=accept_rate_target, 
                                              adapt_factor_exponent=adapt_factor_exponent, 
                                              adapt_factor_numerator=adapt_factor_numerator)
          cov_prop <- adapt_list$cov
          log_scale_prop <- adapt_list$log_scale
          if(adapt_cov_prop) L_cov_prop <- adapt_list$L_cov
          accept_count <- 0L
        }
        
        # Store chain info for current step.
        chain_info[itr,] <- c(llik_curr, lprior_curr, 
                              exp(log_scale_prop) * sqrt(diag(cov_prop)))

        #
        # Gibbs step for sig2. 
        #
        
        if(include_sig2_Gibbs_step) {
          # Update sum of squared residuals (SSR) sample. 
          SSR_curr <- sapply(term_labels_learn_sig2, 
                             function(lbl) emulator_samp_list[[lbl]][curr_prop_idx,])
          
          # TODO: how to deal with parameter ordering here? 
          # SSR_curr, sig2_prior_info, N_obs should all use llik labels. 
          sig2_curr <- sample_NIG_cond_post_sig2(SSR_curr, sig2_prior_info, N_obs)
          sig2_samp[itr,] <- sig2_curr
        }
      }
    },  error = function(cond) {
      err <<- cond
      message("mcmc_gp_noisy() MCMC error; iteration ", itr)
      message(conditionMessage(cond))
    }
  )

  return(list(samp=list(par=par_samp, sig2=sig2_samp),
              info=list(dens=chain_info[,1:2], 
                        prop=chain_info[,3:ncol(chain_info),drop=FALSE]),
              log_scale_prop=log_scale_prop, L_cov_prop=L_cov_prop, 
              par_curr=par_curr, par_prop=par_prop, sig2_curr=sig2_curr,
              itr_curr=itr, par_init=par_init, condition=err))
}


mcmc_gp_unn_post_dens_approx <- function(llik_em, par_prior, approx_type, 
                                         par_init=NULL, sig2_init=NULL, 
                                         sig2_prior=NULL, n_itr=50000L, 
                                         alpha_quantile=0.9, 
                                         cov_prop=NULL, log_scale_prop=NULL, 
                                         adapt_cov_prop=TRUE, adapt_scale_prop=TRUE, 
                                         adapt=adapt_cov_prop||adapt_scale_prop, 
                                         return_prop_sd=FALSE, accept_rate_target=0.24, 
                                         adapt_factor_exponent=0.8, 
                                         adapt_factor_numerator=10, 
                                         adapt_interval=200, ...) {
  # An adaptive Metropolis-Hastings sampler, where the likelihood is given by 
  # a deterministic likelihood approximation defined by a call to 
  # `llik_em$calc_lik_approx()`. Currently this function assumes that 
  # likelihood parameters are fixed, though this will be generalized in the 
  # future.
  #
  # Arguments:
  #    approx_type: character, currently supports "mean", marginal", and 
  #                 "quantile".
  #    alpha_quantile: numeric, number in (0,1), specifying the quantile to use 
  #                    in the quantile approximation. Only used if `approx_type`
  #                    if "quantile".
  #    The remaining arguments are the same as those required by the other 
  #    MCMC functions implemented here.
  #
  # TODO: add check that parameter order in par_prior, llik_em, and optionally 
  #       par_init/cov_prop align.

  # Objects to store samples. 
  d <- llik_em$dim_input
  par_samp <- matrix(nrow=n_itr, ncol=d)
  chain_info <- matrix(nrow=n_itr, ncol=d+2)
  colnames(par_samp) <- llik_em$input_names
  colnames(chain_info) <- c("llik", "lprior", llik_em$input_names)
  
  # At present, this function assumes that the likelihood parameters are fixed.
  learn_sig2 <- FALSE
  sig2_curr <- NULL
  sig2_samp <- NULL

  # Set initial conditions. 
  if(is.null(par_init)) par_init <- sample_prior(par_prior, n=1L)[1,]
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  llik_curr <- llik_em$calc_lik_approx(approx_type, par_curr,
                                       lik_par_val=sig2_curr, 
                                       log_scale=TRUE, 
                                       alpha=alpha_quantile, 
                                       simplify=TRUE, ...)
  lprior <- get_lprior_dens(par_prior)
  lprior_curr <- lprior(par_curr)
  lpost_curr <- lprior_curr + llik_curr
  
  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0L
  times_adapted <- 0L
  
  # Store chain info for initial step.
  chain_info[1,] <- c(llik_curr, lprior_curr, 
                      exp(log_scale_prop) * sqrt(diag(cov_prop)))
  
  # Variable to store error condition, if it occurs.
  err <- NULL
  
  tryCatch(
    {
      for(itr in 2L:n_itr) {
        # Random walk proposal. 
        par_prop <- par_curr + (exp(log_scale_prop) * L_cov_prop %*% matrix(rnorm(d), ncol=1))[,1]
        
        # Compute prior. 
        lprior_prop <- lprior(par_prop)
        
        # Immediately reject if proposal has prior density zero (which will 
        # often happen when the prior has been truncated to stay within the 
        # design bounds). 
        if(is.infinite(lprior_prop)) {
          par_samp[itr,] <- par_samp[itr-1,]
        } else {
          # Compute log-posterior approximation. 
          llik_prop <- llik_em$calc_lik_approx(approx_type, par_prop, 
                                               lik_par_val=sig2_curr, 
                                               log_scale=TRUE, 
                                               alpha=alpha_quantile, 
                                               simplify=TRUE, ...) 
          lpost_prop <- lprior_prop + llik_prop
          
          # Accept-Reject step.
          alpha <- min(1.0, exp(lpost_prop - lpost_curr))
          
          if(runif(1) <= alpha) {
            par_samp[itr,] <- par_prop
            par_curr <- par_prop
            lprior_curr <- lprior_prop
            llik_curr <- llik_prop
            lpost_curr <- lpost_prop
            accept_count <- accept_count + 1L 
          } else {
            par_samp[itr,] <- par_curr
          }
        }
        
        # Adapt proposal covariance matrix and scaling term.
        if(adapt && (((itr-1) %% adapt_interval) == 0)) {
          times_adapted <- times_adapted + 1L
          adapt_list <- adapt_MH_proposal_cov(cov_prop=cov_prop, log_scale_prop=log_scale_prop, 
                                              times_adapted=times_adapted, 
                                              adapt_cov=adapt_cov_prop, 
                                              adapt_scale=adapt_scale_prop,
                                              samp_interval=par_samp[(itr-adapt_interval+1):itr,,drop=FALSE], 
                                              accept_rate=accept_count/adapt_interval, accept_rate_target, 
                                              adapt_factor_exponent, adapt_factor_numerator)
          cov_prop <- adapt_list$cov
          log_scale_prop <- adapt_list$log_scale
          if(adapt_cov_prop) L_cov_prop <- adapt_list$L_cov
          accept_count <- 0
        }
        
        # Store chain info for current step.
        chain_info[itr,] <- c(llik_curr, lprior_curr, 
                              exp(log_scale_prop) * sqrt(diag(cov_prop)))
      }
      
    }, error = function(cond) {
      err <<- cond
      message("mcmc_gp_unn_post_dens_approx() MCMC error; iteration ", itr)
      message(conditionMessage(cond))
    }
  )
  
  return(list(samp=list(par=par_samp, sig2=sig2_samp),
              info=list(dens=chain_info[,1:2], 
                        prop=chain_info[,3:ncol(chain_info),drop=FALSE]),
              log_scale_prop=log_scale_prop, L_cov_prop=L_cov_prop, 
              par_curr=par_curr, par_prop=par_prop, sig2_curr=sig2_curr,
              itr_curr=itr, par_init=par_init, condition=err))
}


mcmc_gp_acc_prob_approx <- function(llik_em, par_prior, par_init=NULL, 
                                    sig2_init=NULL, mode="marginal",
                                    use_joint=TRUE, n_itr=50000L, cov_prop=NULL, 
                                    log_scale_prop=NULL, adapt_cov_prop=TRUE, 
                                    adapt_scale_prop=TRUE, 
                                    adapt=adapt_cov_prop||adapt_scale_prop, 
                                    accept_rate_target=0.24, 
                                    adapt_factor_exponent=0.8, 
                                    adapt_factor_numerator=10, 
                                    adapt_interval=200L, ...) {
  
  # Currently only marginal acceptance probability approximation is supported.
  assert_that(mode %in% "marginal")
  
  # Objects to store samples and auxiliary information.
  d <- llik_em$dim_input
  par_samp <- matrix(nrow=n_itr, ncol=d)
  chain_info <- matrix(nrow=n_itr, ncol=d+3)
  colnames(par_samp) <- llik_em$input_names
  colnames(chain_info) <- c("acc_prob", "lprior", "accept_rate", 
                            llik_em$input_names)
  
  # For algorithm analysis, store GP predictive quantities used in computing 
  # the acceptance probability approximation.
  gp_pred_quantities <- c("mean_curr", "mean_prop", "var_curr", "var_prop",
                          "cov", "acc_prob")
  gp_pred <- matrix(NA, nrow=n_itr, ncol=6L, 
                    dimnames=list(NULL, gp_pred_quantities))

  # At present, this function assumes that the likelihood parameters are fixed.
  learn_sig2 <- FALSE
  sig2_curr <- NULL
  sig2_samp <- NULL
  
  # Set initial conditions. 
  if(is.null(par_init)) par_init <- sample_prior(par_prior, n=1L)[1,]
  par_samp[1,] <- drop(par_init)
  par_curr <- par_samp[1,]
  lprior <- get_lprior_dens(par_prior)
  lprior_curr <- lprior(par_curr)

  # Proposal covariance.
  if(is.null(cov_prop)) cov_prop <- diag(rep(1,d))
  if(is.null(log_scale_prop)) log_scale_prop <- log(2.38) - 0.5*log(d)
  L_cov_prop <- t(chol(cov_prop))
  accept_count <- 0L
  times_adapted <- 0L
  
  # Store information at iteration zero (no accept prob computed yet).
  alpha <- NA
  accept_rate <- NA
  chain_info[1,] <- c(alpha, lprior_curr, accept_rate,
                      exp(log_scale_prop) * sqrt(diag(cov_prop)))
  
  # Variable to store error condition, if it occurs.
  err <- NULL
  
  tryCatch(
    {
      for(itr in 2L:n_itr) {
        # Random walk proposal. 
        par_prop <- par_curr + (exp(log_scale_prop) * L_cov_prop %*% matrix(rnorm(d), ncol=1))[,1]
        
        # Compute prior. 
        lprior_prop <- lprior(par_prop)
        
        # Immediately reject if proposal has prior density zero (which will 
        # often happen when the prior has been truncated to stay within the 
        # design bounds). 
        if(is.infinite(lprior_prop)) {
          par_samp[itr,] <- par_samp[itr-1,]
          alpha <- 0
        } else {
          # Compute acceptance probability approximation.
          alpha_list <- get_gp_mh_acc_prob_approx(par_curr, par_prop, llik_em, 
                                                  mode=mode, use_joint=use_joint,
                                                  lik_par_val=sig2_curr,
                                                  lprior_curr=lprior_curr,
                                                  lprior_prop=lprior_prop, ...)
          
          alpha <- alpha_list$acc_prob
          pred <- alpha_list$llik_pred_list
          gp_pred[itr,] <- c(pred$mean[1], pred$mean[2], pred$var[1], 
                             pred$var[2], pred$cov[1,2], alpha)
          
          # Accept-Reject step.
          if(runif(1) <= alpha) {
            par_samp[itr,] <- par_prop
            par_curr <- par_prop
            lprior_curr <- lprior_prop
            accept_count <- accept_count + 1L 
          } else {
            par_samp[itr,] <- par_curr
          }
        }
        
        # Adapt proposal covariance matrix and scaling term.
        if(adapt && (((itr-1) %% adapt_interval) == 0)) {
          times_adapted <- times_adapted + 1L
          accept_rate <- accept_count/adapt_interval
          adapt_list <- adapt_MH_proposal_cov(cov_prop=cov_prop, log_scale_prop=log_scale_prop, 
                                              times_adapted=times_adapted, 
                                              adapt_cov=adapt_cov_prop, 
                                              adapt_scale=adapt_scale_prop,
                                              samp_interval=par_samp[(itr-adapt_interval+1):itr,,drop=FALSE], 
                                              accept_rate=accept_rate, accept_rate_target, 
                                              adapt_factor_exponent, adapt_factor_numerator)
          cov_prop <- adapt_list$cov
          log_scale_prop <- adapt_list$log_scale
          if(adapt_cov_prop) L_cov_prop <- adapt_list$L_cov
          accept_count <- 0
        }
        
        # Store chain info for current step.
        chain_info[itr,] <- c(alpha, lprior_curr, accept_rate,
                              exp(log_scale_prop) * sqrt(diag(cov_prop)))
      }
      
    }, error = function(cond) {
      err <<- cond
      message("mcmc_gp_acc_prob_approx() MCMC error; iteration ", itr)
      message(conditionMessage(cond))
    }
  )
  
  return(list(samp=list(par=par_samp, sig2=sig2_samp), 
              info=list(dens=chain_info[,1:2], 
                        prop=chain_info[,3:ncol(chain_info),drop=FALSE]),
              gp_pred=gp_pred,
              log_scale_prop=log_scale_prop, L_cov_prop=L_cov_prop, 
              par_curr=par_curr, par_prop=par_prop, sig2_curr=sig2_curr,
              itr_curr=itr, par_init=par_init, condition=err))
}


mcmc_bt_wrapper <- function(llik_em, par_prior, approx_type=NULL, 
                            par_init=NULL, n_itr=50000L, sampler="DEzs", 
                            settings_list=list(), ...) {
  # This function provides a wrapper around the BayesianTools package 
  # `runMCMC()` function, thus providing access to a variety of different 
  # MCMC algorithms. This wrapper is designed so that the argument interface
  # and returned value align with the other MCMC functions in this file.
  # In particular, `mcmc_bt_wrapper()` is set up to run a single chain of 
  # MCMC using a `llik_em` object. Parallelization for multichain runs 
  # should be accomplished using `run_mcmc_chains()`, as for the other 
  # MCMC functions implemented here. The MCMC output from BayesianTools is 
  # transformed to have the same format as the other MCMC functions here (e.g.,
  # `mcmc_noisy_llik()`). Note that the `parallel` element of the BayesianTools
  # settings object refers to within-chain parallelization, for MCMC algorithms
  # for which parallel likelihood evaluations makes sense; it does not control 
  # the running of multiple chains in parallel.
  #
  # This wrapper works by extracting a log-likelihood function from the object
  # `llik_em` via the `get_llik_func()` method. If 
  # `llik_em$exact_llik` is TRUE, then the likelihood will be  
  # exact (not approximated). Otherwise, a deterministic likelihood 
  # approximation is constructed according to the `approx_type` argument.
  # Therefore, `mcmc_bt_wrapper()` provides an interface for both exact and 
  # emulator-accelerated MCMC inference, so long as the approximate likelihood
  # can be evaluated in a pointwise fashion. 
  #
  # Note that BayesianTools appears to require `priorSampler` to be passed to 
  # the bayesianSetup object when the prior is a R function representing a 
  # log prior density. Certain algorithms (e.g., "DEzs") require generating 
  # an initial "population" of samples. In this case, it is recommended to 
  # set `par_init` to NULL and let the BayesianTools code setup the initial 
  # population. 
  #
  # Args:
  #    llik_em: an object for which `is_llik_em(llik_em)` is TRUE.
  #    par_prior: data.frame defining the prior distribution for `par`. 
  #    par_init: numeric, d-dimensional vector, the initial condition. If NULL,
  #              falls back on the default BayesianTools behavior for setting 
  #              the initial condition, which typically implies sampling from 
  #              the prior.
  #    n_itr: integer, the number of MCMC iterations.
  #    sampler: the MCMC algorithm used by `BayesianTools::runMCMC()`. See the
  #             BayesianTools documentation for options. The default "DEzs"
  #             aligns with the BayesianTools default, which refers to a 
  #             specific implementation of a differential evolution MCMC 
  #             algorithm.
  #   settings_list: There are many settings that can be altered in the MCMC 
  #                  algorithms implemented by `runMCMC()`, which are passed
  #                  via the `runMCMC()` "settings" argument. The list provided
  #                  in this argument is passed to this "settings" argument, 
  #                  with a couple caveats: the "n_itr" argument above is 
  #                  used to set `settings_list$iterations` which controls the 
  #                  number of MCMC iterations. Also, `settings_list$nrChains`
  #                  is always set to 1, since this function is designed to run 
  #                  a single MCMC chain. settings_list$startValue` is
  #                  set to `par_init`, which is sampled from the prior if 
  #                  `par_init` is NULL.
  #   ...: additional arguments passed to `llik_em$assemble_llik()`.
  #
  # Returns: 
  # Returns a list with elements "samp" and "bt_output". "samp" follows the 
  # formatting described in e.g. `mcmc_noisy_llik()`. "bt_output" is the object 
  # returned by the `runMCMC()` function. Note that some MCMC algorithms 
  # (e.g., the differential evolution algorithms) run multiple chains in order 
  # to design better proposal distributions. For these algorithms, the returned 
  # sample matrices incorporate samples combined from all of these chains. The 
  # original output from each chain separately can be found in
  # "bt_output$chain". 
  
  # TODO: why does this check fail when executed in parallel?
  # assert_that(is_llik_em(llik_em))

  # Log-likelihood function. Set up to accept numeric vector `par` as an 
  # argument.
  llik <- llik_em$get_llik_func(approx_type=approx_type)

  # Log prior density and sampling function.
  lprior <- get_lprior_dens(par_prior)
  prior_sampler <- get_prior_sampler(par_prior)
  
  # Create BayesianSetup object.
  bayesianSetup <- BayesianTools::createBayesianSetup(likelihood=llik, 
                                                      prior=lprior, 
                                                      parallel=FALSE,
                                                      priorSampler=prior_sampler, 
                                                      names=llik_em$input_names)
  
  # MCMC settings list. Ensure number of MCMC iterations and number of MCMC 
  # chains are in the list; other BayesianTools MCMC settings can optionally be 
  # included.
  settings_list$iterations <- n_itr
  settings_list$nrChains <- 1L
  
  # Set initial condition.
  if(!is.null(par_init) && is.numeric(par_init)) {
    settings_list$startValue <- matrix(par_init, nrow=1L)
  }
  
  # Run MCMC.
  bt_output <- BayesianTools::runMCMC(bayesianSetup=bayesianSetup, 
                                      sampler=sampler, settings=settings_list)
  
  # Convert output into a form that aligns with other MCMC functions. Note that 
  # for certain algorithms (e.g., "DEzs") that return multiple trajectories, 
  # this function will combine the trajectories into a single chain.
  bt_samp <- BayesianTools::getSample(bt_output, parametersOnly=FALSE)

  # Return list with samples stored in matrices as well as the MCMC object 
  # directly returned by BayesianTools.
  par_names <- llik_em$input_names
  dens_info <- bt_samp[,c("Llikelihood", "Lprior"), drop=FALSE]
  colnames(dens_info) <- c("llik", "lprior")
  return(list(samp=list(par=bt_samp[,par_names, drop=FALSE]), 
              info=list(dens=dens_info), 
              bt_output=bt_output, par_init=par_init))
}


# ------------------------------------------------------------------------------
# Proposal Covariance Adaptation 
# ------------------------------------------------------------------------------

adapt_MH_proposal_cov <- function(cov_prop, log_scale_prop, times_adapted, adapt_cov, adapt_scale, 
                                  samp_interval, accept_rate, accept_rate_target=0.24,  
                                  adapt_factor_exponent=0.8, adapt_factor_numerator=10) {
  # Follows the adaptive Metropolis scheme used in the Nimble probabilistic 
  # programming language. 
  # Described here: https://arob5.github.io/blog/2024/06/10/adapt-mh/

  return_list <- list()
  adapt_factor <- 1 / (times_adapted + 3)^adapt_factor_exponent
  
  if(adapt_cov) {
    if(accept_rate > 0) {
      sample_cov_interval <- cov(samp_interval)
      cov_prop <- cov_prop + adapt_factor * (sample_cov_interval - cov_prop)
    } 
    
    # Handle case that `cov_prop` may not be positive definite. 
    L_cov_prop <- try(t(chol(cov_prop)))
    if(inherits(L_cov_prop, "try-error")) {
      message("Problematic proposal cov: ", cov_prop)
      cov_prop <- Matrix::nearPD(cov_prop, ensureSymmetry=TRUE)
      L_cov_prop <- t(chol(cov_prop))
    }
    
    return_list$L_cov <- L_cov_prop
  }
  
  if(adapt_scale) {
    log_scale_factor <- adapt_factor_numerator * adapt_factor * (accept_rate - accept_rate_target)
    log_scale_prop <- log_scale_prop + log_scale_factor
  }
  
  return_list$cov <- cov_prop
  return_list$log_scale <- log_scale_prop
  
  return(return_list)
}

# ------------------------------------------------------------------------------
# MCMC helper/wrapper functions.
#  For running chains in parallel and setting initial conditions.
# ------------------------------------------------------------------------------

sample_mh_llik <- function(llik_em, par_curr, par_prop, llik_curr=NULL, 
                           mode="mcwmh", use_joint=TRUE, ...) {
  # A helper function for `mcmc_noisy` to sample the log-likelihood values 
  # at the current and proposed points for use in computing a 
  # Metropolis-Hastings acceptance ratio.
  #
  # Different possibilities:
  # mcwmh: sample llik at both current and proposed points. Samples jointly
  #        from the bivariate distribution if `joint` is TRUE, otherwise samples
  #        independently.
  # pseudo-marginal: If `use_joint = TRUE`, samples the proposed llik value  
  #                  independently and returns `llik_curr` unchanged for the  
  #                  current value. This targets the "marginal" posterior  
  #                  approximation. If FALSE, `use_joint = TRUE` then samples the 
  #                  proposed llik value conditional on the current 
  #                  log-likelihood value. `llik_curr` is returned unchanged 
  #                  in either case.
  #
  # Args:
  #    llik_em:
  #    par_curr: numeric, the current parameter value (i.e., current MCMC state).
  #    par_prop: numeric, the proposed parameter value.
  #    llik_curr: numeric, the log-likelihood value at the current point. 
  #               Required for all modes except for "mcwmh".
  #    mode: either "mcwmh" or "pseudo-marginal". See description above.
  #    use_joint: logical, see description above.
  
  # The joint pseudo-marginal algorithm requires updating the current GP 
  # emulator, so make a deep copy here to avoid updating in global scope.
  llik_em_copy <- llik_em
  if((mode=="pseudo-marginal") && use_joint) {
    assert_that(!is.null(llik_curr))
    llik_em_copy <- llik_em$copy(shallow=FALSE)
    llik_em_copy$update_emulator(par_curr, llik_curr, update_hyperpar=FALSE)
  }
  
  # mcwmh mode samples llik values at both current and proposed points.
  # pseudo-marg mode samples the llik only at the proposed point, and recycles  
  # the current llik value.
  par <- matrix(par_prop, nrow=1, dimnames=list(NULL, llik_em$input_names))
  if(mode=="mcwmh") par <- rbind(par_curr, par)
  
  llik_samp <- llik_em_copy$sample(par, use_cov=use_joint,
                                   include_nugget=TRUE, ...)
  llik_samp <- drop(llik_samp)
  if(mode=="pseudo-marginal") llik_samp <- c(llik_curr, llik_samp)

  return(llik_samp)
}


run_mcmc_chains <- function(mcmc_func_name, llik_em, n_chain=4L, 
                            par_init=NULL, par_prior=NULL, defer_ic=FALSE,
                            lik_par_init=NULL, lik_par_prior=NULL, 
                            estimate_cov=FALSE, n_samp_cov_est=100L, 
                            try_parallel=TRUE, n_cores=NULL, package_list=NULL, 
                            dll_list=NULL, obj_list=NULL, test_label="default_lbl", 
                            ic_settings=list(), n_itr=50000L, itr_start=1L, ...) {
  # Provides the main interface/entrypoint for running MCMC with llikEmulator
  # objects. This wrapper function runs multiple MCMC chains in parallel. There
  # are three options for specifying initial conditions:
  # 1.) Supply them explicitly via `par_init`.
  # 2.) Generate them via a call to `get_mcmc_ic()`, which by default produces 
  #     a latin hypercube sample with respect to the prior. Supplying additional
  #     arguments to this function via `...` allows the use of different 
  #     initial condition construction schemes.
  # 3.) Set `defer_ic` to TRUE, which simply calls `mcmc_func_name` whether or 
  #     not `par_init` has been specified. This option is designed with 
  #     `mcmc_bt_wrapper()` in mind, which has its own methods for generating 
  #     initial conditions and hence we can "defer" to these methods.
  #
  # This function collects the output from all chains into a single data.table 
  # satisfying the `samp_dt` requirements outlined in `mcmc_helper_functions.r`.
  # Similarly an `info_dt` table is constructed which stores log likelihood and 
  # prior evaluations from MCMC runs, as well as proposal standard deviations.
  # Let d denote the dimension of the parameter space for below.
  #
  # Args:
  #    mcmc_func_name: character, the name of the MCMC function to call.
  #    llik_em: llikEmulator object encoding the exact or approximate 
  #             log-likelihood.
  #    n_chain: integer, number of independent MCMC chains to run.
  #    par_init: matrix of dimension (n_chain,d), each row containing the 
  #              initial condition for one of the MCMC chains.
  #    par_prior: data.frame storing the prior distribution for the parameters.
  #    defer_ic: logical, if FALSE (default), then the initial conditions for
  #              each chain are set up in this function. If TRUE, then defers 
  #              to the underlying MCMC function specified by `mcmc_func_name`.
  #              In this case the initial conditions will be set to NULL in this
  #              function, under the assumption that they will be handled by
  #              `mcmc_func_name`. This argument was motivated by the sampler 
  #              "DEzs" in `mcmc_bt_wrapper()`, which tends to through errors 
  #              when trying to manually set initial conditions.
  #    estimate_cov: If TRUE, estimates the posterior covariance matrix, and 
  #                  uses the estimate to initialize the proposal covariance.
  #                  Currently estimate is produced by simulating from the prior
  #                  and computing an importance sampling estimate. Other methods
  #                  may be implemented in the future; e.g., using estimate
  #                  resulting from ensemble Kalman inversion.
  #   ic_settings: list containing named arguments passed to `get_mcmc_ic()`, 
  #                other than "llik_em", "par_prior", and "n_ic". If `par_init`
  #                is provided then this is not used.
  #   itr_start: integer, passed to `select_mcmc_itr()` to subset the MCMC 
  #              samples. By default returns all samples.
  #   ...: should contain all arguments to be passed to the MCMC function, other 
  #        than the initial conditions and priors.
  #
  # TODO: Need to generalize to handle `lik_par_init`/`lik_par_prior`.

  # If only considering one chain, avoid extra overhead of setting up 
  # cluster for parallel execution.
  if(n_chain==1L) try_parallel <- FALSE
  
  # Sample initial conditions from prior if not provided and if not deferring 
  # to the underlying MCMC function.
  if(!defer_ic) {
    if(is.null(par_init)) {
      ic_settings[c("llik_em","par_prior","n_ic")] <- list(llik_em, par_prior, n_chain)
      par_init <- do.call("get_mcmc_ic", ic_settings)
    } else {
      if(!is.matrix(par_init)) par_init <- matrix(par_init, nrow=1)
      assert_that(nrow(par_init) == n_chain)
    }
    par_init_list <- lapply(1:n_chain, function(i) par_init[i,])
  } else {
    par_init_list <- lapply(1:n_chain, function(i) NULL)
  }
  
  mcmc_chain_list <- tryCatch({
    # Set up MCMC function as a function of the initial condition.
    mcmc_func <- get(mcmc_func_name)
    mcmc_func_ic <- function(ic) mcmc_func(llik_em=llik_em,
                                           par_prior=par_prior,
                                           par_init=ic, n_itr=n_itr, ...)
    
    if(try_parallel) {
      # Prepare for parallel run. 
      if(is.null(n_cores)) {
        n_cores <- parallel::detectCores()
      } else {
        n_cores <- min(n_cores, parallel::detectCores())
      }
      n_cores <- min(n_chain, n_cores)
      cl <- parallel::makeCluster(n_cores)
      
      # Ensure required packages, objects, and DLLs will be available during the 
      # parallel execution.
      export_list <- get_parallel_exports(package_list, dll_list, obj_list)
      parallel::clusterCall(cl, export_list$load_func, export_list$packages, 
                            export_list$dlls)
      parallel::clusterExport(cl, varlist=export_list$objects)
      mcmc_chain_list <- parallel::parLapply(cl=cl, X=par_init_list, fun=mcmc_func_ic)
    } else {
      mcmc_chain_list <- lapply(par_init_list, mcmc_func_ic)
    }
  }, error = function(cond) {
    message("run_mcmc_chains() error:")
    message(conditionMessage(cond))
    err_list <- list(err=cond)
    if(exists("mcmc_chain_list")) err_list$partial_output <- mcmc_chain_list
    return(err_list)
  }, finally = {
    if(try_parallel) try(parallel::stopCluster(cl))
  })
  
  if(isTRUE("err" %in% names(mcmc_chain_list))) return(mcmc_chain_list)

  # Convert to `samp_dt` data.table format
  samp_dt <- format_mcmc_output_multi_chain(lapply(mcmc_chain_list, function(l) l$samp), 
                                            test_label=test_label,
                                            itr_start=itr_start)
  info_dt <- format_mcmc_output_multi_chain(lapply(mcmc_chain_list, function(l) l$info), 
                                            test_label=test_label,
                                            itr_start=itr_start)
  
  # Non-sample output stored in its own list. 
  mcmc_outputs_list <- vector(mode="list", length=length(mcmc_chain_list)) 
  for(i in seq_along(mcmc_outputs_list)) {
    non_samp_elements <- setdiff(names(mcmc_chain_list[[i]]), c("samp","info"))
    mcmc_outputs_list[[i]] <- mcmc_chain_list[[i]][non_samp_elements]
  }
  
  return(list(samp=samp_dt, info=info_dt, output_list=mcmc_outputs_list))
}


run_mcmc <- function(llik_em, par_prior, mcmc_settings) {
  # A higher-level interface to running the MCMC functions that wraps around 
  # `run_mcmc_chains()`. This function requires specifying the log-likelihood 
  # and prior (which are required by all MCMC functions), with all other 
  # arguments stored in a list `mcmc_settings`. This function is primarily 
  # implemented for simulation experiments, where one may want to test many 
  # different MCMC algorithms on the same inverse problem.
  #
  # Args:
  #    llik_em: the llikEmulator object.
  #    par_prior: prior distribution object.
  #    mcmc_settings: list, which must have elements "test_label" and 
  #                   "mcmc_func_name", the former providing an identifier
  #                   for the MCMC execution, and the latter giving the name 
  #                   of the MCMC function to call. All other elements of 
  #                   `mcmc_settings` are interpreted as arguments that will 
  #                   be passed to `run_mcmc_chains()`.
  #
  # Returns:
  # The return value of the call to `run_mcmc_chains()`.
  
  # Ensure MCMC tag and MCMC function name are included in settings.
  required_settings <- "mcmc_func_name"
  missing_settings <- setdiff(required_settings, names(mcmc_settings))
  assert_that(length(missing_settings)==0L, 
              msg=paste0("Missing settings: ", paste0(missing_settings, collapse=", ")))
  assert_that(all(!sapply(mcmc_settings[required_settings], is.null)))
  
  # Assemble full list of arguments for `run_mcmc_chains`.
  mcmc_settings$llik_em <- llik_em
  mcmc_settings$par_prior <- par_prior
  
  # Execute MCMC using `run_mcmc_chains()` wrapper.
  do.call(run_mcmc_chains, args=mcmc_settings)
}


run_mcmc_comparison <- function(llik_em, par_prior, mcmc_settings_list, 
                                save_dir=NULL, return=FALSE, ...) {
  # A wrapper around `run_mcmc()` that runs multiple MCMC algorithms targeting 
  # the same posterior distribution. Each element of `mcmc_settings_list` is 
  # passed to the `mcmc_settings` argument of `run_mcmc()`.
  
  if(return) {
    message("run_mcmc_comparison() currently only supports saving data to file.")
    .NotYetImplemented()
  }
  
  # Store MCMC tags.
  n_algs <- length(mcmc_settings_list)
  tags <- names(mcmc_settings_list)
  if(is.null(tags)) tags <- as.character(1:n_algs)

  for(i in 1:n_algs) {
    # Tag used to identify the algorithm.
    tag <- tags[i]
    print(paste0("Running MCMC alg: ", tag))
    
    tryCatch(
      {
        # Run MCMC using the ith algorithm.
        mcmc_settings_list[[i]]$test_label <- tag
        mcmc_output <- run_mcmc(llik_em, par_prior, mcmc_settings_list[[i]])
        
      }, error=function(cond) {
        message("Error with MCMC run: ", tag)
        message(conditionMessage(cond))
      }, finally = {
        # Save output.
        if(!is.null(out_dir)) {
          filename <- paste0("mcmc_samp_", tag, ".rds")
          saveRDS(mcmc_output, file.path(out_dir, filename))
        }
      }
    )
    
  }
}


get_parallel_exports <- function(package_list=NULL, dll_list=NULL, obj_list=NULL) {
  # In using `parallel` apply functions, objects and packages must explicitly 
  # be made available in the cluster environment. This function returns lists 
  # containing R objects, packages, and dynamic-linked libraries (DLLs) that 
  # can then be passed to `parallel::clusterCall()` and 
  # `parallel::clusterExport()` to make all of these available in the cluster 
  # environment. This function is closely modeled after the function 
  # `generateParallelExecutor()` in the R `BayesianTools` package. By default, 
  # this function will return all loaded packages, DLLs, and objects in the 
  # global environment. This can be wasteful when the function being executed 
  # in parallel only relies on a small subset of these. The arguments of this 
  # function can thus explicitly pass the objects to include when the 
  # dependencies of the function are explicitly known. 
  #

  # Select all loaded packages if not explicitly passed.
  if(is.null(package_list)) package_list <- .packages()
  
  # Select all DLLs if not explicitly passed. Returns object of class 
  # `DLLInfoList`, which can be iterated and DLLs extracted as strings.
  if(is.null(dll_list)) {
    dll_list <- getLoadedDLLs()
    dll_list <- sapply(dll_list, function(x) x[[2]])
  }
  
  # Select all objects in global environment if not explicitly passed.
  if(is.null(obj_list)) obj_list <- ls(envir=.GlobalEnv)

  # Helper function to load packages and DLLs in the parallel cluster.
  load_package_dll <- function(package_list=NULL, dll_list=NULL) {
    # Load packages. 
    if(!is.null(package_list)) {
      for(package_name in package_list) library(package_name, character.only=TRUE)
    }
    
    # Load DLLs.
    if(!is.null(dll_list)) {
      for(dll_name in dll_list) try(dyn.load(dll_name), silent=TRUE)
    }
  }
  
  return_list <- list(packages=package_list, dlls=dll_list, objects=obj_list,
                      load_func=load_package_dll)
  return(return_list)
}


get_mcmc_ic <- function(llik_em, par_prior, n_ic, 
                        n_ic_by_method=setNames(n_ic,"prior"),
                        design_method="LHS", approx_type="mean",  
                        test_info=NULL, n_test_inputs=200L, ...) {
  # Generates a matrix of parameter values intended to be used as initial 
  # conditions for a set of MCMC chains. `n_ic` is the total number of 
  # initial conditions, but different methods can be used to generate 
  # different subsets of the total number of initial conditions. For example, 
  # if `n_ic = 4` then setting 
  # `n_ic_by_method = c(prior=1, design_max=2, approx_max=1)` will generate 
  # one initial condition using the "prior" method, 2 using the "design_max"
  # method and one using the "approx_max" method. These methods are described 
  # below. The latter two methods are especially tailored to the approximate
  # MCMC algorithms that tend to be very multimodal, and for which we 
  # typically hope to initialize the MCMC chains such that they cover the 
  # dominant modes of the distribution.
  #
  # Currently supported methods include:
  #  "prior": generates the ICs via a call to `get_batch_design()`, which 
  #           essentially samples from the prior. The `design_method`
  #           is passed to the "method" argument of `get_batch_design()`. The 
  #           default produces a Latin Hypercube Sample with respect to the
  #           prior `par_prior`.
  #  "design_max": extracts the inputs associated with the largest observed 
  #                true likelihood values. The design inputs and log-likelihood
  #                values are extracted from `llik_em`.
  #  "approx_max": extracts the inputs associated with the largest predicted 
  #                likelihood values at a set of test points. If predictions
  #                have already been made external to this function, then 
  #                they can be provided in a list `test_info` with named 
  #                elements "input" and "llik". Otherwise, `n_test_inputs` are 
  #                sampled from the prior using Latin Hypercube Sampling, then 
  #                `em_llik` is used to predict at these test points using the 
  #                likelihood approximation given by `approx_type`.
  #
  # Returns:
  # matrix, of dimension (n_ic, d), where d is the dimension of the parameter 
  # space. Column names are set to the corresponding parameter names, as given
  # by `llik_em$input_names`.
  
  assert_that(sum(n_ic_by_method)==n_ic)
  methods <- names(n_ic_by_method[n_ic_by_method >= 1])
  par_names <- llik_em$input_names
  ic <- matrix(nrow=0, ncol=length(par_names), dimnames=list(NULL, par_names))

  if("prior" %in% methods) {
    ic_prior <- get_batch_design(design_method, N_batch=n_ic_by_method["prior"], 
                                 prior_params=par_prior, ...)[,par_names]
    ic <- rbind(ic, ic_prior)
  }
  
  if("design_max" %in% methods) {
    inputs_design <- llik_em$get_design_inputs()
    llik_design <- llik_em$get_design_llik()
    llik_obs_order <- order(drop(llik_design), decreasing=TRUE)
    input_sel <- llik_obs_order[1:n_ic_by_method["design_max"]]
    ic_design_max <- inputs_design[input_sel, par_names]
    ic <- rbind(ic, ic_design_max)
  }
  
  if("approx_max" %in% methods) {
    
    # Extract inputs and associated approximate log-likelihood values, if 
    # passed in arguments.
    if(!is.null(test_info)) {
      assert_that(!is.null(test_info$llik) && !is.null(test_info$input))
      inputs <- test_info$input
      llik_approx <- test_info$llik
    } else {
      # Otherwise sample set of test inputs and compute log-likelihood 
      # predictions.
      inputs <- get_batch_design("LHS", N_batch=n_test_inputs, 
                                 prior_params=par_prior)
      llik_approx <- llik_em$calc_lik_approx(approx_type, input=inputs,
                                             simplify=TRUE, log_scale=TRUE, ...)
    }

    # Extract inputs corresponding to largest predicted log-likelihood values.
    llik_approx_order <- order(drop(llik_approx), decreasing=TRUE)
    input_sel <- llik_approx_order[1:n_ic_by_method["approx_max"]]
    ic_approx_max <- inputs[input_sel, par_names]
    ic <- rbind(ic, ic_approx_max)
  }

  return(ic)
}


estimate_post_cov <- function(lpost_dens, par_prior, n_samp) {
  # Produces an estimate of the posterior covariance matrix by sampling from
  # the prior and computing an importance sampling estimate. This requires 
  # the ability to evaluate the unnormalized log posterior density.
  #
  # Args:
  #    lpost_dens: function, the (unnormalized) log posterior density. Assumed
  #                that `lpost_dens` is vectorized so that it can accept a 
  #                matrix of inputs with each row a parameter vector.
  #    par_prior: the data.frame encoding the prior distribution.
  #    n_samp: the number of samples to draw from the prior to compute the 
  #            importance sampling estimate. Note that this is also the number 
  #            of log posterior density evaluations required.
  #
  # Returns:
  #    matrix, the estimated covariance matrix.
  
  .NotYetImplemented()
  
  # # Sample from prior.
  # samp <- sample_prior(par_prior, n=n_samp)
  # 
  # # Evaluate unnormalized log posterior density.
  # lpost <- lpost_dens(samp)
  # 
  # # Compute log of the importance weights.
  # log_w <- lpost_dens(samp) - calc_lprior_dens(samp, par_prior)
}


# ------------------------------------------------------------------------------
# Metropolis-Hastings acceptance probability approximations using 
# log-likelihood emulators. 
# ------------------------------------------------------------------------------

calc_lik_ratio_approx <- function(input_num, input_denom, llik_em, 
                                  lik_par_val=NULL, methods=c("mean"),
                                  log_scale=TRUE, quantile_probs=NULL, 
                                  em_pred_list=NULL, ...) {
  # Computes approximations to the likelihood ratio L(u')/L(u) using a log
  # likelihood emulator. Here `input_num` and `input_denom` are u' and u,
  # respectively. Only one pair of inputs is accepted, but multiple 
  # approximation methods may be provided. The function 
  # `calc_lik_ratio_comparison()` wraps around `calc_lik_ratio_approx()` in 
  # order to compute likelihood ratios over multiple pairs of points, as well 
  # as to include computation of the exact likelihood ratio for comparison.
  # 
  # NOTE: this function is intended for testing purposes. It is not intended 
  #       for use with any of the MCMC functions.
  #
  # Args:
  #    input_num: numeric, or matrix with one row representing `u`.
  #    input_denom: numeric, or matrix with one row representing `u'`.
  #    llik_em: llikEmulator object.
  #    lik_par_val: Likelihood parameter to use. Only a single likelihood 
  #                 parameter is allowed.
  #    method: character vector of likelihood ratio approximation methods to 
  #            compute. See function body for valid options. Setting this 
  #            argument to "all" will use all available approximations.
  #    log_scale: logical, if TRUE returns the log of the likelihood ratio 
  #               approximation.
  #    quantile_prob: only used if `method = "quantile"`. The probability 
  #                   determining which quantile will be used. Can be a 
  #                   vector of probabilities to compute quantile approximations
  #                   using different quantiles.
  #    em_pred_list: optional emulator prediction list, as returned by 
  #                  `llik_em$predict_emulator()`, assumed to correspond to 
  #                   predictions for `rbind(input_num, input_denom)` 
  #                   (in that order).
  
  # TODO: need to add ability to pass quantile prob to `llik_em$calc_lik_approx()`.
  
  input_num <- matrix(input_num, nrow=1)
  input_denom <- matrix(input_denom, nrow=1)
  
  valid_methods <- c("mean", "marginal", "quantile", "joint-marg",
                     "joint-marg-ind")
  if(methods=="all") methods <- valid_methods
  assert_that(all(methods %in% valid_methods))
  require_cov <- ("joint-marg" %in% methods)
  pointwise_methods <- setdiff(methods, c("joint-marg", "joint-marg-ind"))
  
  # If not provided, compute predictive quantities.
  input_pair <- rbind(input_num, input_denom)
  if(is.null(em_pred_list)) {
    colnames(input_pair) <- llik_em$input_names
    em_pred_list <- llik_em$predict_emulator(input_pair, lik_par_val=lik_par_val,
                                             return_mean=TRUE, return_var=TRUE, 
                                             return_cov=require_cov, ...)
  }
  
  # Assemble log-likelihood predictions.
  llik_pred_list <- llik_em$predict(input, lik_par_val=lik_par_val, 
                                    return_mean=TRUE, return_var=TRUE, 
                                    return_cov=require_cov, 
                                    em_pred_list=em_pred_list, ...)
  if(require_cov) assert_that(!is.null(llik_pred_list$cov))
  else llik_pred_list$cov <- NA
  
  # Vector to store results.
  return_vec <- c(llik_num_mean=llik_pred_list$mean[[1]], 
                  llik_denom_mean=llik_pred_list$mean[[2]],
                  llik_num_sd=sqrt(llik_pred_list$var[[1]]),
                  llik_denom_sd=sqrt(llik_pred_list$var[[2]]),
                  llik_cross_cov=llik_pred_list$cov[1,2])

    # `joint-marg` is the only method that requires the cross-covariance.
    # Setting the priors equal here implies we are calculating the likelihood
    # ratio.
  if("joint-marg" %in% methods) {
    l <- gp_acc_ratio_marginal(input_curr=input_denom, 
                               input_prop=input_num, 
                               llik_em, use_joint=TRUE, 
                               lik_par_val=lik_par_val, 
                               llik_pred_list=llik_pred_list, 
                               lprior_curr=1, 
                               lprior_prop=1, 
                               log_scale=log_scale, ...)
    return_vec["joint-marg"] <- l$accept_ratio
  }
    
  if("joint-marg-ind" %in% methods) {
    l <- gp_acc_ratio_marginal(input_curr=input_denom, 
                               input_prop=input_num, 
                               llik_em, use_joint=FALSE, 
                               lik_par_val=lik_par_val, 
                               llik_pred_list=llik_pred_list, 
                               lprior_curr=1, 
                               lprior_prop=1, 
                               log_scale=log_scale, ...)
    return_vec["joint-marg-ind"] <- l$accept_ratio
  } 
    
  # Pointwise methods compute approximations for the numerator and denominator
  # separately, then take the ratio of the approximations.
  if(length(pointwise_methods) > 0) {
    pw_llik_approx <- llik_em$calc_lik_approx(pointwise_methods, 
                                              input=input_pair, 
                                              em_pred_list=em_pred_list,
                                              lik_par_val=lik_par_val, 
                                              log_scale=TRUE, 
                                              return_type="matrix", ...)
    log_lik_approx <- pw_llik_approx[1,] - pw_llik_approx[2,]
    if(!log_scale) log_lik_approx <- exp(log_lik_approx)
    return_vec <- c(return_vec, log_lik_approx)
  }
    
  return(return_vec)
}


calc_lik_ratio_comparison <- function(input_num, input_denom, llik_em,
                                      lik_par_val=NULL, llik_exact=NULL, 
                                      methods="mean", log_scale=TRUE,
                                      quantile_probs=c(0.9), ...) {
  # This function wraps around `calc_lik_ratio_approx()` in order to facililate 
  # computation of approximate likelihood ratios over multiple pairs of points,
  # and to optionally compute the exact likelihood ratio as a basis for 
  # comparison. `input_num` and `input_denom` are matrices of input points, 
  # with each row representing a single input point. The matrices must be the
  # same length. The (approximate) ikelihood ratios are computed for each 
  # pair `input_num[i,]`, `input_denom[i,]`.
  # Only a single `lik_par_val` is allowed. `input_curr` and `input_prop`
  # can both be vectors, and they must be of equal length. `methods`
  # can be "all" to include all tags. See `calc_lik_ratio_approx()` for 
  # specifics. If `llik_exact` is provided then the exact likelihood ratio will 
  # be computed as well. `llik_exact` is currently assumed to be a llikEmulator
  # object.
  # 
  # TODO: will be interesting to run this: (i.) for values sampled from prior; 
  # (ii.) for values sampled from true posterior, and (iii.) for current values
  # from true posterior and proposed values from prior (or for this latter one
  # could consider the hyperrectangle defined by the extent of the exact 
  # posterior samples, and only sample proposed values falling outside of this 
  # hyperrectangle). 

  assert_that(nrow(input_num) == nrow(input_denom))

  # Compute approximate likelihood ratios one input pair at a time.
  results <- NULL
  for(i in 1:nrow(input_num)) {
    results_pair <- calc_lik_ratio_approx(input_num[i,], input_denom[i,], 
                                          llik_em, lik_par_val=lik_par_val, 
                                          methods=methods, log_scale=log_scale, 
                                          quantile_probs=quantile_probs, ...)
    results <- rbind(results, results_pair)
  }
  
  # Optionally compute exact likelihood ratio.
  if(!is.null(llik_exact)) {
    llik_num <- llik_exact$assemble_llik(input_num, 
                                         lik_par_val=lik_par_val, ...)
    llik_denom <- llik_exact$assemble_llik(input_denom, 
                                           lik_par_val=lik_par_val, ...)
    
    llik_diff <- llik_num-llik_denom
    if(!log_scale) llik_diff <- exp(llik_diff)
    results <- cbind(results, llik_exact_num=llik_num)
    results <- cbind(results, llik_exact_denom=llik_denom)
    results <- cbind(results, exact=llik_diff)
  }
  
  rownames(results) <- NULL
  return(results)
}


get_gp_mh_acc_prob_approx <- function(input_curr, input_prop, llik_em, 
                                      mode="marginal", use_joint=TRUE,
                                      lik_par_val=NULL, 
                                      conditional=llik_em$default_conditional, 
                                      normalize=llik_em$default_normalize, 
                                      par_prior=NULL, lprior_curr=NULL, 
                                      lprior_prop=NULL, llik_pred_list=NULL, ...) {
  # TODO: currently assumes symmetric proposal density - can generalize this 
  #       when needed. 
  
  # Compute log-prior evaluations, if not already provided. 
  if(is.null(lprior_curr)) lprior_curr <- calc_lprior_dens(input_curr, par_prior)
  if(is.null(lprior_prop)) lprior_prop <- calc_lprior_dens(input_prop, par_prior)

  if(mode == "mean") {
    .NotYetImplemented()
  } else if(mode == "marginal") {
    alpha_list <- gp_acc_prob_marginal(input_curr, input_prop, llik_em, use_joint, 
                                       lik_par_val=lik_par_val, 
                                       conditional=conditional, 
                                       normalize=normalize, 
                                       llik_pred_list=llik_pred_list, 
                                       lprior_curr=lprior_curr, 
                                       lprior_prop=lprior_prop, ...)
  } else {
    stop("Invalid `approx_type` ", alpha_approx_type)
  }
  
  return(alpha_list)
}


gp_acc_ratio_marginal <- function(input_curr, input_prop, llik_em, 
                                  use_joint=TRUE, lik_par_val=NULL, 
                                  conditional=llik_em$default_conditional, 
                                  normalize=llik_em$default_normalize, 
                                  llik_pred_list=NULL, par_prior=NULL, 
                                  lprior_curr=NULL, lprior_prop=NULL, 
                                  log_scale=FALSE, ...) {
  # For a Metropolis-Hastings acceptance probability of the form 
  # min{1, r(input_curr, input_prop)}, this function returns 
  # log E[r(input_curr, input_prop)], where the expectation is with respect to 
  # the `llik_em` distribution. Currently this function only works when the 
  # emulator distribution implies that r(input_curr, input_prop) ~ LN(m, s^2).
  # Currently assumes a symmetric proposal distribution, as the proposal 
  # density is not considered here.
  
  # Compute log-prior evaluations, if not already provided. 
  if(is.null(lprior_curr)) lprior_curr <- calc_lprior_dens(input_curr, par_prior)
  if(is.null(lprior_prop)) lprior_prop <- calc_lprior_dens(input_prop, par_prior)
  
  # Compute predictive distribution at inputs.   
  if(is.null(llik_pred_list)) {
    input <- rbind(input_curr, input_prop)
    rownames(input) <- llik_em$input_names
    llik_pred_list <- llik_em$predict(input, lik_par_val=lik_par_val, 
                                      return_mean=TRUE, return_var=TRUE, 
                                      return_cov=use_joint, 
                                      conditional=conditional, 
                                      normalize=normalize, ...)
  } else {
    assert_that(!is.null(llik_pred_list$mean) && !is.null(llik_pred_list$var))
    if(use_joint) assert_that(!is.null(llik_pred_list$cov))
  }
  
  # Predictive mean and variance of the log acceptance ratio. This is the m and 
  # s^2 in r ~ LN(m, s^2). The variance depends on whether or not the joint 
  # dist is used. 
  # TODO: update the cov indexing below once this is changed. 
  m <- lprior_prop - lprior_curr + llik_pred_list$mean[2] - llik_pred_list$mean[1]
  s2 <- sum(llik_pred_list$var)
  if(use_joint) s2 <- s2 - 2*llik_pred_list$cov[1,2]

  # Marginal acceptance ratio approximation.
  log_ratio <- m + 0.5*s2
  return_list <- list(llik_pred_list=llik_pred_list)
  return_list$accept_ratio <- ifelse(log_scale, log_ratio, exp(log_ratio))
  
  return(return_list)
}


gp_acc_prob_marginal <- function(input_curr, input_prop, llik_em,  
                                 use_joint=TRUE, lik_par_val=NULL, 
                                 conditional=llik_em$default_conditional, 
                                 normalize=llik_em$default_normalize, 
                                 llik_pred_list=NULL, par_prior=NULL, 
                                 lprior_curr=NULL, lprior_prop, ...) {
  # For a Metropolis-Hastings acceptance probability of the form 
  # alpha(input_curr, input_prop) min{1, r(input_curr, input_prop)}, 
  # this function returns E[alpha(input_curr, input_prop)], where the 
  # expectation is with respect to  the `llik_em` distribution.
  
  # Compute log-prior evaluations, if not already provided. 
  if(is.null(lprior_curr)) lprior_curr <- calc_lprior_dens(input_curr, par_prior)
  if(is.null(lprior_prop)) lprior_prop <- calc_lprior_dens(input_prop, par_prior)
  
  # Compute predictive distribution at inputs.   
  if(is.null(llik_pred_list)) {
    input <- rbind(input_curr, input_prop)
    colnames(input) <- llik_em$input_names
    llik_pred_list <- llik_em$predict(input, lik_par_val=lik_par_val, 
                                      return_mean=TRUE, return_var=TRUE, 
                                      return_cov=use_joint, 
                                      conditional=conditional, 
                                      normalize=normalize, ...)
  } else {
    assert_that(!is.null(llik_pred_list$mean) && !is.null(llik_pred_list$var))
    if(use_joint) assert_that(!is.null(llik_pred_list$cov))
  }
  
  # Predictive mean and variance of the log acceptance ratio. This is the  
  # m and s^2 in r ~ LN(m, s^2). The variance depends on whether or not the 
  # joint dist is used.
  # TODO: update the cov indexing below once this is changed. 
  m <- lprior_prop - lprior_curr + llik_pred_list$mean[2] - llik_pred_list$mean[1]
  s2 <- sum(llik_pred_list$var)
  if(use_joint) s2 <- s2 - 2*llik_pred_list$cov[1,2]
  
  # Marginal acceptance probability is a linear combination of 1 and the 
  # marginal ratio approximation.
  log_acc_ratio_marginal <- m + 0.5*s2
  lw1 <- log(pnorm(m/sqrt(s2)))
  lw2 <- log(pnorm(-(m+s2)/sqrt(s2)))
  acc_prob_marginal <- exp(matrixStats::logSumExp(c(lw1, lw2+log_acc_ratio_marginal)))
  
  return(list(acc_prob=acc_prob_marginal, llik_pred_list=llik_pred_list))
}


# -----------------------------------------------------------------------------
# Sample-Based Posterior Approximation: 
# The algorithms here attempt to approximately sample the so-called "sample-
# based" approximate posterior. Samples from this distribution are generated 
# by first sampling a random distribution (via sampling a GP trajectory) and 
# then sampling a parameter value from this distribution. In general, this is 
# impractical due to the requirement of sampling a random probability measure. 
# The functions here attempt to approximate this sampling procedure. 
# -----------------------------------------------------------------------------


sim_sample_based_approx_1d_grid <- function(llik_em, par_prior_params, bounds, 
                                            N_grid, N_samp, lik_par_val=NULL, 
                                            max_rejects=1000L, ...) {
  # This function (approximately) generates samples via: 
  #    (1) Sample a random log-likelihood L. 
  #    (2) Sample a value u from the posterior implied by L. 
  # To make this practical, L is sampled at a finite number of locations (an equally
  # spaced grid) and then linearly interpolated to produce an approximation to L^hat(u).
  # Then u is sampled from the distribution with unnormalized density 
  # prior(u)*exp{L^hat(u)} via rejection sampling. 
  
  assert_that(llik_em$dim_input==1L)
  
  # Construct grid. 
  par_grid <- seq(bounds[1], bounds[2], length.out=N_grid)
  par_grid_mat <- matrix(par_grid,ncol=1)
  colnames(par_grid_mat) <- llik_em$input_names
  
  # Iterate until desired number of samples is achieved. 
  par_samp <- vector(mode="numeric", length=N_samp)

  for(i in 1:N_samp) {
    
    # Sample log-likelihood values over dense grid. 
    llik_grid <- drop(llik_em$sample(par_grid_mat, lik_par_val=lik_par_val, N_samp=1, ...))
    lpost_grid <- llik_grid + calc_lprior_theta(par_grid_mat, par_prior_params)
    
    # Linearly interpolate to obtain approximate log-likelihood sample/trajectory. 
    # llik_interp <- get_linear_interpolation_1d(par_grid, llik_grid)
    # lpost_interp <- function(par) llik_interp(par) + calc_lprior_theta(par, par_prior_params)
    lik_interp <- get_linear_interpolation_1d(par_grid, exp(llik_grid))
    post_interp <- function(par) lik_interp(par) * exp(calc_lprior_theta(par, par_prior_params))
     
    # Upper bound for (unnormalized) density. 
    M <- max(exp(lpost_grid))
    assert_that(M > 0)
    assert_that(!is.infinite(M))
    
    # Rejection sampling until a sample is obtained. 
    N_rejects <- 0L 
    accepted <- FALSE
    while(!accepted && (N_rejects < max_rejects)) {
      par_prop <- runif(n=1, min=bounds[1], max=bounds[2])
      u <- runif(n=1, min=0, max=M)
      if(u <= post_interp(par_prop)) {
      # if(log(u) <= lpost_interp(par_prop)) {
        accepted <- TRUE
        par_samp[i] <- par_prop
      } else {
        N_rejects <- N_rejects+1L
        if(N_rejects==max_rejects) stop(paste0("Reached max rejects when attempting to simulate sample ", i))
      }
    }
  }
  
  return(par_samp)
}

estimate_sample_based_density_1d_grid <- function(llik_em, par_prior_params, input_grid, 
                                                  N_monte_carlo=1000L, lik_par_val=NULL, ...) {
  # This function is similar to `sim_sample_based_approx_1d_grid()` but its goal is not 
  # to draw parameter samples `u`, but instead to estimate the expectation of the random 
  # density of `u` implied by the algorithm 
  #    (1) Sample a random log-likelihood L. 
  #    (2) Sample a value u from the posterior implied by L.
  # This function estimates the density of the samples `u` returned by this algorithm. 
  # It does so by sampling density values on a grid, numerically normalizing each,
  # and then computing the average over the normalized values for each grid point. 
  # Note that this (approximately) computes the expected normalized density, whereas
  # the marginal approximation computes the expected unnormalized density and then 
  # normalizes after-the-fact. This approach has the benefit over 
  # `sim_sample_based_approx_1d_grid()` of not requiring rejection sampling, but
  # requires more care regarding numerical stability. I am completely ignoring this 
  # for now and working on the original (non-log) scale. If this function is to be 
  # used in cases where log-likelihood values might be very small, then this code 
  # should be updated. 
  
  assert_that(llik_em$dim_input==1L)
  
  # Simulate log-likelihood values. Return shape is (N_grid, N_monte_carlo).
  llik_samp <- llik_em$sample(input_grid, lik_par_val=lik_par_val, N_samp=N_monte_carlo, ...)

  # Assuming numbers are reasonable for now so exponentiating is fine. 
  lprior_dens <- calc_lprior_theta(input_grid, par_prior_params)
  post_dens_samp <- exp(add_vec_to_mat_cols(lprior_dens, llik_samp))
  
  # Normalize.
  dx <- abs(drop(input_grid)[2] - drop(input_grid)[1])
  for(j in 1:ncol(post_dens_samp)) {
    post_dens_samp[,j] <- post_dens_samp[,j] / int_trap(post_dens_samp[,j], dx)
  }
  
  # Average over Monte Carlo samples for each input.
  density_estimates <- rowMeans(post_dens_samp)
  
  return(density_estimates)
}


get_linear_interpolation_1d <- function(input_points, output_points) {
  
  point_order <- order(input_points)
  input_points <- input_points[point_order]
  output_points <- output_points[point_order]
  
  f <- function(x) {
    idx_l <- max(which(input_points <= x))
    idx_u <- idx_l+1L
    x_l <- input_points[idx_l]
    x_u <- input_points[idx_u]
    y_l <- output_points[idx_l]
    y_u <- output_points[idx_u]
    
    y_l + (y_u-y_l)/(x_u-x_l)*(x-x_l)
  }
  
  return(f)
  
}


# ---------------------------------------------------------------------
# Helper functions 
# ---------------------------------------------------------------------

run_llik_em_samplers <- function(sampler_settings_list, llik_em_default=NULL, 
                                       par_prior=NULL,  lik_par_prior_default=NULL, N_mcmc_default=NULL) {
  # Need to think about how to implement this. I want the ability to easily specify common settings that 
  # are shared across all algorithms, but also want the ability to have complete control over the 
  # settings for each algorithm, if desired. I'm currently thinking that there should be arguments 
  # to specify "defaults", which will be overwritten if they are also specified in "sampler_settings_list". 
  
  .NotYetImplemented()
}











