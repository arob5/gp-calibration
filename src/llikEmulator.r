#
# gp.r
# Class definitions for the `llikEmulator`, `llikSumEmulator`, and classes which 
# inherit from these classes. All classes defined using reference classes. 
#
# Andrew Roberts
# 

library(assertthat)
library(abind)


# -----------------------------------------------------------------------------
# llikEmulator: Class encapsulating a (typically stochastic) approximation to 
#               the log-likelihood. 
# 
# This is the base class for a log likelihood approximation. It is intentionally 
# quite general and intended to be the parent class for more specialized 
# log likelihood approximation classes (e.g. multiplicative Gaussian llik 
# with sum of squares emulated by GPs). This class is primarily intended 
# to encapsulate a stochastic approximation to the log-likelihood, but 
# it can be set up to function with a deterministic approximation (or 
# just the exact likelihood, which is useful for testing code). Given the 
# generality, not all of the core methods may apply to all approximate 
# likelihoods. For example, the `mean_lik()` method is intended to return 
# the expectation of the likelihood (i.e. the expectation of the 
# exponentiated log-likelihood). In certain cases, there may not be a closed
# form for this expression and the user may just implement this method 
# to return an error if it is called. Alternatively, the user may want to 
# implement this function to return a numerical approximation of this 
# quantity. The field `model` is intended to be a model object of some 
# sort (e.g. a gpWrapper instance) which is used to define the llik 
# approximation. `lik_par` is stores the likelihood parameters. These 
# can be stored in advance if the likelihood parameters will be fixed
# throughout an analysis; otherwise, this field can be left NULL. Note 
# that the likelihood approximation may be a function of the likelihood 
# parameters (e.g. the likelihood parameters are included as inputs to 
# a GP emulator). 
#
# This class organizes the pieces comprising the likelihood into 
# (1) "inputs": the primary parameters of interest; e.g. the parameters 
# characterizing the forward model in an inverse problem; (2) "lik_par", 
# other likelihood parameters (e.g. variance parameters) that define 
# the observation model, may be unknown, but are not of primary interest; 
# and (3) "lik_components", which are fixed quantities such as a 
# sample size. 
# 
# Implementing a deterministic approximation: 
# There are different ways to approach this. One is to just have the 
# methods `mean_log()`, `mean_lik()`, `sample()` all return the 
# determistic approximation. It may make sense to implement 
# `var_log()` and `var_lik()` to return 0 or throw an error, depending 
# on the specific needs. Note that this same method could be used 
# to implement the exact (not approximate) likelihood as well, which 
# could then be used for testing purposes. 
# -----------------------------------------------------------------------------

#
# TODO: update this so that it is defined by default to be a CONDITIONAL likelihood, 
# viewed as a function of `input` with `lik_par` being the other parameters being 
# conditioned on, and `lik_component` being other fixed stuctural information (such 
# as sample size, number of outputs, etc.). Thus, by default all functions should 
# output unnormalized values, dropping terms that only depend on `lik_par`. Overriding 
# the default would then include these terms. An alternative is to have a class 
# attribute that stores the default behavior.

llikEmulator <- setRefClass(
  Class = "llikEmulator", 
  fields = list(emulator_model="ANY", lik_description="character", 
                emulator_description="character", default_conditional="logical",
                default_normalize="logical", lik_par_fixed="logical", lik_par="ANY", 
                input_names="character", lik_par_names="character")
)

llikEmulator$methods(
  
  initialize = function(lik_description, emulator_description, emulator_model=NULL, 
                        default_conditional=FALSE, default_normalize=FALSE, 
                        lik_par_fixed=FALSE, lik_par=NULL, ...) {
    
    if(lik_par_fixed && !is.null(lik_par)) {
      message("`lik_par_fixed` currently does not formally lock the variable - it can be modified.")
      # .self$lock("lik_par")
    }
    
    initFields(lik_description=lik_description, emulator_description=emulator_description,
               emulator_model=emulator_model, default_conditional=default_conditional,
               default_normalize=default_normalize, lik_par_fixed=lik_par_fixed, lik_par=lik_par)
    
    if(lik_par_fixed && is.null(lik_par)) {
      message("Fixed `lik_par` not passed. May be set once and then field will be locked.")
    }
  }, 
  
  get_lik_par = function (lik_par_val=NULL) {
    
    if(lik_par_fixed) return(lik_par)
    return(lik_par_val)
    
  },
  
  sample = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                    normalize=default_normalize, return_list=FALSE, ...) {
    stop("`sample()` method not implemented.")
  },
  
  mean_log = function(input, lik_par_val=NULL, conditional=default_conditional, 
                      normalize=default_normalize, ...) {
    stop("`mean_log()` method not implemented.")
  }, 
  
  var_log = function(input, lik_par_val=NULL, conditional=default_conditional, 
                     normalize=default_normalize, ...) {
    stop("`var_log()` method not implemented.")
  }, 
  
  mean_lik = function(input, lik_par_val=NULL, conditional=default_conditional, 
                      normalize=default_normalize, ...) {
    stop("`mean_lik()` method not implemented.")
  }, 
  
  var_lik = function(input, lik_par_val=NULL, conditional=default_conditional, 
                     normalize=default_normalize, ...) {
    stop("`var_lik()` method not implemented.")
  },
  
  log_norm_constant = function(lik_par_val=NULL, conditional=default_conditional, 
                               normalize=default_normalize, ...) {
    stop("`log_norm_constant()` method not implemented.")
  }
  
)


# -----------------------------------------------------------------------------
# llikSumEmulator Class 
# -----------------------------------------------------------------------------

llikSumEmulator <- setRefClass(
  Class = "llikSumEmulator", 
  contains = "llikEmulator",
  fields = list(llik_emulator_terms="list", term_labels="character", N_terms="integer")
)

llikSumEmulator$methods(
  
  initialize = function(llik_emulator_list, default_conditional=FALSE, default_normalize=FALSE, 
                        lik_description="Sum of llik terms.", emulator_description="", 
                        llik_emulator_labels=NULL, ...) {
    
    assert_that(all(sapply(testlist, function(obj) inherits(obj, "llikEmulator"))), 
                msg="`llikEmulator_list` must be a list of llikEmulator objects.")
    assert_that(!is.null(llik_emulator_labels) || !is.null(names(llik_emulator_list)), 
                msg="`llik_emulator_list` must be named list OR `llik_emulator_labels` must be passed.")
    if(!is.null(llik_emulator_labels)) {
      assert_that(is.character(llik_emulator_labels) && (length(llik_emulator_labels)==length(llik_emulator_list)), 
                  msg="`llik_emulator_labels` must be character vector of equal length to `llik_emulator_list`.")
    }
    
    if(is.null(llik_emulator_labels)) llik_emulator_labels <- names(llik_emulator_list)
    initFields(llik_emulator_terms=llik_emulator_list, term_labels=llik_emulator_labels, 
               N_terms=length(llik_emulator_list))
    
    callSuper(lik_description=lik_description, emulator_description=emulator_description,
              emulator_model=NULL, default_conditional=default_conditional, 
              default_normalize=default_normalize, 
              lik_par_fixed=all(sapply(llik_emulator_list, function(x) x$lik_par_fixed)), lik_par=NULL)
  },
  
  sample = function(input, lik_par=NULL, N_samp=1, conditional=default_conditional, 
                    normalize=default_normalize, sum_terms=TRUE, labels=term_labels, 
                    return_list=FALSE, ...) {
    
    samp <- array(dim=c(nrow(input), N_samp, length(term_labels)))
    samp_emulator <- list()
    
    for(i in seq_along(labels)) {
      lbl <- labels[i]
      samp_term <- llik_emulator_terms[[lbl]]$sample(input, lik_par[[lbl]], N_samp=N_samp, 
                                                     conditional=conditional, normalize=normalize, 
                                                     return_list=TRUE)
      samp[,,i] <- samp_term$llik
      if(return_list) samp_emulator[[lbl]] <- samp_term$emulator
    }
    
    if(sum_terms) samp <- rowSums(samp, dims=2)
    if(return_list) return(list(llik=samp, emulator=samp_emulator))
    return(samp)
  }
  
)


# -----------------------------------------------------------------------------
# llikEmulatorMultGausGP class
# -----------------------------------------------------------------------------

llikEmulatorMultGausGP <- setRefClass(
  Class = "llikEmulatorMultGausGP", 
  contains = "llikEmulator",
  fields = list(N_obs="integer")
)

llikEmulatorMultGausGP$methods(
  
  initialize = function(gp_model, N_obs, sig2=NULL, default_conditional=FALSE, default_normalize=FALSE, 
                        lik_par_fixed=FALSE, ...) {
    assert_that(inherits(gp_model, "gpWrapper"), msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(is.integer(N_obs) && (length(N_obs) == 1) && (N_obs > 0),
                msg="`N_obs` must be an integer greater than 0.")
    assert_that(gp_model$Y_dim==1, msg="`llikEmulatorMultGausGP` only supports single-output GP emulator.")
    
    if(!is.null(sig2)) {
      assert_that(is.numeric(sig2) && (sig2>0), msg="`sig2` must be NULL or numeric positive value.")
    }
    
    initFields(N_obs=N_obs)
    callSuper(lik_description="Multiplicative Gaussian.",
              emulator_description="GP emulating sum of squared error function.",
              default_conditional=default_conditional, default_normalize=default_normalize,
              emulator_model=gp_model, lik_par=sig2, lik_par_fixed=lik_par_fixed, ...)
  }, 
  
  assemble_llik = function(SSR, sig2=NULL, conditional=default_conditional, normalize=default_normalize) {
    # SSR should be N_input x N_samp. 
    
    # Fetch the variance parameter. 
    sig2 <- get_lik_par(sig2)
    
    # Construct likelihood using the quadratic error. 
    llik <- -0.5 * SSR / sig2
    if(normalize || !conditional) llik <- llik - 0.5*N_obs*log(sig2)
    if(normalize) llik <- llik - 0.5*N_obs*log(2*pi)
    
    return(llik)
    
  },
  
  sample = function(input, sig2=NULL, N_samp=1, use_cov=FALSE, include_nugget=TRUE,
                    conditional=default_conditional, normalize=default_normalize, 
                    return_list=FALSE, ...) {
    # Sample quadratic error.  
    samp <- emulator_model$sample(input, use_cov=use_cov, include_nugget=include_nugget, N_samp=N_samp)[,,1]
    
    # Set negative error samples to 0. 
    if(any(samp < 0)) {
      message("Warning: Setting negative GP quadratic error samples to 0.")
      samp[samp < 0] <- 0
    }
    
    # Compute unnormalized or normalized log-likelihood. 
    llik_samp <- assemble_llik(samp, sig2, conditional, normalize)
    
    # If `return_list` is TRUE, include emulator model samples as well. 
    if(return_list) return(list(llik=llik_samp, emulator=samp))
    return(llik_samp)

  }

)



# -----------------------------------------------------------------------------
# llikEmulatorMultGausGP: Encapsulates an approximation of a multiplicative 
# Gaussian likelihood where the "sum of squared error" functions have been 
# approximated by independent GPs. In particular, considers log-likelihood 
# of the form: 
#   log p(Y|u, Sig) = C - 0.5 * sum_{p=1}^{P} sum_{t=1}^{T_p} (Y_{tp} - G(u)_{tp})^2 / sig2_p
# where C is a constant and G a function. The mappings 
#   Phi(u) := sum_{t=1}^{T_p} (Y_{tp} - G(u)_{tp})^2
# are assumed to be emulated by independent GPs. Hence the `model` field of this 
# class is required to inherit from the `gpWrapper` class. This class interprets
# the `input` argument as `u` (typically the forward model parameters, which 
# are of primary interest). This class also requires a second argument `sig2` to its 
# core functions, which is the vector of variance parameters. 
#
# By default, the class represents a "conditional" llik (meaning conditional 
# on some current value of `sig2`), which practically means that setting 
# `normalize=FALSE` will drop additive constants that depend on `sig2`. 
# By setting `default_conditional=FALSE`, this behavior is reversed and 
# terms involving `sig2` are no longer interpreted as normalizing constants. 
# -----------------------------------------------------------------------------

#
# TODO: should probably also have method for computing the exact likelihood. i.e., 
#       pass in either the sufficient statistic or model output and it will compute. 
#

llikSumEmulatorMultGausGP <- setRefClass(
  Class = "llikSumEmulatorMultGausGP", 
  contains = "llikEmulator",
  fields = list(N_obs="integer", N_output="integer", sum_obs="integer")
)

llikSumEmulatorMultGausGP$methods(
  
  initialize = function(gp_model, N_obs, N_output, sig2=NULL, default_conditional=TRUE, 
                        default_normalize=FALSE, ...) {
    assert_that(inherits(gp_model, "gpWrapper"), msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(is.integer(N_output) && N_output>0, msg="`N_output` must be an integer greater than 0.")
    assert_that(is.integer(N_obs) && (length(N_obs)==N_output), 
                msg="`N_obs` must be integer vector of length equal to `N_output`.")
    assert_that(gp_model$Y_dim==N_output, msg="Number of independent GP emulators must equal `N_output`.")
    if(!is.null(sig2)) {
      assert_that(is.numeric(sig2) && (length(sig2)==N_output) && all(sig2>0), 
                  msg="`sig2` must be NULL or numeric vector of length `N_output`.")
    }
    
    initFields(N_obs=N_obs, N_output=N_output, default_conditional=default_conditional, 
               default_normalize=default_normalize, sum_obs=sum(N_obs))
    callSuper(lik_description="Multiplicative Gaussian.",
              emulator_description="Independent GPs emulating sum of squared error functions.",
              emulator_model=gp_model, ...)
  },
  
  
  sample = function(input, sig2, N_samp=1, use_cov=FALSE, include_nugget=TRUE, sum_output_llik=TRUE,
                    conditional=default_conditional, normalize=FALSE, ...) {
    samp <- emulator_model$sample(input, use_cov=use_cov, include_nugget=include_nugget, N_samp=N_samp)
    
    for(j in 1:N_output) {
      samp[,,j] <- -0.5 * samp[,,j] / sig2[j]
      if(normalize || !conditional) samp[,,j] <- samp[,,j] - 0.5*N_obs[j]*log(sig2[j])
    }
    
    if(normalize) samp <- samp - 0.5*sum_obs*log(2*pi)
    if(sum_output_llik) return(rowSums(samp, dims=2))
    return(samp)
  }
  
)

