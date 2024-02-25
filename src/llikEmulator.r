#
# gp.r
# Class definitions for the `llikEmulator`, `llikSumEmulator`, and classes which 
# inherit from these classes. All classes defined using reference classes. 
#
# Andrew Roberts
# 

library(assertthat)
library(ggmatplot)
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
  fields = list(emulator_model="ANY", lik_description="character", llik_label="character",
                emulator_description="character", default_conditional="logical",
                default_normalize="logical", use_fixed_lik_par="logical", lik_par="ANY",
                input_names="character", dim_input="integer")
                 
)

llikEmulator$lock("llik_label")

llikEmulator$methods(
  
  initialize = function(llik_label, input_names, lik_description, emulator_description, dim_input,  
                        emulator_model=NULL, default_conditional=FALSE, default_normalize=FALSE,
                        use_fixed_lik_par=FALSE, lik_par=NULL, ...) {

    initFields(llik_label=llik_label, input_names=input_names, lik_description=lik_description, 
               dim_input=dim_input, emulator_description=emulator_description,
               emulator_model=emulator_model, default_conditional=default_conditional,
               default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, lik_par=lik_par)  
        
    if(use_fixed_lik_par && is.null(lik_par)) message("Fixed `lik_par` not yet passed.")
  }, 
  
  get_lik_par = function(lik_par_val=NULL, ...) {
    if(use_fixed_lik_par) return(lik_par)
    assert_that(!is.null(lik_par_val), msg="`lik_par_val` arg and `lik_par` field are both NULL.")
    return(lik_par_val)
  },
  
  get_input = function(input, ...) {
    assert_that(is.matrix(input) && (ncol(input)==dim_input), 
                msg="`input` must be matrix with ncol equal to `dim_input`.")
    
    return(input[,input_names, drop=FALSE])
  },
  
  sample_emulator = function(input, N_samp=1, ...) {
    .NotYetImplemented()
  },
  
  assemble_llik = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                           normalize=default_normalize, ...) {
    .NotYetImplemented()
  },
  
  sample = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                    normalize=default_normalize, ...) {
    .NotYetImplemented()
  },
  
  plot_llik_samp_1d = function(input_new, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                               normalize=default_normalize, true_llik_new=NULL, ...) {
    
    assert_that(input_dim==1, msg=paste0("plot_llik_samp_1d() requires 1d input space. input_dim = ", input_dim))
    
    input_new <- get_input(input_new)
    llik_samp <- .self$sample(input_new, lik_par=lik_par_val, N_samp=N_samp, ...)
    
    matplot(input_new, llik_samp, type="l", col="gray", main="llik Samples", xlab=input_names, 
            ylab=paste0("Log Likelihood: ", llik_label))
    if(!is.null(true_llik_new)) matlines(input_new, true_llik_new, col="red")

    return(plts)
    
  },
  
  plot_llik_samp = function(input, lik_par_val=NULL, N_samp=1, conditional=default_conditional, 
                            normalize=default_normalize, ...) {
    llik_samp <- .self$sample(input, lik_par_val=lik_par_val, N_samp=N_samp, 
                              conditional=conditional, normalize=normalize, ...)
  },
  
  mean_log = function(input, lik_par_val=NULL, conditional=default_conditional, 
                      normalize=default_normalize, ...) {
    .NotYetImplemented()
  }, 
  
  var_log = function(input, lik_par_val=NULL, conditional=default_conditional, 
                     normalize=default_normalize, ...) {
    .NotYetImplemented()
  }, 
  
  mean_lik = function(input, lik_par_val=NULL, conditional=default_conditional, 
                      normalize=default_normalize, ...) {
    .NotYetImplemented()
  }, 
  
  var_lik = function(input, lik_par_val=NULL, conditional=default_conditional, 
                     normalize=default_normalize, ...) {
    .NotYetImplemented()
  },
  
  log_norm_constant = function(lik_par_val=NULL, conditional=default_conditional, 
                               normalize=default_normalize, ...) {
    .NotYetImplemented()
  }
  
)


# -----------------------------------------------------------------------------
# llikSumEmulator Class 
# -----------------------------------------------------------------------------

llikSumEmulator <- setRefClass(
  Class = "llikSumEmulator", 
  contains = "llikEmulator",
  fields = list(llik_emulator_terms="list", N_terms="integer")
)

llikSumEmulator$methods(
  
  initialize = function(llik_emulator_list, default_conditional=FALSE, default_normalize=FALSE, 
                        lik_description="Sum of llik terms.", emulator_description="", ...) {
    
    assert_that(all(sapply(llik_emulator_list, function(obj) inherits(obj, "llikEmulator"))), 
                msg="`llikEmulator_list` must be a list of llikEmulator objects.")
    assert_that(length(unique(sapply(llik_emulator_list, function(obj) obj$dim_input)))==1, 
                msg="Inputs of llik terms have differing dimensions.")
    
    # Set labels for each llik emulator term, ensuring they are unique. The names attribute of 
    # the list in the field `llik_emulator_terms` is also assigned the vector of labels, ensuring
    # the order of the two is the same. 
    term_lbls <- sapply(llik_emulator_list, function(obj) obj$llik_label)
    assert_that(length(unique(term_lbls))==length(llik_emulator_list), 
                msg="Found duplicate `llik_label` attributes in `llik_emulator_list`.")
    names(llik_emulator_list) <- term_lbls
    initFields(llik_emulator_terms=llik_emulator_list, N_terms=length(llik_emulator_list))
    
    # Ensure that the `input_names` are the same for all llik terms. They are allowed 
    # to be ordered differently for each term, but viewed as sets they must be equal. 
    equal_set <- function(x,y) if(setequal(x,y)) x else FALSE
    input_names_list <- lapply(llik_emulator_list, function(obj) obj$input_names)
    assert_that(!isFALSE(Reduce(equal_set, input_names_list)), msg="Different `input_names` found for different llik terms.")
  
    callSuper(llik_label=term_lbls, lik_description=lik_description, emulator_description=emulator_description,
              emulator_model=NULL, default_conditional=default_conditional, 
              default_normalize=default_normalize, dim_input=llik_emulator_list[[1]]$dim_input,
              input_names=llik_emulator_list[[1]]$input_names,
              use_fixed_lik_par=all(sapply(llik_emulator_list, function(x) x$use_fixed_lik_par)), lik_par=NULL)
  },
  
  get_llik_term_attr = function(attr_name, labels=llik_label) {
    
    attr_list <- list()
    for(i in seq_along(labels)) {
      lbl <- labels[i]
      attr_list[[lbl]] <- llik_emulator_terms[[lbl]]$field(attr_name)
    }
    
    return(attr_list)
    
  },
  
  # TODO: fix creating of lists with names in other methods. 
  get_lik_par = function(lik_par_val=NULL, labels=llik_label) {
    lik_par_list <- list()
    for(i in seq_along(labels)) {
      lbl <- labels[i]
      lik_par_list[[lbl]] <- llik_emulator_list[[lbl]]$get_lik_par(lik_par_val[[lbl]])
    }
    
    return(lik_par_list)
  },
  
  sample_emulator = function(input, N_samp=1, labels=llik_label, ...) {
    
    samp_list <- list()
    for(i in seq_along(labels)) {
      lbl <- labels[i]
      samp_list[[lbl]] <- llik_emulator_terms[[lbl]]$sample_emulator(input, N_samp=N_samp, ...)
    }   
  
    return(samp_list)
  },
  
  assemble_llik = function(emulator_vals_list, lik_par=NULL, conditional=default_conditional, 
                           normalize=default_normalize, sum_terms=TRUE, labels=names(emulator_vals_list), ...) {
    
    llik_list <- list()
    for(i in seq_along(labels)) {
      lbl <- labels[i]
      llik_list[[lbl]] <- llik_emulator_terms[[lbl]]$assemble_llik(emulator_vals_list[[lbl]], lik_par[[lbl]], 
                                                                   conditional=conditional, normalize=normalize)
    }
    
    # Combine into single array. Term label names assigned to third dimension. 
    llik_vals <- abind(llik_list, along=3)
    
    if(sum_terms) llik_vals <- rowSums(llik_vals, dims=2)
    return(llik_vals)
  },
  
  sample = function(input, lik_par=NULL, N_samp=1, conditional=default_conditional, 
                    normalize=default_normalize, sum_terms=TRUE, labels=llik_label, ...) {
  
    emulator_samp_list <- sample_emulator(input, N_samp=N_samp, labels=term_labels, ...)
    asemble_llik(emulator_samp_list, lik_par, conditional=conditional, normalize=normalize,
                 sum_terms=sum_terms, labels=labels)
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
  
  initialize = function(gp_model, llik_lbl, N_obs, sig2=NULL, default_conditional=FALSE, 
                        default_normalize=FALSE, use_fixed_lik_par=FALSE, ...) {
    assert_that(inherits(gp_model, "gpWrapper"), msg="`gp_model` must inherit from `gpWrapper` class.")
    assert_that(is.integer(N_obs) && (length(N_obs) == 1) && (N_obs > 0),
                msg="`N_obs` must be an integer greater than 0.")
    assert_that(gp_model$Y_dim==1, msg="`llikEmulatorMultGausGP` only supports single-output GP emulator.")
    assert_that(!is.null(gp_model$X_names) && noNA(gp_model$X_names), 
                msg="`llikEmulatorMultGausGP` requires that `gp_model` has `X_names` field set.")
    if(!is.null(sig2)) assert_that(is.numeric(sig2) && (sig2>0), msg="`sig2` must be NULL or numeric positive value.")
    
    initFields(N_obs=N_obs)
    callSuper(emulator_model=gp_model, llik_label=llik_lbl, lik_par=sig2, input_names=gp_model$X_names,
              dim_input=gp_model$X_dim, default_conditional=default_conditional, 
              default_normalize=default_normalize, use_fixed_lik_par=use_fixed_lik_par, 
              lik_description="Multiplicative Gaussian.",
              emulator_description="GP emulating sum of squared error function.", ...)
  }, 
  
  assemble_llik = function(SSR, lik_par=NULL, conditional=default_conditional, normalize=default_normalize) {
    # SSR should be N_input x N_samp. 
    
    # Fetch the variance parameter. 
    sig2 <- get_lik_par(lik_par)
    
    # Construct likelihood using SSR.  
    llik <- -0.5 * SSR / sig2
    if(normalize || !conditional) llik <- llik - 0.5*N_obs*log(sig2)
    if(normalize) llik <- llik - 0.5*N_obs*log(2*pi)
    
    return(llik)
    
  },
  
  sample_emulator = function(input, N_samp=1, use_cov=FALSE, include_nugget=TRUE, ...) {
                             
    # Sample SSR. 
    samp <- emulator_model$sample(get_input(input), use_cov=use_cov, include_nugget=include_nugget, 
                                  N_samp=N_samp)[,,1]
    
    # Set negative samples to 0. 
    if(any(samp < 0)) {
      message("Warning: Setting negative GP quadratic error samples to 0.")
      samp[samp < 0] <- 0
    }
    
    return(samp)
  },
  
  sample = function(input, lik_par=NULL, N_samp=1, use_cov=FALSE, include_nugget=TRUE,
                    conditional=default_conditional, normalize=default_normalize, ...) {
    # Sample SSR. 
    samp <- sample_emulator(get_input(input), N_samp, use_cov, include_nugget)
    
    # Compute unnormalized or normalized log-likelihood. 
    assemble_llik(samp, lik_par, conditional, normalize)
  }

)

















# -----------------------------------------------------------------------------
# llikSumEmulatorMultGausGP: Encapsulates an approximation of a multiplicative 
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

