# numerical_experiment_functions.r
#
# Functions related to numerical tests for emulator-assisted model calibration. This 
# includes functions related to specific forward models (e.g. VSEM) and well as 
# helper functions for generating priors and plotting data. 


get_IG_priors_numerical_test <- function(sig2_true, bias_frac = 0, coef_var = 0.1, return_prior_plots = FALSE, 
                                         output_variables = NULL, bins = 50) {
  # A convenience function to generate inverse gamma priors on the variance parameters in the product Gaussian 
  # likelihood. The priors are constructed by setting the mean of the inverse gamma priors to 
  # sig2_true * (1 + bias_frac), so that a `bias_frac` of 0 implies the prior will be centered on the true value. 
  # The spread of the prior is determined by passing the coefficient of variation `coef_var`. Together, the
  # mean and coefficient of variation determine the shape and scale parameters of the inverse gamma distributions. 
  # Optionally, this function also returns histograms of samples from the prior distributions. 
  #
  # Args:
  #    sig2_true: numeric(p), vector of length p containing the true variance parameters for each of the p outputs. 
  #    bias_frac: numeric, described above. Can be a vector of length p, or length 1, in which case the same 
  #               value will be used for each of the p outputs. 
  #    coef_var: numeric, the coefficient of variation. Can be a vector of length p, or length 1, in which case the same 
  #              value will be used for each of the p outputs.
  #    return_prior_plots: logical(1), if TRUE returns list of prior histogram plots in addition to the prior specifications. 
  #    output_variables: character(p), vector of output variable names. Used to label the plots. 
  #    bins: integer(1), number of bins in the histogram plots. 
  #
  # Returns:
  #    If `return_prior_plots` is FALSE, returns list containing the prior information. Otherwise, returns a list 
  #    where the prior info list is the first element, and the second element is a list of prior histogram plots. 
  
  # Parameters for inverse gamma priors. 
  prior_means <- (1 + bias_frac) * sig2_true
  IG_shape_params <- 2 + (1 / coef_var^2)
  IG_scale_params <- sig2_true * (IG_shape_params - 1)
  
  if(any(IG_shape_params < 2)) stop("Inverse Gamma shape parameter(s) less than 2, meaning variance is infinite.")
  
  sig_eps_prior_params <- list(dist = "IG", 
                               IG_shape = IG_shape_params, 
                               IG_scale = IG_scale_params)
  
  if(!return_prior_plots) return(sig_eps_prior_params)
  
  # Sample from prior and produce histograms. 
  N_samp <- 10000
  p <- length(sig2_true)
  prior_sig_eps_samples <- matrix(NA, nrow = N_samp, ncol = p)
  if(is.null(output_variables)) output_variables <- paste0("output", 1:p)
  colnames(prior_sig_eps_samples) <- output_variables
  
  for(i in seq_len(N_samp)) {
    prior_sig_eps_samples[i,] <- sample_prior_Sig_eps(sig_eps_prior_params)
  }
  
  plts <- vector(mode = "list", length = p)
  for(j in seq_len(p)) {
    x_col <- output_variables[j]
    plts[[j]] <- ggplot(as.data.frame(prior_sig_eps_samples), aes(x = .data[[x_col]])) + 
                  geom_histogram(bins = bins) + 
                  ggtitle(paste0("Prior Samples: sig2, ", x_col)) + 
                  ylab(x_col) + 
                  geom_vline(xintercept = sig2_true[j], color = "red")
  }
  
  return(list(prior = sig_eps_prior_params, plots = plts))
  
}



