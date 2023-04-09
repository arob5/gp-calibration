#
# gp_emulator_functions.r
#


# ------------------------------------------------------------------------------
# Data Pre and Post-processing
# ------------------------------------------------------------------------------

prep_GP_training_data <- function(X = NULL, Y = NULL, scale_X = FALSE, normalize_Y = FALSE) {
  # Preps training data inputs X and outputs Y for fitting Gaussian Process (GP) model by 
  # scaling each input variable to the unit interval, and normalizing the output variable 
  # by subtracting its mean and dividing by its standard deviation. The bounds used to 
  # standardize X are computed using the range of each input variable in the training set X.
  # If Y has multiple columns, each column is treated as a different output and each column
  # is normalized independently. 
  #
  # Source for matrix operations code to normalize X: find_reps() function of hetGP package. 
  #
  # Args:
  #    X: matrix of dimension N x d, where N is the number of training points and d is the dimension 
  #       of the input space. 
  #    Y: matrix of dimension N x p where p is the number of output variables.
  #    scale_X: logical, if TRUE, maps X to d-dimensional unit hypercube. 
  #    normalize_Y: logical, if TRUE, transforms y via Z-score. 
  #
  # Returns:
  #    list, with named elements "X", "Y", "input_bounds", and "output_stats". X and Y are the 
  #    training inputs and outputs, and may or may not be standardized/normalized depending on 
  #    the arguments `scale_X` and `normalize_Y`. `input_bounds` is a 2 x d matrix summarizing the 
  #    range of each input variable used to standardize X. The first row contains the minimum 
  #    value of each of the respective input variables in the training set, and similarly the 
  #    second row stores the maxima. `output_stats` is a named vector with names 
  #    "mean_Y" and "var_Y" storing the mean and variance of Y used to compute the Z-scores. 
  #    If `scale_X` is FALSE then "input_bounds" will be NULL and likewise with 
  #    `normalize_Y` and "output_stats".
  
  if(!is.null(X) && scale_X) {
    input_bounds <- apply(X, 2, range)
    X <- scale_input_data(X, input_bounds)
  } else {
    input_bounds <- NULL
  }
  
  if(normalize_Y) {
    output_stats <- rbind(apply(Y, 2, mean), apply(Y, 2, var))
    rownames(output_stats) <- c("mean_Y", "var_Y")
    Y <- normalize_output_data(Y, output_stats)
  } else {
    output_stats <- NULL
  }
  
  return(list(X = X, Y = Y, input_bounds = input_bounds, output_stats = output_stats))
  
}


scale_input_data <- function(X, input_bounds, inverse = FALSE) {
  # Transforms input data X (typically points at which to predict) by performing a linear scaling encoded
  # in `input_bounds`. Can also perform the inverse transformation. 
  #
  # Args:
  #    X: matrix of dimension N_pred x d, where N_pred is the number of input points and d is the dimension of the 
  #       input space. 
  #    input_bounds: 2 x d matrix summarizing the range of each input variable used to standardize the training inputs. 
  #                  The first row contains the minimum value of each of the respective input variables in the training set, 
  #                  and similarly the second row stores the maxima. This object is returned by prep_GP_training_data(). 
  #    inverse: logical, if TRUE, treats X as already scaled and reverses the scaling. Otherwise performs the forward 
  #             transformation. Default is FALSE. 
  #
  # Returns:
  #    matrix of dimension N_pred x d; the matrix X whose columns have been scaled according to `input_bounds`. 
  #
  
  if(inverse) {
    X <- X %*% diag(input_bounds[2,] - input_bounds[1,], ncol(X)) + matrix(input_bounds[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  } else {
    X <- (X - matrix(input_bounds[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(input_bounds[2,] - input_bounds[1,]), ncol(X))
  }
  
  return(X)
}


normalize_output_data <- function(Y, output_stats, inverse = FALSE) {
  # Transforms output data by transforming Y to Z-scores, or performing the reverse operation. 
  # Can handle multiple output variables, where each column of `Y` is treated as an output 
  # variable and is normalized independently. 
  #
  # Args:
  #    Y: matrix of dimension N x p where p is the number of output variables. 
  #    output_stats: a matrix of dimensions 2xp, where p is the number of output variables. The matrix 
  #                  must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of each output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). 
  #
  # Returns: 
  #    matrix of dimension N x p, the transformed version of Y. 
  
  if(inverse) {
    Y <- Y * matrix(sqrt(output_stats["var_Y",]), nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE) + 
      matrix(output_stats["mean_Y",], nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  } else {
    Y <- (Y - matrix(output_stats["mean_Y",], nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)) / 
      matrix(sqrt(output_stats["var_Y",]), nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  }
  
  return(Y)
  
}



# ------------------------------------------------------------------------------
# Fitting GPs
# ------------------------------------------------------------------------------

fit_GP <- function(X_train, y_train, gp_lib, gp_kernel) {
  # Estimate kernel and mean function hyperparameters for a univariate GP regression 
  # using specified Gaussian Process (GP) package. 
  #
  # Args:
  #    X_train: matrix of shape N x d, where N is the number of design (training) points
  #             and d is the dimension of the input space. 
  #    y_train: matrix of shape N x 1, the training outputs corresponding to the training inputs. 
  #    gp_lib: character, string specifying the GP package to use. Currently supports 
  #            "mlegp" and "hetGP". 
  #    gp_kernel: character, string specifying the kernel family to use. Potential options are 
  #               limited by the specified GP library. Currently only supports "Gaussian" (i.e. the 
  #               squared exponential/exponentiated quadratic/radial basis function kernel). 
  #
  # Returns:
  #    Returns, list with named elements "fit" and "time". The former stores the fit GP object 
  #    returned by the fitting function in the specified GP package. The latter contains the 
  #    time taken by the fitting procedure. 
  
  eps_nug <- sqrt(.Machine$double.eps)
  
  tic <- proc.time()[3]
  if(gp_lib == "mlegp") {
    gp_fit <- mlegp(X_train, y_train, nugget.known = 1, nugget = eps_nug, constantMean = 1)
  } else if(gp_lib == "hetGP") {
    gp_fit <- mleHomGP(X_train, y_train, covtype = gp_kernel, known = list(g = eps_nug))
  } else {
    stop("Invalid GP library: ", gp_lib)
  }
  
  toc <- proc.time()[3]
  gp_fit_time <- toc - tic
  
  return(list(fit = gp_fit, time = gp_fit_time))
  
}

 
# TODO: Parallelize independent GP fitting
fit_independent_GPs <- function(X_train, Y_train, gp_lib, gp_kernel) {
  # Builds on top of `fit_GP()` to generalize to multivariate GP regression. Simply 
  # fits independent GPs for each output (where the same input points are used for each GP).
  #
  # Args:
  #    X_train: matrix of shape N x d, where N is the number of design (training) points
  #             and d is the dimension of the input space. 
  #    y_train: matrix of shape N x p, with jth column containing the training outputs 
  #             for the jth output variable corresponding to the training inputs. 
  #    gp_lib: character, string specifying the GP package to use. Currently supports 
  #            "mlegp" and "hetGP". 
  #    gp_kernel: character, string specifying the kernel family to use. Potential options are 
  #               limited by the specified GP library. Currently only supports "Gaussian" (i.e. the 
  #               squared exponential/exponentiated quadratic/radial basis function kernel).
  #
  # Returns:
  #    list, with named elements "fits" and "times". The first is itself a list, 
  #    with length `ncol(Y_train)`; i.e. length equal to the number of outputs. The jth 
  #    element stores the fit GP object returned by `fit_GP()`, containing the fit GP for the jth output. 
  #    The "times" element is a numeric vector storing the times required to fit each GP. 
  
  p <- ncol(Y_train)
  GP_objects <- vector(mode = "list", length = p)
  GP_fit_times <- vector(mode = "numeric", length = p)
  
  for(j in seq_len(p)) {
    GP_fit_list <- fit_GP(X_train, Y_train[,j, drop = FALSE], gp_lib, gp_kernel)
    GP_objects[[j]] <- GP_fit_list[["fit"]]
    GP_fit_times[j] <- GP_fit_list[["time"]]
  }
  
  return(list(fits = GP_objects, times = GP_fit_times))
  
}


# ------------------------------------------------------------------------------
# Predicting with GPs
# ------------------------------------------------------------------------------

# TODO: 
#    - Currently sd2_nug should not be trusted; e.g. it is not corrected in the case of truncation/rectification. Might just be easier to 
#      drop this or combine it with the pointwise variance predictions. 
predict_GP <- function(X_pred, gp_obj, gp_lib, cov_mat = FALSE, denormalize_predictions = FALSE,
                       output_stats = NULL, exponentiate_predictions = FALSE, transformation_method = NA_character_) {
  # Calculate GP predictive mean, variance, and optionally covariance matrix at specified set of 
  # input points. 
  #
  # Args:
  #    X_pred: matrix of dimension N_pred x d, where d is the dimension of the input space. Each row is an input 
  #            point at which to predict. 
  #    gp_obj: An object representing a GP fit, which will differ based on the library used to fit the GP. 
  #    gp_lib: character(1), the library used to fit the GP. Currently supports "mlegp" or "hetGP". 
  #    cov_mat: logical(1), if TRUE, calculates and returns the N_pred x N_pred predictive covariance matrix 
  #             over the set of input points. Otherwise, only calculates the pointwise predictive variances. 
  #    denormalize_predictions: logical(1), if TRUE, applies linear transformation to predictions, 
  #                             inverting the Z-score transformation.
  #    output_stats: If not NULL, then a matrix of dimensions 2x1. The matrix 
  #                  must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of the output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). Only required if `transform_predictions` is TRUE. 
  #    exponentiate_predictions: logical(1), if TRUE converts Gaussian predictive mean/variance to log-normal 
  #                              predictive mean variance. This is useful if the GP model was fit to a log-transformed
  #                              response variable, but one wants predictions on the original scale. Default is FALSE. 
  #    sampling_method: character(1), if "default" treats the predictive GP distribution as Gaussian. If "rectified", treats
  #                     it as rectified Gaussian and if "truncated" treats it as truncated Gaussian (truncated at 0).
  #
  # Returns:
  #    list, with named elements "mean", "sd2", "sd2_nug", "cov". 
  
  pred_list_names <- c("mean", "sd2", "sd2_nug", "cov")
  pred_list <- vector(mode = "list", length = length(pred_list_names))
  names(pred_list) <- pred_list_names
  
  # Generate GP Predictions. 
  if(gp_lib == "mlegp") {
    mlegp_pred <- predict(gp_obj, newData = X_pred, se.fit = TRUE)
    pred_list[["mean"]] <- mlegp_pred[["fit"]]
    pred_list[["sd2"]] <- mlegp_pred[["se.fit"]]^2
  } else if(gp_lib == "hetGP") {
    # Second matrix for computing predictive covariance
    if(cov_mat) {
      X_prime <- X_pred
    } else {
      X_prime <- NULL
    }
    hetGP_pred <- predict(gp_obj, X_pred, xprime = X_prime)
    pred_list[pred_list_names] <- hetGP_pred[c("mean", "sd2", "nugs", "cov")]
  }
  
  
  # Make adjustments for truncated or rectified Gaussian predictive distribution. 
  if(!exponentiate_predictions && ((sampling_method == "truncated") || (sampling_method == "rectified"))) {
    n <- length(pred_list[["mean"]])
    mu_trunc <- sapply(seq(1, n), function(i) etruncnorm(a=0, b=Inf, mean = pred_list[["mean"]][i], sd = sqrt(pred_list[["sd2"]][i])))
    var_trunc <- sapply(seq(1, n), function(i) vtruncnorm(a=0, b=Inf, mean = pred_list[["mean"]][i], sd = sqrt(pred_list[["sd2"]][i])))
    
    if(sampling_method == "truncated") {
      pred_list[["mean"]] <- mu_trunc
      pred_list[["sd2"]] <- var_trunc
    } else if(sampling_method == "rectified") {
      rect_Gaussian_moments <- transform_GP_to_rectified_GP(pred_list[["mean"]], pred_list[["sd2"]], mu_trunc = mu_trunc, sig2_trunc = var_trunc)
      
      pred_list[["mean"]] <- rect_Gaussian_moments[["mean"]]
      pred_list[["sd2"]] <- rect_Gaussian_moments[["var"]]
    }
    
  }
  
  # Invert Z-score transformation of response variable. 
  if(denormalize_predictions) {
    pred_list[["mean"]] <- output_stats["mean_Y",1] + sqrt(output_stats["var_Y",1]) * pred_list[["mean"]]
    pred_list[["sd2"]] <- output_stats["var_Y",1] * pred_list[["sd2"]]
    pred_list[["sd2_nug"]] <- output_stats["var_Y",1] * pred_list[["sd2_nug"]]
    if(cov_mat) {
      pred_list[["cov"]] <- output_stats["var_Y",1] * pred_list[["cov"]]
    }
  }
  
  
  # Apply transformation to GP predictions. 
  if(!is.na(transformation_method)) {
    
  }
  
  
  
  
  # Transform log-transformed predictions back to original scale (for log-normal process predictions). 
  if(exponentiate_predictions) {
    pred_list_exp <- transform_GP_to_LNP(pred_list$mean, pred_list$sd2, pred_list$cov)
    pred_list_exp[["sd2_nug"]] <- transform_GP_to_LNP(gp_mean = pred_list$mean, gp_var = pred_list$sd2_nug)$var
    pred_list <- pred_list_exp
  } 
  
  return(pred_list)
  
}


predict_independent_GPs <- function(X_pred, gp_obj_list, gp_lib, cov_mat = FALSE, 
                                    denormalize_predictions = FALSE, output_stats = NULL, 
                                    exponentiate_predictions = FALSE, sampling_method = "default") {
  # A wrapper function for predict_GP() that generalizes the latter to generating predictions for 
  # multi-output GP regression using independent GPs. 
  #
  # Args:
  #    X_pred: matrix of dimension N_pred x d, where d is the dimension of the input space. Each row is an input 
  #            point at which to predict. 
  #    gp_obj_list: A list of GP objects, each of which represents a GP fit to one of the outputs. The objects will 
  #                differ based on the specific GP library used for fitting. All of the objects in the list must have 
  #                been fit using the same library. 
  #    gp_lib: character(1), the library used to fit the GP. Currently supports "mlegp" or "hetGP". 
  #    cov_mat: logical(1), if TRUE, calculates and returns the N_pred x N_pred predictive covariance matrix 
  #             over the set of input points. Otherwise, only calculates the pointwise predictive variances. 
  #    denormalize_predictions: logical(1), if TRUE, applies linear transformation to predictions, 
  #                             inverting the Z-score transformation.
  #    output_stats: If not NULL, then a matrix of dimensions 2xp, where p is the number of output variables. The matrix 
  #                  must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of each output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). Only required if `transform_predictions` is TRUE. 
  #    exponentiate_predictions: logical(1), if TRUE converts Gaussian predictive mean/variance to log-normal 
  #                              predictive mean variance. This is useful if the GP model was fit to a log-transformed
  #                              response variable, but one wants predictions on the original scale. Default is FALSE. 
  #    sampling_method: character(1), if "default" treats the predictive GP distribution as Gaussian. If "rectified", treats
  #                     it as rectified Gaussian and if "truncated" treats it as truncated Gaussian (truncated at 0).
  # 
  # Returns:
  #    list, with length equal to the length of `gp_obj_list`. Each element of this list is itself a list, 
  #    with named elements "mean", "sd2", "sd2_nug", "cov" (the output of the function `predict_GP()` applied 
  #    to each GP in `gp_obj_list`). 
  
  lapply(seq_along(gp_obj_list), function(j) predict_GP(X_pred, gp_obj_list[[j]], gp_lib, cov_mat, denormalize_predictions, 
                                                        output_stats[,j,drop=FALSE], exponentiate_predictions, sampling_method = sampling_method))
  
}


# ------------------------------------------------------------------------------
# Transforming GPs
# ------------------------------------------------------------------------------

transform_GP_predictions <- function(gp_mean, gp_var, transformation_method, ...) {
  
  if(transformation_method == "LNP") {
    return(transform_GP_to_LNP(gp_mean, gp_var))
  } else if(tranformation_method == "truncated") {
    return(transform_GP_to_truncated_GP(gp_mean, gp_var))
  } else if(transformation_method == "rectified") {
    return(transform_GP_to_rectified_GP(gp_mean, gp_var))
  }
  
}


transform_GP_to_truncated_GP <- function(gp_mean, gp_var) {
  # Compute the mean and variance of the distribution of the conditional random variable X|X > 0, 
  # where X ~ N(mu, sig2). This is vectorized so that the arguments `mu` and `sig2` can be vectors 
  # where the ith element pertains to a random variable X_i ~ N(mu_i, sig2_i). 
  #
  # Args:
  #    gp_mean: numeric(), vector of means of the Gaussian variables X_i. 
  #    gp_var: numeric(), vector of equal length as `mu` of variances of the Gaussian variables X_i. 
  #
  # Returns:
  #    list, with elements "mean", "var". The elements contain vectors storing the mean and 
  #    variance, respectively, of the truncated Gaussian distribution.
  
  n <- length(gp_mean)
  gp_mean_trunc<- sapply(seq(1, n), function(i) etruncnorm(a=0, b=Inf, mean = gp_mean[i], sd = sqrt(gp_var[i])))
  gp_var_trunc <- sapply(seq(1, n), function(i) vtruncnorm(a=0, b=Inf, mean = gp_mean[i], sd = sqrt(gp_var[i])))
  
  return(list(mean = gp_mean_trunc, 
              var = gp_var_trunc))
  
}


compute_zero_truncated_Gaussian_moments <- function(gp_mean, gp_var) {
  # Compute the mean and variance of the distribution of the conditional random variable X|X > 0, 
  # where X ~ N(mu, sig2). This is vectorized so that the arguments `mu` and `sig2` can be vectors 
  # where the ith element pertains to a random variable X_i ~ N(mu_i, sig2_i). This function plays 
  # the same exact role as `transform_GP_to_truncated_GP()` but implements the transformations 
  # manually, rather than relying on the "truncnorm" package. 
  #
  # Args:
  #    gp_mean: numeric(), vector of means of the Gaussian variables X_i. 
  #    gp_var: numeric(), vector of equal length as `mu` of variances of the Gaussian variables X_i. 
  #
  # Returns:
  #    list, with elements "mean", "var", and "Z". The first two contain vectors storing the mean and 
  #    variance, respectively, of the truncated Gaussian distribution. The third is a vector of 
  #    evaluations of the form P(X_i > 0). This is returned primarily as it can be used to compute 
  #    the moments of the rectified Gaussian. 
  
  sig <- sqrt(gp_var)
  alpha <- -gp_mean / sig
  phi_alpha <- dnorm(alpha)
  Z <- 1 - pnorm(alpha)
  
  mu_trunc <- gp_mean + sig * phi_alpha / Z
  sig2_trunc <- sig2 * (1 + (alpha * phi_alpha) / Z - (phi_alpha / Z)^2)
  
  return(list(mean = mu_trunc, var = sig2_trunc, Z = Z))
  
}


transform_GP_to_rectified_GP <- function(gp_mean, gp_var, gp_mean_trunc = NULL, gp_var_trunc = NULL, Z = NULL) {
  # Compute the mean and variance of the distribution of the random max(X, 0), 
  # where X ~ N(mu, sig2). This is vectorized so that the arguments `mu` and `sig2` can be vectors 
  # where the ith element pertains to a random variable X_i ~ N(mu_i, sig2_i). The rectified Gaussian 
  # moments are closely related to the truncated Gaussian moments, so optionally computations from the 
  # latter can be passed in to avoid additional computation. In this case, the arguments 
  # `mu_trunc`, `sig2_trunc` must both be passed. The argument `Z` can additionally be passed, or it can 
  # be computed separately. 
  #
  # Args:
  #    gp_mean: numeric(), vector of means of the Gaussian variables X_i. 
  #    gp_var: numeric(), vector of equal length as `gp_mean` of variances of the Gaussian variables X_i. 
  #    gp_mean_trunc: numeric(), vector of equal length as `gp_mean` of the means of the truncated Gaussian 
  #              variables X_i|X_i > 0. If NULL, will not use partial computations from the truncated 
  #              Gaussian computation. 
  #    gp_var_trunc: numeric(), same as `gp_mean_trunc` but containing the truncated Gaussian variances. 
  #    Z: numeric(), same as `gp_mean_trunc` but containing vector of evaluations P(X_i > 0). 
  #
  # Returns:
  #    list, with elements "mean" and "var"; the vectors storing the mean and 
  #    variance, respectively, of the rectified Gaussian distribution.
  
  if(is.null(gp_mean_trunc) || is.null(gp_var_trunc)) {
    gp_trunc_moments <- transform_GP_to_truncated_GP(gp_mean, gp_var)
    gp_mean_trunc <- gp_trunc_moments$mean
    gp_var_trunc <- gp_trunc_moments$var
    Z <- gp_trunc_moments$Z
  }
  
  if(is.null(Z)) {
    Z <- 1 - pnorm(-gp_mean / sqrt(gp_var))
  }
  
  gp_mean_rect = gp_mean_trunc * Z
  return(list(mean = gp_mean_rect, 
              var = (gp_var_trunc + gp_mean_trunc^2) * Z - gp_mean_rect^2))
}


transform_GP_to_LNP <- function(gp_mean = NULL, gp_var = NULL, gp_cov = NULL) {
  # Transforms distribution y ~ GP to exp(y) ~ LNP. The optional arguments are predictive means, 
  # variances, and covariance matrix from a Gaussian Process. This function will transform 
  # these into the means, variances, and covariance matrix of the log-normal process obtained
  # by exponentiating the Gaussian Process. The mean and either the variances or covariance matrix 
  # must be provided, since the log-normal mean and variance depend on both the GP mean and variance. 
  #
  # Args:
  #    gp_mean: numeric, vector of GP mean predictions. 
  #    gp_var: numeric, vector of GP variance predictions. Must be ordered to correspond to `gp_mean`.
  #    gp_cov: matrix, GP predicted covariance matrix. Must be ordered to correspond to `gp_mean`.
  #                  
  # Returns:
  #    list, containing the predictions obeying the distribution of the exponentiated GP. 
  #    The list will contain named elements "mean", "var", and "cov". The 
  #     element "mean" will be non-NULL, containing the vector transformed means. 
  #    If `gp_cov` is non-NULL then the element "cov" will also be non-NULL. If
  #    `gp_var` is non-NULL then the element "var" will also be non-NULL. 
  
  include_var <- TRUE
  if(is.null(gp_var)) {
    include_var <- FALSE
    if(is.null(gp_pred_cov)) {
      stop("Either gp_mean or gp_var must be provided.")
    }
    gp_var <- diag(gp_cov)
  }
  
  output_list <- list(mean = NULL, var = NULL, cov = NULL)
  output_list[["mean"]] <- exp(gp_mean + 0.5 * gp_var)
  
  if(!is.null(gp_cov)) {
    N_obs <- length(gp_mean)
    mu_mat <- matrix(gp_mean, nrow = N_obs, ncol = N_obs, byrow = TRUE)
    sd2_mat <- matrix(gp_var, nrow = N_obs, ncol = N_obs, byrow = TRUE)
    output_list[["cov"]] <- exp(mu_mat + t(mu_mat) + 0.5 * (sd2_mat + t(sd2_mat))) * (exp(gp_cov) - 1)
  }
  
  if(include_var) {
    if(!is.null(gp_cov)) {
      output_list[["var"]] <- diag(output_list[["cov"]])
    } else {
      output_list[["var"]] <- (exp(gp_var) - 1) * exp(2*gp_mean + gp_var)
    }
  }
  
  return(output_list)
  
}


# ------------------------------------------------------------------------------
# Plotting GPs
# ------------------------------------------------------------------------------

plot_gp_fit_1d <- function(X_test, y_test, X_train, y_train, gp_mean_pred, gp_var_pred, 
                           exponentiate_predictions = FALSE, log_scale = FALSE, vertical_line = NULL,
                           xlab = "", ylab = "", main_title = "", CI_prob = 0.9, transformation = "default") {
  # Core function for producing plots for GP predictions with one-dimensional input space. The function plots
  # the true, known latent function values at the design inputs and test inputs. It also plots the 
  # GP predictive mean and confidence bands at the test inputs. This function assunmes that `gp_mean_pred` and 
  # `gp_var_pred` are the means and variances of a Gaussian process, but allows these predictions to be transformed
  # by 1.) exponentiating the predictions, resulting in a log-normal process, 2.) converting the Gaussian predictive 
  # distribution to a left-truncated (at zero) Gaussian distribution, or 3.) converting the Gaussian predictive 
  # distribution to a rectified Gaussian distribution. The plot can also optionally include
  # a vertical line corresponding to some "true" parameter in the input space. By setting `log_scale` to TRUE, the y-axis
  # will be set to a log base 10 sdale. 
  # 
  # Args:
  #    X_test: matrix, of dimension M x 1 where M is the number of test input points. 
  #    y_test: numeric(M), the vector of true outputs at the M test points. 
  #    X_train: matrix, of dimension N x 1, where N is the number of design/training points. 
  #    y_train: numeric(N), the vector of true outputs at the N design points (the training response values). 
  #    gp_mean_pred: numeric(M), vector of GP predictive mean at the test input points. 
  #    gp_var_pred: numeric(M), vector of GP predictive variance at the test input points. 
  #    exponentiate_predictions: logical(1), if TRUE, produces a log-normal process plot by exponentiating the 
  #                              GP. Default is FALSE. Note that if `exponentiate_predictions` is TRUE then the 
  #                              assumption is that the GP predictions must be exponentiated, but the baseline 
  #                              data (X_test, y_test, X_train, y_train) is on the right scale and is not modified 
  #                              fpr the plot. 
  #    log_scale: logical(1), if TRUE the y-axis of the plot will be on a log base 10 scale. In particular, the 
  #               plotted y values will be transformed as log10(y) and the y-axis labels will be printed as 
  #               10^1, 10^2, etc. This has nothing to do with the log-normal process described above; the plot can 
  #               be displayed on a log-scale whether or not `exponentiate_predictions` is TRUE. To avoid taking the 
  #               log of non-positive values, the actual transformation performed is log10(y + C), for a constant 
  #               C determined by the data to plot. 
  #    vertical_line: numeric(1), if non-NULL this is the x-intercept at which a vertical dashed line will be printed.
  #                   This typically represents some true/baseline parameter in the input space. If NULL, no vertical 
  #                   line will be included. 
  #    xlab: character(1), the label/title for the x-axis. 
  #    ylab: character(1), the label/title for the y-axis.
  #    main_title: character(1), the main title for the plot. 
  #    CI_prob: numeric(1), value in (0, 1) determining the confidence intervals that will be plotted. e.g. 0.9 means
  #             90% confidence intervals (the default). 
  #    transformation: character(1), if "truncated", will convert the GP predictive distribution to truncated Gaussian distribution. 
  #                    If "rectified", will instead transform to rectified Gaussian. If "default", will not transform the predictive 
  #                    distribution. 
  #
  # Returns:
  #    A ggplot2 plot object. This can be used to display or to further modify the plot. 
  
  order_pred <- order(X_test)
  order_train <- order(X_train)
  gp_sd_pred <- sqrt(gp_var_pred)
  n <- length(gp_mean_pred)
  
  # Confidence intervals
  CI_tail_prob <- 0.5 * (1 - CI_prob)
  CI_plot_label <- paste0(CI_prob * 100, "% CI")
  if(exponentiate_predictions) {
    CI_lower <- qlnorm(CI_tail_prob, gp_mean_pred, sqrt(gp_var_pred), lower.tail = TRUE)
    CI_upper <- qlnorm(CI_tail_prob, gp_mean_pred, sqrt(gp_var_pred), lower.tail = FALSE)
    transformed_predictions <- transform_GP_to_LNP(gp_mean = gp_mean_pred, gp_var = gp_var_pred)
    gp_mean_pred <- transformed_predictions$mean
    gp_var_pred <- transformed_predictions$sd2
  } else if(transformation == "truncated") {
    CI_lower <- qtruncnorm(CI_tail_prob, a = 0, b = Inf, mean = gp_mean_pred, sd = sqrt(gp_var_pred))
    CI_upper <- qtruncnorm(1 - CI_tail_prob, a = 0, b = Inf, mean = gp_mean_pred, sd = sqrt(gp_var_pred))
    gp_mean_pred_trunc <- sapply(seq(1, n), function(i) etruncnorm(a=0, b=Inf, mean = gp_mean_pred[i], sd = sqrt(gp_var_pred[i])))
    gp_var_pred <- sapply(seq(1, n), function(i) vtruncnorm(a=0, b=Inf, mean = gp_mean_pred[i], sd = sqrt(gp_var_pred[i])))
    gp_mean_pred <- gp_mean_pred_trunc
  } else if(transformation == "rectified") {
    # TODO: CIs for rectified Gaussian
    transformed_predictions <- transform_GP_to_rectified_GP(gp_mean_pred, gp_var_pred)
    gp_mean_pred <- transformed_predictions$mean
    gp_var_pred <- transformed_predictions$var
  } else {
    CI_lower <- qnorm(CI_tail_prob, gp_mean_pred, sqrt(gp_var_pred), lower.tail = TRUE)
    CI_upper <- qnorm(CI_tail_prob, gp_mean_pred, sqrt(gp_var_pred), lower.tail = FALSE)
  }
  
  if(log_scale) {
    shift <- -min(CI_lower) + 1
    y_test <- log10(y_test + shift)
    CI_lower <- log10(CI_lower + shift)
    CI_upper <- log10(CI_upper + shift)
    y_train <- log10(y_train + shift)
    gp_mean_pred <- log10(gp_mean_pred + shift)
  }
  
  df_test <- data.frame(x_test = X_test[order_pred,1], 
                        y_test = y_test[order_pred], 
                        y_test_pred = gp_mean_pred[order_pred],
                        CI_lower = CI_lower[order_pred], 
                        CI_upper = CI_upper[order_pred])
  df_train <- data.frame(x_train = X_train[order_train,1], 
                         y_train = y_train[order_train])
  
  gp_plot <- ggplot(data=df_test, aes(x = x_test, y = y_test_pred)) + 
    geom_line(color = "blue") + 
    geom_line(aes(y = y_test), color = "red") + 
    geom_line(aes(y = CI_lower), color = "gray") + 
    geom_line(aes(y = CI_upper), color = "gray") + 
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), fill = "gray", alpha = 0.5) + 
    geom_point(data = df_train, aes(x = x_train, y = y_train), color = "red") + 
    xlab(xlab) + 
    ylab(paste0(ylab, " (", CI_plot_label, ")")) + 
    ggtitle(main_title)
  
  # Add vertical line
  if(!is.null(vertical_line)) {
    gp_plot <- gp_plot + geom_vline(xintercept = vertical_line, linetype = 2, color = "pink1")
  }
  
  # Adjust y-axis labels if on log scale
  if(log_scale) {
    y_axis_labels <- ggplot_build(gp_plot)$layout$panel_params[[1]]$y$get_labels()
    y_axis_labels[!is.na(y_axis_labels)] <- paste0("10^", y_axis_labels[!is.na(y_axis_labels)])
    gp_plot <- gp_plot + scale_y_continuous(labels = y_axis_labels)
  }
  
  return(gp_plot)
  
}














