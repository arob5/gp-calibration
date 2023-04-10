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
  
  if(!is.null(Y) && normalize_Y) {
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
                       output_stats = NULL, transformation_method = NA_character_) {
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
  #    transformation_method: character(1), if "default" treats the predictive GP distribution as Gaussian. If "rectified", treats
  #                           it as rectified Gaussian and if "truncated" treats it as truncated Gaussian (truncated at 0). If 
  #                           "LNP", exponentiates the GP, resulting in a log-normal process. 
  #
  # Returns:
  #    list, with named elements "mean", "var", "var_nug", "cov". 
  
  pred_list_names <- c("mean", "var", "var_nug", "cov")
  pred_list <- vector(mode = "list", length = length(pred_list_names))
  names(pred_list) <- pred_list_names
  
  # Generate GP Predictions. 
  if(gp_lib == "mlegp") {
    mlegp_pred <- predict(gp_obj, newData = X_pred, se.fit = TRUE)
    pred_list[["mean"]] <- mlegp_pred[["fit"]]
    pred_list[["var"]] <- mlegp_pred[["se.fit"]]^2
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
  
  # Invert Z-score transformation of response variable. 
  if(denormalize_predictions) {
    pred_list[["mean"]] <- output_stats["mean_Y",1] + sqrt(output_stats["var_Y",1]) * pred_list[["mean"]]
    pred_list[["var"]] <- output_stats["var_Y",1] * pred_list[["var"]]
    pred_list[["var_nug"]] <- output_stats["var_Y",1] * pred_list[["var_nug"]]
    if(cov_mat) {
      pred_list[["cov"]] <- output_stats["var_Y",1] * pred_list[["cov"]]
    }
  }
  
  # Apply transformation to GP predictions. 
  if(!is.na(transformation_method)) {
    pred_list_transformed <- transform_GP_predictions(pred_list$mean, pred_list$var, transformation_method = transformation_method, gp_cov = pred_list$cov)
    pred_list_transformed[["var_nug"]] <- transform_GP_predictions(pred_list$mean, pred_list$var_nug, transformation_method = transformation_method)$var
    pred_list <- pred_list_transformed
  }
  
  return(pred_list)
  
}


# TODO: gp_lib, transformation_method, etc. should be allowed to be vectors here. 
predict_independent_GPs <- function(X_pred, gp_obj_list, gp_lib, cov_mat = FALSE, denormalize_predictions = FALSE,
                                    output_stats = NULL, transformation_method = NA_character_) {
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
  #    transformation_method: character(1), specifies the transformation to apply. The current valid options
  #                           are "LNP" (exponentiates the GP, resulting in a log-normal process), "truncated"
  #                           (transforms the GP to a left-truncated (at 0) GP) and "rectified" (transforms the 
  #                           GP to a rectified Gaussian). 
  # 
  # Returns:
  #    list, with length equal to the length of `gp_obj_list`. Each element of this list is itself a list, 
  #    with named elements "mean", "sd2", "sd2_nug", "cov" (the output of the function `predict_GP()` applied 
  #    to each GP in `gp_obj_list`). 
  
  lapply(seq_along(gp_obj_list), function(j) predict_GP(X_pred, gp_obj_list[[j]], gp_lib, cov_mat, denormalize_predictions, 
                                                        output_stats[,j,drop=FALSE], transformation_method = transformation_method))
  
}


# ------------------------------------------------------------------------------
# Transforming GPs
# ------------------------------------------------------------------------------

transform_GP_predictions <- function(gp_mean, gp_var, transformation_method, ...) {
  # A convenience function for applying a transformation to GP distribution. In particular, this 
  # function takes GP means and variances, and computes the means and variances of the transformed 
  # distribution. 
  #
  # Args:
  #    gp_mean: numeric, vector of GP mean predictions. 
  #    gp_var: numeric, vector of GP variance predictions. Must be ordered to correspond to `gp_mean`.
  #    transformation_method: character(1), specifies the transformation to apply. The current valid options
  #                           are "LNP" (exponentiates the GP, resulting in a log-normal process), "truncated"
  #                           (transforms the GP to a left-truncated (at 0) GP) and "rectified" (transforms the 
  #                           GP to a rectified Gaussian). 
  #    ...: optional arguments that can be passed along to the specific transformation functions. 
  #
  # Returns:
  #    list, typically with elements "mean" and "var" containing vectors of the means and variances of the 
  #    transformed distributions. The elements of the list can differ slightly depending on the transformation. 
  #    See the transformation-specific functions for details. 
  
  if(transformation_method == "LNP") {
    return(transform_GP_to_LNP(gp_mean, gp_var, ...))
  } else if(transformation_method == "truncated") {
    return(transform_GP_to_truncated_GP(gp_mean, gp_var))
  } else if(transformation_method == "rectified") {
    return(transform_GP_to_rectified_GP(gp_mean, gp_var, ...))
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
  sig2_trunc <- gp_var * (1 + (alpha * phi_alpha) / Z - (phi_alpha / Z)^2)
  
  return(list(mean = mu_trunc, var = sig2_trunc, Z = Z))
  
}


transform_GP_to_rectified_GP <- function(gp_mean, gp_var, gp_mean_trunc = NULL, gp_var_trunc = NULL, Z = NULL, ...) {
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


transform_GP_to_LNP <- function(gp_mean = NULL, gp_var = NULL, gp_cov = NULL, ...) {
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

plot_gp_fit_1d <- function(X_test, y_test, X_train, y_train, gp_mean_pred, gp_var_pred, log_scale = FALSE, vertical_line = NULL,
                           xlab = "", ylab = "", main_title = "", CI_prob = 0.9, transformation_method = NA_character_) {
  # Core function for producing plots for GP predictions with one-dimensional input space. The function plots
  # the true, known latent function values at the design inputs and test inputs. It also plots the 
  # GP predictive mean and confidence bands at the test inputs. This function assumes that `gp_mean_pred` and 
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
  #    log_scale: logical(1), if TRUE the y-axis of the plot will be on a log base 10 scale. In particular, the 
  #               plotted y values will be transformed as log10(y) and the y-axis labels will be printed as 
  #               10^1, 10^2, etc. This has nothing to do with the log-normal process transformation (see 
  #               `transformation_method` below); the plot can be displayed on a log-scale whether or not
  #               the log-normal process transformation is applied. To avoid taking the 
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
  #    transformation_method: character(1), if "truncated", will convert the GP predictive distribution to truncated Gaussian distribution. 
  #                           If "rectified", will instead transform to rectified Gaussian. If "LNP" will exponentiate the GP, resulting in a 
  #                           log-normal process. The default is not to perform any transformation. 
  #
  # Returns:
  #    A ggplot2 plot object. This can be used to display or to further modify the plot. 
  
  order_pred <- order(X_test)
  order_train <- order(X_train)
  gp_sd_pred <- sqrt(gp_var_pred)
  n <- length(gp_mean_pred)
  
  # Confidence intervals. 
  CI_plot_label <- paste0(CI_prob * 100, "% CI")
  confidence_intervals <- get_GP_confidence_intervals(CI_prob, gp_mean_pred, gp_var_pred, transformation_method = transformation_method)
  CI_lower <- confidence_intervals$lower
  CI_upper <- confidence_intervals$upper
  
  # Transform GP mean predictions. 
  if(!is.na(transformation_method)) {
    transformed_predictions <- transform_GP_predictions(gp_mean = gp_mean_pred, gp_var = gp_var_pred, transformation_method = transformation_method)
    gp_mean_pred <- transformed_predictions$mean
    gp_var_pred <- transformed_predictions$var
  }
  
  # Plot y-axis on log base 10 scale. 
  if(log_scale) {
    shift <- -min(CI_lower) + 1
    y_test <- log10(y_test + shift)
    CI_lower <- log10(CI_lower + shift)
    CI_upper <- log10(CI_upper + shift)
    y_train <- log10(y_train + shift)
    gp_mean_pred <- log10(gp_mean_pred + shift)
  }
  
  # Generate plot. 
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


get_GP_confidence_intervals <- function(CI_prob, gp_mean, gp_var, transformation_method = NA_character_) {
  # Computes pointwise upper and lower confidence bounds for a GP or transformed GP. Note that the mean and variance 
  # arguments always correspond to the GP, not the transformed GP. 
  #
  # Args:
  #    CI_prob: numeric(1), value in (0, 1) determining the confidence intervals that will be plotted. e.g. 0.9 means
  #             90% confidence intervals (the default). 
  #    gp_mean: numeric(), vector of means of Gaussian distributions (not the means of the transformed GPs). 
  #    gp_var: numeric(), vector of variances of the Gaussian distributions. 
  #    transformation_method: character(1), if "truncated", will convert the GP distribution to truncated Gaussian distribution. 
  #                           If "rectified", will instead transform to rectified Gaussian. If "LNP" will exponentiate the GP, resulting in a 
  #                           log-normal process. The default is not to perform any transformation. 
  #
  # Returns:
  #    list, with named elements "lower" and "upper". Each corresponds to a vector of length equal to the length of `gp_mean` and `gp_var` 
  #    containing the lower and upper 100*`CI_prob`% confidence bounds. 
  
  CI_tail_prob <- 0.5 * (1 - CI_prob)
  gp_sd <- sqrt(gp_var)
  
  if(is.na(transformation_method)) {
    CI_lower <- qnorm(CI_tail_prob, gp_mean, gp_sd, lower.tail = TRUE)
    CI_upper <- qnorm(CI_tail_prob, gp_mean, gp_sd, lower.tail = FALSE)
  } else if(transformation_method == "LNP") {
    CI_lower <- qlnorm(CI_tail_prob, gp_mean, gp_sd, lower.tail = TRUE)
    CI_upper <- qlnorm(CI_tail_prob, gp_mean, gp_sd, lower.tail = FALSE)
  } else if(transformation_method == "truncated") {
    CI_lower <- qtruncnorm(CI_tail_prob, a = 0, b = Inf, mean = gp_mean, sd = gp_sd)
    CI_upper <- qtruncnorm(1 - CI_tail_prob, a = 0, b = Inf, mean = gp_mean, sd = gp_sd)
  } else if(transformation_method == "rectified") {
    CI_lower <- quantile_rectified_norm(CI_tail_prob, mean = gp_mean, sd = gp_sd)
    CI_upper <- quantile_rectified_norm(1 - CI_tail_prob, mean = gp_mean, sd = gp_sd)
  } else {
    stop("Invalid transformation method: ", transformation_method)
  }
 
  return(list(lower = CI_lower, 
              upper = CI_upper))
   
}

quantile_zero_trunc_norm <- function(p, mean = 0, sd = 1, lower.tail = TRUE) {
  # This performs the same calculation as `qtruncnorm(p, a = 0, b = Inf, mean = mean, sd = sd)`
  # in the `truncnorm` package. Computes the quantile for the zero-truncated Gaussian. 
  #
  # Args:
  #    p: numeric(1), probability between 0 and 1. A probability of the form p = P(X <= t), for which this 
  #       function will return the value t. 
  #    mean: numeric(), vector of means of the underlying Gaussian distributions (not the means of the truncated
  #          Gaussians). 
  #    sd: numeric(), vector of standard deviations of the underlying Gaussian distributions. 
  #    lower.tail: logical(1), if TRUE the probability p is of the form P(X <= t). Otherwise, it is interpreted as 
  #                P(X > t). Note that: 
  #                quantile_zero_trunc_norm(p, lower.tail = TRUE) == quantile_zero_trunc_norm(1 - p, lower.tail = FALSE). 
  #
  # Returns: 
  #    numeric(), vector of equal length as `mean` and `sd` containing the evaluations of the truncated Gaussian quantile function. 
  
  if((p < 0) || (p > 1)) stop("p = ", p, " must be in range [0, 1].")
  
  if(!lower.tail) {
    p <- 1 - p
  }
  
  P <- pnorm(-mean/sd)
  sd * qnorm((1 - P) * p + P) + mean
  
}


quantile_rectified_norm <- function(p, mean = 0, sd = 1, lower.tail = TRUE) {
  # The quantile function for the rectified Gaussian distribution. That is, 
  # the quantile distribution for the random variable Y = max(0, X), where 
  # X ~ N(mean, sd). 
  #
  # Args:
  #    p: numeric(1), probability between 0 and 1. A probability of the form p = P(X <= t), for which this 
  #       function will return the value t. 
  #    mean: numeric(), vector of means of the underlying Gaussian distributions (not the means of the truncated
  #          Gaussians). 
  #    sd: numeric(), vector of standard deviations of the underlying Gaussian distributions. 
  #    lower.tail: logical(1), if TRUE the probability p is of the form P(X <= t). Otherwise, it is interpreted as 
  #                P(X > t). 
  #
  # Returns:
  #    numeric(), vector of equal length as `mean` and `sd` containing the evaluations of the rectified Gaussian quantile function. 
  
  if((p < 0) || (p > 1)) stop("p = ", p, " must be in range [0, 1].")
  
  if(!lower.tail) {
    p <- 1 - p
  }
  
  # Probability that the Gaussian is negative. 
  prob_neg <- pnorm(0, mean = mean, sd = sd)
  zero_selector <- (p <= prob_neg)
  
  # Compute quantiles.
  quantiles_rect_norm <- vector(mode = "numeric", length = length(mean))
  quantiles_rect_norm[zero_selector] <- 0
  quantiles_rect_norm[!zero_selector] <- qnorm(p, mean = mean[!zero_selector], sd = sd[!zero_selector])
  
  return(quantiles_rect_norm)
  
}


# ------------------------------------------------------------------------------
# Design Points
# ------------------------------------------------------------------------------

get_input_output_design <- function(N_points, theta_prior_params, ref_pars, pars_cal_sel, data_obs, PAR_data, output_vars, 
                                    scale_inputs, normalize_response, transformation_method = NA_character_, design_method = "LHS", 
                                    order_1d = FALSE) {
  #    normalize_response: logical(1), if TRUE will compute a Z-score transformation of the response/output variable. Will return the normalized 
  #                        variable in addition to the unnormalized one and will return the information used for the normalization so that the 
  #                        transformation can be inverted.
  #    
  
  
}


get_input_design <- function(N_points, theta_prior_params, design_method, scale_inputs, param_ranges = NULL, order_1d = FALSE, tail_prob_excluded = 0.01) {
  # Generates a set of points in parameter space to be used in GP emulation. This can either be the training inputs or 
  # a validation/test set. Optionally scales the points to the unit hypercube. If `scale_inputs` is TRUE and `param_ranges` is NULL
  # then uses the range of the input points in each dimension to perform this scaling. If `scale_inputs` is TRUE and `param_ranges` is provided, then 
  # `param_ranges` is used to perform the scaling. The common use case for this is to scale training data using the ranges of the data, then use these 
  # same ranges when scaling validation data. 
  #
  # Args:
  #    N_points: integer(1), the number of input points to generate. 
  #    theta_prior_params: data.frame containing the prior distribution information of the input 
  #                        parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                        for the requirements of this data.frame. 
  #    design_method: character(1), the algorithm used to generate the inputs. Currently supports "LHS" or "grid". 
  #    scale_inputs: logical(1), if TRUE will linearly scale the samples to the unit hypercube. Will return the scaled 
  #             sample in addition to the un-scaled one, and will also return the information used for the scaling 
  #             so that the transformation can be inverted. 
  #    param_ranges: matrix, of shape 2 x d, where d is the dimension of the parameter space. The rows correspond 
  #                  to the respective rows in `theta_prior_params`. The first and second columns are the lower and
  #                  upper bounds on the sample ranges of each parameter, respectively. Default is NULL, in which case 
  #                  no additional bounds are imposed on the samples, other than those already provided by the prior 
  #                  distributions.
  #    order_1d: logical(1), only relevant if the dimension of the input space (i.e. the number of 
  #              calibration parameters) is one-dimensional. In this case, if `order_1d` is TRUE then 
  #              the design points will be sorted in increasing order, with the training response 
  #              values ordered accordingly. This is convenient when plotting 1d GP plots. 
  #    tail_prob_excluded: numeric(), this is only relevant in certain cases, such as a Gaussian prior with "grid" design method. In this case 
  #                        the Gaussian has infinite support, but the grid method requires bounded support. If `tail_prob_excluded` is 0.01 then 
  #                        these bounds will be set to the .5% and 99.5% quantiles of the Gaussian. 
  #
  # Returns: 
  #    list, at the minimum this will have one named element "inputs", which is a N_points x d matrix containing the 
  #    input points that were generated. If `scale_inputs` is true, the list will also contain "inputs_scaled", the 
  #    scaled version of `inputs`, and `input_bounds`, a d x 2 matrix storing the range of each input dimension prior 
  #    to scaling. 
  
  # Generate train and test input datasets. 
  if(design_method == "LHS") {
    X_design <- get_LHS_design(N_points, theta_prior_params, param_ranges = param_ranges, order_1d = order_1d)
  } else if(design_method == "grid") {
    X_design <- get_grid_design(N_points, theta_prior_params, param_ranges = param_ranges)
  } else {
    stop("Invalid design method: ", design_method)
  }
  output_list <- list(inputs = X_design)
  
  # Scale inputs. 
  if(scale_inputs && is.null(param_ranges)) {
    output_list[c("inputs_scaled", "input_bounds")] <- prep_GP_training_data(X = X_design, scale_X = TRUE)[c("X", "input_bounds")]
  } else if(scale_inputs && !is.null(param_ranges)) {
    output_list[["inputs_scaled"]] <- scale_input_data(X_design, param_ranges)
  }

  return(output_list)
  
}


get_LHS_design <- function(N_points, theta_prior_params, param_ranges = NULL, order_1d = FALSE) {
  # Returns a Latin Hypercube Design, with distributions given by `theta_prior_params`. Can optionally provide 
  # upper and lower bounds on the samples. Note that if bounds are specified, then the points may not be 
  # exactly distributed as specified in `theta_prior_params`; e.g. Gaussian priors will become truncated Gaussian 
  # priors, and the lower/upper bounds on uniform priors may change. 
  # 
  # Args:
  #    N_points: integer(1), the number of Latin Hypercube samples. 
  #    theta_prior_params: data.frame containing the prior distribution information of the input 
  #                        parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                        for the requirements of this data.frame. 
  #    param_ranges: matrix, of shape d x 2, where d is the dimension of the parameter space. The rows correspond 
  #                  to the respective rows in `theta_prior_params`. The first and second columns are the lower and
  #                  upper bounds on the sample ranges of each parameter, respectively. Default is NULL, in which case 
  #                  no additional bounds are imposed on the samples, other than those already provided by the prior 
  #                  distributions. 
  #    order_1d: logical(1), only relevant if the dimension of the input space (i.e. the number of 
  #              calibration parameters) is one-dimensional. In this case, if `order_1d` is TRUE then 
  #              the design points will be sorted in increasing order, with the training response 
  #              values ordered accordingly. This is convenient when plotting 1d GP plots. 
  #
  # Returns:
  #    matrix, of dimension N_points x d, the Latin Hypercube sample. 
  
  # The dimension of the input space.
  d <- nrow(prior_params)
  
  if(is.null(param_ranges)) {
    param_ranges <- matrix(NA, nrow = 2, ncol = d)
  }
  
  # Generate LHS design. 
  X_LHS <- randomLHS(N_points, d)
  
  # Apply inverse CDS transform using prior distributions.
  for(j in seq_len(d)) {
    if(theta_prior_params[j, "dist"] == "Uniform") {
      lower_bound <- ifelse(is.na(param_ranges[1, j]), theta_prior_params[j, "param1"], param_ranges[1, j])
      upper_bound <- ifelse(is.na(param_ranges[2, j]), theta_prior_params[j, "param2"], param_ranges[2, j])
      X_LHS[,j] <- qunif(X_LHS[,j], lower_bound, upper_bound) 
    } else if(theta_prior_params[j, "dist"] == "Gaussian") {
      lower_bound <- ifelse(is.na(param_ranges[1, j]), -Inf, param_ranges[1, j])
      upper_bound <- ifelse(is.na(param_ranges[2, j]), Inf, param_ranges[2, j])
      X_LHS[,j] <- qtruncnorm(X_LHS[,j], a = lower_bound, b = upper_bound, mean = theta_prior_params[j, "param1"], sd = theta_prior_params[j, "param2"])
    } else {
      stop("Unsupported prior distribution: ", theta_prior_params[j, "dist"])
    }
  }
  
  # For 1 dimensional data, order samples in increasing order. 
  if(order_1d && (d == 1)) {
    X_LHS <- X_LHS[order(X_LHS),,drop=FALSE]
  }
  
  return(X_LHS)
  
}


get_grid_design <- function(N_points, theta_prior_params, param_ranges = NULL, tail_prob_excluded = 0.01) {
  # Generates a grid of inputs points. The prior distributions defined in `theta_prior_params` are only used to define the bounds of the grid 
  # in each dimension; they do not affect the distribution of input points generated, unlike in `get_LHS_design()`. Can optionally provide 
  # upper and lower bounds on the input points, which will override the parameter ranges defined in `theta_prior_params`. For cases where the 
  # prior has infinite support, see `tail_prob_excluded` in `Args` below.   
  #
  # Args:
  #    N_points: integer(1), the number of Latin Hypercube samples. 
  #    theta_prior_params: data.frame containing the prior distribution information of the input 
  #                        parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                        for the requirements of this data.frame. 
  #    param_ranges: matrix, of shape d x 2, where d is the dimension of the parameter space. The rows correspond 
  #                  to the respective rows in `theta_prior_params`. The first and second columns are the lower and
  #                  upper bounds on the sample ranges of each parameter, respectively. Default is NULL, in which case 
  #                  no additional bounds are imposed on the samples, other than those already provided by the prior 
  #                  distributions. 
  #    tail_prob_excluded: numeric(), this is only relevant in certain cases, such as a Gaussian prior with "grid" design method. In this case 
  #                        the Gaussian has infinite support, but the grid method requires bounded support. If `tail_prob_excluded` is 0.01 then 
  #                        these bounds will be set to the .5% and 99.5% quantiles of the Gaussian. 
  #
  # Returns:
  #    matrix, of dimension N_points x d, the grid points. 
  
  # The dimension of the input space.
  d <- nrow(prior_params)
  
  if(is.null(param_ranges)) {
    param_ranges <- matrix(NA, nrow = 2, ncol = d)
  }
  
  # Number of points marginally for each dimension.
  N_per_dim <- N_points^(1/d)
  if(round(N_per_dim) != N_per_dim) {
    stop("N_points must have an integer d^th root; d = ", d)
  }
  
  # Grid of inputs
  X_marginals <- matrix(nrow = N_per_dim, ncol = d)
  
  for(j in seq_len(d)) {
    if(theta_prior_params[j, "dist"] == "Uniform") {
      lower_bound <- ifelse(is.na(param_ranges[1, j]), theta_prior_params[j, "param1"], param_ranges[1, j])
      upper_bound <- ifelse(is.na(param_ranges[2, j]), theta_prior_params[j, "param2"], param_ranges[2, j])
      X_marginals[,j] <- seq(lower_bound, upper_bound, length.out = N_per_dim)
    } else if(theta_prior_params[j, "dist"] == "Gaussian") {
      lower_bound <- ifelse(is.na(param_ranges[1, j]), qnorm(tail_prob_excluded/2, theta_prior_params[j, "param1"], theta_prior_params[j, "param2"]), 
                            param_ranges[1, j])
      upper_bound <- ifelse(is.na(param_ranges[2, j]), qnorm(tail_prob_excluded/2, theta_prior_params[j, "param1"], theta_prior_params[j, "param2"], lower.tail = FALSE), 
                            param_ranges[2, j])
      X_marginals[,j] <- seq(lower_bound, upper_bound, length.out = N_per_dim)
    } else {
      stop("Unsupported prior distribution: ", theta_prior_params[j, "dist"])
    }
  }

  # Create grid.
  X_grid <- as.matrix(expand.grid(lapply(seq(1, d), function(j) X_marginals[,j])))
  colnames(X_grid) <- rownames(theta_prior_params)

  return(X_grid)
  
}


get_train_test_data <- function(N_train, N_test, prior_params, extrapolate, ref_pars, pars_cal_sel,
                                data_obs, PAR, output_vars, scale_X, normalize_Y, log_SSR, 
                                joint_LHS = FALSE, method = "LHS", order_1d = FALSE, true_theta_samples = NULL) {
  # Generates training (design) dataset and test dataset via Latin Hypercube sampling or 
  # via a grid-based approach. Note that the LHS method takes into account the prior, while 
  # the grid-based approach simply creates a grid of evenly spaced points within the bounds
  # defined by the prior. 
  # Input points are sampled from the space of calibration parameters, while the outputs 
  # are the squared L2 errors between the observed data and VSEM outputs, or the log 
  # of this quantity. Optionally pre-processes the data by scaling the inputs and 
  # normalizing the outputs. 
  #
  # Args:
  #    N_train: Number of design points. 
  #    N_test: Number of test points. 
  #    prior_params: data.frame containing the prior distribution information of the input 
  #                  parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                  for the requirements of this data.frame. 
  #    extrapolate: logical, if TRUE no truncation is performed for the test samples. Otherwise 
  #                 truncation is performed to guarantee the samples are within the extent of the 
  #                 training set. Default is TRUE. 
  #    ref_pars: data.frame, rownames should correspond to parameters of computer model. 
  #             Must contain column named "best". Parameters that are fixed at nominal 
  #             values (not calibrated) are set to their values given in the "best" column. 
  #    par_cal_sel: integer vector, selects the rows of 'ref_pars' that correspond to 
  #                 parameters that will be calibrated. 
  #    data_obs: matrix, dimension n x p (n = length of time series, p = number outputs).
  #              Colnames set to output variable names. 
  #    output_vars: character vector, used to the select the outputs to be considered in 
  #                 the likelihood; e.g. selects the correct sub-matrix of 'Sig_eps' and the 
  #                 correct columns of 'data_obs'. 
  #    PAR: numeric vector, time series of photosynthetically active radiation used as forcing
  #         term in VSEM.
  #    log_SSR: logical(1), if TRUE includes log-transformed outputs as well, in addition to the
  #             outputs on the original scale. 
  #    joint_LHS: logical, only relevant if `method` is LHS. If TRUE jointly samples train/test in a single LHS 
  #               sample. Othwerwise, samples the train and test data separately. Default is FALSE. 
  #    method: character(1), character string indicating the sampling or deterministic method used to 
  #            generate the train/test points. Valid options are "LHS" (latin hypercube sampling, the default), 
  #            or "grid". 
  #    order_1d: logical(1), only relevant if the dimension of the input space (i.e. the number of 
  #              calibration parameters) is one-dimensional. In this case, if `order_1d` is TRUE then 
  #              the design points will be sorted in increasing order, with the training response 
  #              values ordered accordingly. This is convenient when plotting 1d GP plots. 
  #    true_theta_samples: matrix, containing samples from the true posterior over the calibration parameters. If this is provided 
  #                        (non-NULL) it overrides the other test data settings. The train data is unaffected. 
  #
  # Returns:
  #    list, with names "X_train", "X_test", "X_train_preprocessed", "X_test_preprocessed", 
  #    "Y_train", "Y_train_preprocessed", "Y_test", "input_bounds", and "output_stats".
  #    The latter two encode the transformations applied to the inputs and outputs, respectively 
  #    (see prep_GP_training_data()). 
  
  # Generate train and test input datasets
  if(method == "LHS") {
    X_train_test <- LHS_train_test(N_train, N_test, prior_params, joint_LHS, extrapolate, order_1d, true_theta_samples = true_theta_samples)
  } else if(method == "grid") {
    X_train_test <- grid_train_test(N_train, N_test, prior_params, extrapolate, order_1d, true_theta_samples = true_theta_samples)
  }
  X_train <- X_train_test$X_train
  X_test <- X_train_test$X_test
  
  # Run VSEM on train and test sets to obtain outputs
  model_outputs_list_train <- lapply(X_train, function(theta) run_VSEM(theta, ref_pars, pars_cal_sel, PAR, output_vars))
  model_outputs_list_test <- lapply(X_test, function(theta) run_VSEM(theta, ref_pars, pars_cal_sel, PAR, output_vars))
  Y_train <- calc_SSR(data_obs[, output_vars], model_outputs_list_train, na.rm = TRUE)
  Y_test <- calc_SSR(data_obs[, output_vars], model_outputs_list_test, na.rm = TRUE)
  
  # Pre-process training data: scale inputs and normalize outputs
  GP_data_preprocessed <- prep_GP_training_data(X_train, Y_train, scale_X, normalize_Y)
  X_train_preprocessed <- GP_data_preprocessed$X
  Y_train_preprocessed <- GP_data_preprocessed$Y
  
  # Pre-process testing input data
  X_test_preprocessed <- scale_input_data(X_test, GP_data_preprocessed$input_bounds)
  
  # Pre-process testing output data
  Y_test_preprocessed <- normalize_output_data(Y = Y_test, output_stats = GP_data_preprocessed$output_stats)
  
  output_list <- list(X_train = X_train, X_test = X_test, X_train_preprocessed = X_train_preprocessed, 
                      X_test_preprocessed = X_test_preprocessed, Y_train = Y_train, 
                      Y_train_preprocessed = Y_train_preprocessed, Y_test = Y_test,
                      Y_test_preprocessed = Y_test_preprocessed,
                      input_bounds = GP_data_preprocessed$input_bounds, 
                      output_stats = GP_data_preprocessed$output_stats)
  
  # Include log-transformed L2 error
  if(log_SSR) {
    output_list[["log_Y_train"]] <- log(output_list$Y_train)
    output_list[["log_Y_test"]] <- log(output_list$Y_test)
    output_list[c("log_Y_train_preprocessed", "log_output_stats")] <- prep_GP_training_data(Y = output_list$log_Y_train, normalize_Y = TRUE)[c("Y", "output_stats")]
    output_list[["log_Y_test_preprocessed"]] <- normalize_output_data(Y = output_list$log_Y_test, output_stats = output_list$log_output_stats)
  }
  
  return(output_list)
  
}


LHS_train_test <- function(N_train, N_test, prior_params, joint = TRUE, extrapolate = TRUE, order_1d = FALSE, true_theta_samples = NULL) {
  # Generate a set training (design) points and test/validation points for evaluating a Gaussian 
  # process fit using Latin Hypercube Sampling (LHS). Can generate the train and test sets jointly, 
  # meaning they are sampled together in the LHS procedure. Or they can be sampled independently 
  # using two separate LHS schemes. Also allows the option to ensure that the test points are 
  # within the extent of the training points, if only the interpolation accuracy of the GP 
  # is of interest. In this case, the truncation means that the test points will not exactly
  # be distributed as specified in `prior_params`.
  #
  # Args:
  #    N_train: integer, the number of training points to sample. 
  #    N_test: integer, the number of test points to sample. 
  #    prior_params: data.frame containing the prior distribution information of the input 
  #                  parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                  for the requirements of this data.frame. 
  #    joint: logical, if TRUE jointly samples train/test. Default is TRUE. 
  #    extrapolate: logical, if TRUE no truncation is performed for the test samples. Otherwise 
  #                 truncation is performed to guarantee the samples are within the extent of the 
  #                 training set. Default is TRUE. 
  #    order_1d: logical(1), only relevant if the dimension of the input space (i.e. the number of 
  #              calibration parameters) is one-dimensional. In this case, if `order_1d` is TRUE then 
  #              the design points will be sorted in increasing order. This is convenient when 
  #              plotting 1d GP plots. 
  #    true_theta_samples: matrix, containing samples from the true posterior over the calibration parameters. If this is provided 
  #                        (non-NULL) it overrides the other test data settings. The train data is unaffected. 
  #
  # Returns: 
  #    list, with named elements "X_train" and "X_test", storing the N_train x d and 
  #    N_test x d matrices of samples, respectively. 
  
  # The dimension of the input space
  d <- nrow(prior_params)
  
  # Generate train and test sets
  if(joint) {
    X_combined <- randomLHS(N_train + N_test, d)
    X_train <- X_combined[1:N_train,,drop=FALSE]
    X_test <- X_combined[-(1:N_train),,drop=FALSE]
  } else {
    X_train <- randomLHS(N_train, d)
    X_test <- randomLHS(N_test, d)
  }
  
  # Apply inverse CDS transform using prior distributions.
  for(j in seq_len(d)) {
    
    if(prior_params[j, "dist"] == "Uniform") {
      
      # Train
      X_train[,j] <- qunif(X_train[,j], prior_params[j, "param1"], prior_params[j, "param2"])
      
      # Test
      if(extrapolate) {
        lower_bound_test <- prior_params[j, "param1"]
        upper_bound_test <- prior_params[j, "param2"]
      } else {
        lower_bound_test <- min(X_train[,j])
        upper_bound_test <- max(X_train[,j])
      }
      X_test[,j] <- qunif(X_test[,j], lower_bound_test, upper_bound_test)
      
    } else if(prior_params[j, "dist"] == "Gaussian") {
      
      # Train
      X_train[,j] <- qnorm(X_train[,j], mean = prior_params[j, "param1"], sd = prior_params[j, "param2"])
      
      # Test
      if(extrapolate) {
        lower_bound_test <- -Inf
        upper_bound_test <- Inf
      } else {
        lower_bound_test <- min(X_train[,j])
        upper_bound_test <- max(X_train[,j])
      }
      X_test[,j] <- qtruncnorm(X_test[,j], a = lower_bound_test, b = upper_bound_test, mean = prior_params[j, "param1"], sd = prior_params[j, "param2"])
    }
    
  }
  
  if(!is.null(true_theta_samples)) {
    X_test <- true_theta_samples
  }
  
  if(order_1d && (d == 1)) {
    X_train <- X_train[order(X_train),,drop=FALSE]
    X_test <- X_test[order(X_test),,drop=FALSE]
  }
  
  return(list(X_train = X_train, X_test = X_test))
  
}


grid_train_test <- function(N_train, N_test, prior_params, extrapolate, tail_prob_excluded = 0.01, true_theta_samples = NULL) {
  # Note that currently the "extrapolate" argument makes no difference, since grid points are placed
  # at the edges of the defined bounds. This can be changed if needed, by adding an "extrapolate_prob" 
  # argument or something like that. 
  
  # The dimension of the input space.
  d <- nrow(prior_params)
  
  # Number of points marginally for each dimension.
  N_train_per_dim <- N_train^(1/d)
  N_test_per_dim <- N_test^(1/d)
  if((round(N_train_per_dim) != N_train_per_dim) || (round(N_test_per_dim) != N_test_per_dim)) {
    stop("N_train and N_test must have an integer d^th root; d = ", d)
  }
  
  # Grid of inputs
  X_train_marginals <- matrix(nrow = N_train_per_dim, ncol = d)
  X_test_marginals <- matrix(nrow = N_test_per_dim, ncol = d)
  
  for(j in seq_len(d)) {
    
    if(prior_params[j, "dist"] == "Uniform") {
      
      # Train
      X_train_marginals[,j] <- seq(prior_params[j, "param1"], prior_params[j, "param2"], length.out = N_train_per_dim)
      
      # Test
      if(extrapolate) {
        lower_bound_test <- prior_params[j, "param1"]
        upper_bound_test <- prior_params[j, "param2"]
      } else {
        lower_bound_test <- min(X_train_marginals[,j])
        upper_bound_test <- max(X_train_marginals[,j])
      }
      X_test_marginals[,j] <- seq(lower_bound_test, upper_bound_test, length.out = N_test_per_dim)
      
    } else if(prior_params[j, "dist"] == "Gaussian") {
      
      # Train
      left_bound <- qnorm(tail_prob_excluded/2, prior_params[j, "param1"], prior_params[j, "param2"])
      right_bound <- left_bound + 2*abs(left_bound)
      X_train_marginals[,j] <- seq(left_bound, right_bound, length.out = N_train_per_dim)
      
      # Test
      if(extrapolate) {
        lower_bound_test <- left_bound
        upper_bound_test <- right_bound
      } else {
        lower_bound_test <- min(X_train_marginals[,j])
        upper_bound_test <- max(X_train_marginals[,j])
      }
      X_test_marginals[,j] <- seq(lower_bound_test, upper_bound_test, length.out = N_test_per_dim)
    }
    
  }
  
  # Create grids
  X_train <- as.matrix(expand.grid(lapply(seq(1, d), function(j) X_train_marginals[,j])))
  
  if(is.null(true_theta_samples)) {
    X_test <- as.matrix(expand.grid(lapply(seq(1, d), function(j) X_test_marginals[,j])))
  } else {
    X_test <- true_theta_samples
  }
  colnames(X_train) <- rownames(theta_prior_params)
  colnames(X_test) <- rownames(theta_prior_params)
  
  return(list(X_train = X_train, X_test = X_test))
  
}








