#
# gp_emulator_functions.r
#

library(scoringRules)

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
    rownames(input_bounds) <- c("min", "max")
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
  
  cols <- colnames(X)
  
  if(inverse) {
    X <- X %*% diag(input_bounds[2,] - input_bounds[1,], ncol(X)) + matrix(input_bounds[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  } else {
    X <- (X - matrix(input_bounds[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(input_bounds[2,] - input_bounds[1,]), ncol(X))
  }
  
  colnames(X) <- cols
  
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

 
# TODO: Parallelize independent GP fitting.
fit_independent_GPs <- function(X_train, Y_train, gp_lib, gp_kernel) {
  # Builds on top of `fit_GP()` to generalize to multivariate GP regression. Simply 
  # fits independent GPs for each output (where the same input points are used for each GP).
  #
  # Args:
  #    X_train: matrix of shape N x d, where N is the number of design (training) points
  #             and d is the dimension of the input space. 
  #    Y_train: matrix of shape N x p, with jth column containing the training outputs 
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


update_GP <- function(gp_obj, gp_lib, X_new, y_new = NULL, update_hyperparameters = FALSE) {
  # Updates a GP fit by conditioning on new data (X_new, y_new). Currently this is only 
  # supported for `gp_lib == "hetGP"`. Currently, this updates GP hyperparameter fit and 
  # then adds (X_new, y_new) to the design set, so that future predictive mean/variance 
  # calculations will condition on (X_new, y_new) in addition to the previous design. 
  # Note that `hetGP` also supports conditioning on the new data without updating 
  # GP hyperparameters. This function could be expanded to allow for this if desired. 
  # The new design data `X_new` and `y_new` are assumed to be properly scaled/normalized. 
  # The function `update_independent_GPs()` offers the option to perform pre-processing on 
  # the new data. 
  #
  # Args:
  #    gp_obj: An object representing a GP fit, which will differ based on the library used to fit the GP. 
  #    gp_lib: character(1), the library used to fit the GP. Currently supports "hetGP". 
  #    X_new: matrix, the new design inputs, of shape M x D where M is the number of new inputs and 
  #           D is the dimension of the input space. Note that these inputs should be properly scaled,
  #           as this function does not pre-process the new design data. 
  #    y_new: numeric, vector of length M. The outputs corresponding to the M new design points. 
  #           Note that these inputs should be properly normalized, as this function does not pre-process 
  #           the new design data. If NULL, then only updates the GP predictive variance. This is equivalent 
  #           to conditioning the GP on (X_new, mu(X_new)), where mu() is the GP predictive mean. 
  #    update_hyperparameters: If TRUE, runs optimization to re-fit GP hyperparameters. Otherwise, 
  #                            the hyperparamters are not changed. 
  #
  # Returns:
  #    The updated GP object. 
  
  if(gp_lib == "hetGP") {
    if(is.null(y_new)) y_new <- rep(NA, nrow(X_new))
    
    if(update_hyperparameters) gp_obj <- update(object = gp_obj, Xnew = X_new, Znew = y_new)
    else gp_obj <- update(object = gp_obj, Xnew = X_new, Znew = y_new, maxit = 0)
    
  } else {
    stop("Currently `update_GP` only supports GP library `hetGP`")
  }
  
  return(gp_obj)
  
}


update_independent_GPs <- function(gp_fits, gp_lib, X_new, Y_new = NULL, update_hyperparameters = FALSE, input_bounds = NULL, output_stats = NULL) {
  # A wrapper around `update_GP()` which udpates a set of independent GPs given new design 
  # data `(X_new, Y_new)`. See `update_GP()` for more details. This function allows the option 
  # to pre-process the data (i.e. scale inputs and normalize outputs). 
  #
  # Args:
  #    gp_obj_list: A list of GP objects, each of which represents a GP fit to one of the outputs. The objects will 
  #                differ based on the specific GP library used for fitting. All of the objects in the list must have 
  #                been fit using the same library. 
  #    gp_lib: character(1), the library used to fit the GP. Currently supports "hetGP". 
  #    X_new: matrix, the new design inputs, of shape M x D where M is the number of new inputs and 
  #           D is the dimension of the input space. Note that these inputs should be properly scaled,
  #           as this function does not pre-process the new design data. 
  #    Y_new: matrix of shape M x P, with pth column containing the training outputs for the pth output variable
  #           corresponding to the new inputs `X_new`. If NULL, only updates GP predictive variances; see `update_GP()`.
  #    update_hyperparameters: If TRUE, runs optimization to re-fit GP hyperparameters. Otherwise, 
  #                            the hyperparamters are not changed. 
  #    input_bounds: matrix, of shape 2 x D. The columns correspond to the respective rows in `theta_prior_params`.
  #                  The first and second rows are the lower and upper bounds on the sample ranges of each parameter,
  #                  respectively. If not NULL, used to scale `X_new` prior to updating GPs. 
  #    output_stats: matrix of dimensions 2x1. The matrix must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of the output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). Note that `output_stats` 
  #                  is assumed to be on the correct scale; e.g. for the log-normal process, 
  #                  `output_stats` should be on the log-scale. If not NULL, used to normalize the outputs 
  #                  `Y_new` prior to updating GPs. 
  #
  # Returns:
  #    The updated list of GP objects. 
  
  # Pre-process data. 
  if(!is.null(input_bounds)) {
    X_new <- scale_input_data(X_new, input_bounds = input_bounds)
  }
  if(!is.null(output_stats) && !is.null(Y_new)) {
    Y_new <- normalize_output_data(Y_new, output_stats)
  }
  
  for(j in seq_along(gp_fits)) {
    gp_fits[[j]] <- update_GP(gp_fits[[j]], gp_lib, X_new, Y_new[,j], update_hyperparameters = update_hyperparameters)
  }
  
  return(gp_fits)
  
}


fit_emulator_design_list <- function(emulator_setting, design_list) {
  # A function for fitting a single GP emulator to different designs. Returns 
  # an "emulator_info_list", of length equal to the length of `design_list`. 
  # Each element of the list is itself a list with elements "gp_fits", "input_bounds", 
  # "output_stats", and "settings". 
  #
  # Args:
  #    emulator_setting: list, with elements "gp_lib", "kernel", "transformation_method" 
  #                      providing the GP specifications. 
  #    design_list: list of designs for emulator, as returned by `get_design_list()`. 
  #
  # Returns:
  #    list, the emulator info list as described above. The names of the list are set 
  #    as the names of `design_list`. 
  
  emulator_info_list <- list()
  is_LNP <- (emulator_setting$transformation_method == "LNP")
  output_stats_sel <- ifelse(is_LNP, "log_output_stats", "output_stats")
  outputs_normalized_sel <- ifelse(is_LNP, "log_outputs_normalized", "outputs_normalized")
  
  for(j in seq_along(design_list)) {
    label <- names(design_list)[j]
    
    gp_fits <- fit_independent_GPs(X_train = design_list[[j]]$inputs_scaled, 
                                   Y_train = design_list[[j]][[outputs_normalized_sel]], 
                                   gp_lib = emulator_setting$gp_lib, 
                                   gp_kernel = emulator_setting$kernel)$fits
    names(gp_fits) <- colnames(design_list[[j]]$outputs_normalized)
    
    emulator_info_list[[label]] <- list(gp_fits = gp_fits, 
                                        input_bounds = design_list[[j]]$input_bounds, 
                                        output_stats = design_list[[j]][[output_stats_sel]], 
                                        settings = emulator_setting)

  }
  
  return(emulator_info_list)
  
}


# ------------------------------------------------------------------------------
# Predicting with GPs
# ------------------------------------------------------------------------------

# TODO: 
#    - Currently sd2_nug should not be trusted; e.g. it is not corrected in the case of truncation/rectification. Might just be easier to 
#      drop this or combine it with the pointwise variance predictions. 
predict_GP <- function(X_pred, gp_obj, gp_lib, include_cov_mat = FALSE, denormalize_predictions = FALSE,
                       output_stats = NULL, transformation_method = NA_character_, return_df = FALSE, 
                       X_pred_2 = NULL) {
  # Calculate GP predictive mean, variance, and optionally covariance matrix at specified set of 
  # input points. 
  #
  # Args:
  #    X_pred: matrix of dimension N_pred x d, where d is the dimension of the input space. Each row is an input 
  #            point at which to predict. 
  #    gp_obj: An object representing a GP fit, which will differ based on the library used to fit the GP. 
  #    gp_lib: character(1), the library used to fit the GP. Currently supports "mlegp" or "hetGP". 
  #    include_cov_mat: logical(1), if TRUE, calculates and returns the N_pred x N_pred predictive covariance matrix 
  #                     over the set of input points. Otherwise, only calculates the pointwise predictive variances. 
  #    denormalize_predictions: logical(1), if TRUE, applies linear transformation to predictions, 
  #                             inverting the Z-score transformation.
  #    output_stats: If not NULL, then a matrix of dimensions 2x1. The matrix 
  #                  must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of the output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). Only required if `denormalize_predictions` is TRUE. Note that `output_stats` 
  #                  is assumed to be on the correct scale; e.g. for the log-normal process, `output_stats` should be on the log-scale. 
  #    transformation_method: character(1), if "default" treats the predictive GP distribution as Gaussian. If "rectified", treats
  #                           it as rectified Gaussian and if "truncated" treats it as truncated Gaussian (truncated at 0). If 
  #                           "LNP", exponentiates the GP, resulting in a log-normal process. 
  #    return_df: logical(1), if TRUE then the predictive means, variance, and nugget variances are returned in a single data.frame
  #               with column names "mean", "var", and "var_nug". The return type of this function is then a list with the first 
  #               element being this data.frame, and the second element "cov" being the predictive covariance matrix (NULL if 
  #               `include_cov_mat` is FALSE). Default is FALSE (return list with named elements "mean", "var", "var_nug", "cov").
  #    X_pred_2: matrix of dimension N_pred_2 x 2. A second set of input locations used to compute the predictive covariance 
  #              between  `X_pred` and `X_pred_2`. Only required if `include_cov_mat` is TRUE. To produce the predictive 
  #              covariance matrix at the inputs `X_pred`, then simply set this argument equal to `X_pred`. 
  #
  # Returns:
  #    list, potentially with named elements "mean", "var", "cov". If no transformation is applied, will also contain elements 
  #    "var_nug" and "var_comb" containing the nugget variance and combined variance (var + var_nug). If a transformation is applied, 
  #    "var_nug" will not be present, but "var_comb" will be. Depending on the tranformation, "var_comb" may no longer be the sum of 
  #    var and var_nug. "cov" will be NULL if `include_cov_mat` is FALSE. Otherwise, it will include the predictive covariance matrix. 
  #    The diagonal of this matrix equals "var", not "var_comb". To take into account the nugget variance as well, the diagonal of this 
  #    matrix can simply be replaced with "var_comb". Due to the potential for transformations, it is recommended to work with 
  #    "var" and "var_comb" instead of "var_nug". 
  #
  #    Alternatively returns a list of two elements, "df" and "cov" if `return_df` is TRUE. 
  #    See above argument `return_df` for details. 
  
  if(!denormalize_predictions && !is.na(transformation_method)) {
    message("Applying transformation to normalized predictions may not make sense; e.g. for truncated Gaussian.")
  }
  
  pred_list_names <- c("mean", "var", "var_nug", "cov")
  pred_list <- vector(mode = "list", length = length(pred_list_names))
  names(pred_list) <- pred_list_names
  
  # Generate GP Predictions. 
  if(gp_lib == "mlegp") {
    mlegp_pred <- predict(gp_obj, newData = X_pred, se.fit = TRUE)
    pred_list[["mean"]] <- mlegp_pred[["fit"]]
    pred_list[["var"]] <- mlegp_pred[["se.fit"]]^2
    pred_list[["var_nug"]] <- gp_obj$nugget
  } else if(gp_lib == "hetGP") {
    # Second matrix for computing predictive covariance
    if(include_cov_mat) {
      X_prime <- X_pred_2
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
    if(include_cov_mat) {
      pred_list[["cov"]] <- output_stats["var_Y",1] * pred_list[["cov"]]
    }
  }
  
  # Sum predictive variance of latent function value and nugget to obtain combined variance. 
  pred_list[["var_comb"]] <- pred_list$var + pred_list$var_nug
  
  # Apply transformation to GP predictions. 
  if(!is.na(transformation_method)) {
    pred_list_transformed <- transform_GP_predictions(pred_list$mean, pred_list$var, transformation_method = transformation_method, gp_cov = pred_list$cov)
    pred_list_transformed[["var_comb"]] <- transform_GP_predictions(pred_list$mean, pred_list$var_comb, transformation_method = transformation_method)$var
    pred_list <- pred_list_transformed
  }
  
  if(return_df) {
    pred_df <- data.frame(mean = pred_list[["mean"]], 
                          var = pred_list[["var"]], 
                          var_comb = pred_list[["var_comb"]])
    return(list(df = pred_df, cov = pred_list$cov))
  }
  
  return(pred_list)
  
}


# TODO: 
#    - should be able to return only mean or variance. 
#    - update predict_GP so that it can return data.frame in wide or long format. 
predict_independent_GPs <- function(X_pred, gp_obj_list, gp_lib, include_cov_mat = FALSE, denormalize_predictions = FALSE,
                                    output_stats = NULL, transformation_method = NA_character_, return_df = FALSE, 
                                    output_variables = NULL, X_pred_2 = NULL) {
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
  #    include_cov_mat: logical(1), if TRUE, calculates and returns the N_pred x N_pred predictive covariance matrix 
  #             over the set of input points. Otherwise, only calculates the pointwise predictive variances. 
  #    denormalize_predictions: logical(1), if TRUE, applies linear transformation to predictions, 
  #                             inverting the Z-score transformation.
  #    output_stats: If not NULL, then a matrix of dimensions 2xp, where p is the number of output variables. The matrix 
  #                  must have rownames "mean_Y" and "var_Y" storing the mean and 
  #                  variance of each output variable used to compute the Z-scores. This object is 
  #                  returned by prep_GP_training_data(). Only required if `denormalize_predictions` is TRUE. 
  #                  Note that `output_stats` is assumed to be on the correct scale; i.e. for log-normal process,
  #                  `output_stats` should be on log-scale. 
  #    transformation_method: character(1), specifies the transformation to apply. The current valid options
  #                           are "LNP" (exponentiates the GP, resulting in a log-normal process), "truncated"
  #                           (transforms the GP to a left-truncated (at 0) GP) and "rectified" (transforms the 
  #                           GP to a rectified Gaussian). 
  #    return_df: logical(1), if TRUE then the predictive means, variance, and nugget variances are returned in 
  #               a single data.frame The columns of the resulting data.frame are of the form  
  #               "mean_<output>" and "var_<output>", where `<output>` identifies the different 
  #               response variables/data constraints. The values of `<output>` will be set using the argument 
  #               `output_variables`. The return type of the function is then a list with two elements, the first 
  #               "df" being the data.frame just described, and the second "cov_list" a list of the predictive 
  #               covariance matrices (is `include_cov_mat` is FALSE then these will be NULL). Note that each element
  #               of the predictive covariance list is itself a list of predictive covariance matrices, with each 
  #               matrix corresponding to a single output of a multi-output GP. 
  #    output_variables: character(p), vector of the p output variable names used to label the predictions. 
  #                      If NULL, will use labels "output1", ..., "outputp". 
  # 
  # Returns:
  #    list, with length equal to the length of `gp_obj_list`. Each element of this list is itself a list, 
  #    which is the return value of the function `predict_GP()` applied 
  #    to each GP in `gp_obj_list`). Or if `return_df` is TRUE, a list of two elements "df" and "cov_list". 
  #    See `return_df` above for details. 
  
  if(is.null(output_variables)) {
    output_variables <- paste0("output", seq_along(gp_obj_list))
  }
  
  pred_list <- lapply(seq_along(gp_obj_list), function(j) predict_GP(X_pred, gp_obj_list[[j]], gp_lib, include_cov_mat, denormalize_predictions, 
                                                                     output_stats[,j,drop=FALSE], transformation_method = transformation_method, 
                                                                     return_df = return_df, X_pred_2 = X_pred_2))

  if(return_df) {
    pred_df_list <- lapply(pred_list, function(l) l$df)
    for(j in seq_along(pred_df_list)) colnames(pred_df_list[[j]]) <- paste(colnames(pred_df_list[[j]]), output_variables[j], sep = "_")
    df <- do.call("cbind", pred_df_list)
    cov_list <- lapply(pred_list, function(l) l$cov)
    return(list(df = df, cov_list = cov_list))
  }
  
  names(pred_list) <- output_variables
  return(pred_list)
  
}

# TODO:
#    - Need to be more careful that sig2_eps is ordered the same as the cols returned by predict_independent_GPs(). 
#    - Having to select cols "var_comb_output" or "var_output" is awkward. 
#    - Need to write function to compute sig_eps prior log density. 
predict_lpost_GP_approx <- function(theta_vals_scaled = NULL, theta_vals_unscaled = NULL, emulator_info_list,
                                    sig2_eps, theta_prior_params = NULL, N_obs = NULL, include_nugget = TRUE,  
                                    gp_pred_list = NULL, include_sig_eps_prior = FALSE, theta_vals_scaled_2 = NULL,
                                    return_vals = c("mean", "var")) {
  # In the loss emulation setting, where the squared error maps (SSR) are modeled as GPs, these GPs induce a 
  # random field approximation on the unnormalized log (conditional) posterior 
  # log pi(theta|Sigma) := log p(Y|theta, Sigma) + log pi_0(theta). This random field approximation is also 
  # a Gaussian process over the input space of calibration parameters `theta`. This function computes the predictive 
  # mean and variance variance of this random field representation of the unnormalized log (conditional) posterior  
  # at a discrete set of input locations `theta_vals_unscaled`. Note that this function assumes the predictive 
  # distribution of the GP emulators  is Gaussian; thus, it may differ slightly from a Monte Carlo estimate 
  # given that in reality the distribution  is truncated or rectified Gaussian due to the fact that negative 
  # sampled SSR values are not allowed.
  #
  # Args:
  #    theta_vals_scaled: matrix of dimension M x D of input locations (scaled to lie in unit hypercube) at which to sample
  #                       the log density values. Each row is a location. 
  #    theta_vals_unscaled: matrix, of dimension M x D; the same as `theta_vals_scaled` but the inputs are unscaled. The unscaled 
  #                         inputs are only needed if "mean" is in `return_vals` in order to compute the prior density. However, 
  #                         in this case, only one of either scaled or unscaled need be passed, as the scaling 
  #                         can be performed using the information in `emulator_info_list$input_bounds`. 
  #    emulator_info_list: list, the emulator info list. 
  #    sig2_eps: numeric, vector of length P containing the observation variance parameters. 
  #    theta_prior_params: data.frame containing the prior distribution information of the input 
  #                        parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                        for the requirements of this data.frame. Only required if "mean" is in `return_vals`. 
  #    include_nugget: logical, if TRUE includes the GP nugget variance in the variance calculations. 
  #    gp_pred_list: A list, as returned by `predict_independent_GPs()`. This allows the predictive means and 
  #                  variances of the underlying GP emulators evaluated at the input locations at the to be passed 
  #                  if they have already been computed. If NULL, they are computed here. 
  #    include_sig_eps_prior: If TRUE, sample from the unnormalized log joint posterior, instead of the conditional; 
  #                           see description above for clarification. 
  #    return_vals: character, either "mean", "var", or c("mean", "var") depending on whether the predictive mean
  #                 or variance is desired, or both. 
  #
  # Returns:
  #    numeric, vector of length M containing the predictive variance evaluations at the M inputs. 
  #
  
  if(is.null(theta_vals_scaled) && is.null(theta_vals_unscaled)) {
    stop("Either `theta_vals_scaled` or `theta_vals_unscaled` must be non-NULL.")
  } else if(is.null(theta_vals_scaled)) {
    theta_vals_scaled <- scale_input_data(theta_vals_unscaled, input_bounds = emulator_info_list$input_bounds)
  }
  
  # Mean and variance predictions for the underlying GP fits to SSR. 
  if(is.null(gp_pred_list)) {
    
    include_cov_mat <- FALSE
    if("cov" %in% return_vals) {
      if(is.null(theta_vals_scaled_2)) stop("`theta_vals_scaled_2` is required to compute predictive covariance.")
      include_cov_mat <- TRUE
    }
    
    gp_pred_info <- predict_independent_GPs(X_pred = theta_vals_scaled, 
                                            gp_obj_list = emulator_info_list$gp_fits,  
                                            gp_lib = emulator_info_list$settings$gp_lib, 
                                            denormalize_predictions = TRUE, 
                                            output_stats = emulator_info_list$output_stats, 
                                            include_cov_mat = include_cov_mat, 
                                            X_pred_2 = theta_vals_scaled_2, 
                                            return_df = TRUE)
    
    gp_pred_list <- gp_pred_info$df
    gp_pred_cov_list <- gp_pred_info$cov_list
    
  }

  return_list <- list(mean = NULL, var = NULL, cov = NULL)
  
  # Compute predictive mean of lpost approximation.
  if("mean" %in% return_vals) {
    if(any(is.null(theta_prior_params), is.null(N_obs))) stop("`theta_prior_params` and `N_obs` are required for predictive mean calculation.")
    if(is.null(theta_vals_unscaled)) theta_vals_unscaled <- scale_input_data(theta_vals_scaled, input_bounds = emulator_info_list$input_bounds, inverse = TRUE)
    
    scaled_means <- as.matrix(gp_pred_list[, grep("mean_output", colnames(gp_pred_list))]) %*% diag(1/sig2_eps)
    llik_pred_mean <- -0.5 * sum(N_obs * log(2*pi*sig2_eps)) - 0.5 * rowSums(scaled_means)
    return_list$mean <- llik_pred_mean + calc_lprior_theta(theta_vals_unscaled, theta_prior_params)
  }
  
  # Compute predictive variance of lpost approximation. 
  if("var" %in% return_vals) {
    col_pattern_var <- ifelse(include_nugget, "var_comb_output", "var_output")
    scaled_vars <- as.matrix(gp_pred_list[, grep(col_pattern_var, colnames(gp_pred_list))]) %*% diag(1/sig2_eps^2)
    lpost_pred_var <- 0.25 * rowSums(scaled_vars)
    return_list$var <- lpost_pred_var
  }
  
  # Compute predictive covariance of lpost values between inputs `theta_vals_scaled` and `theta_vals_scaled_2`
  if("cov" %in% return_vals) {
    gp_pred_cov_list <- lapply(seq_along(gp_pred_cov_list), function(p) gp_pred_cov_list[[p]] / sig2_eps[p]^2)
    return_list$cov <- 0.25 * Reduce("+", gp_pred_cov_list)
  }
  
  # Optionally include prior on likelihood variance parameters, in which case `lpost` refers to the point unnormalized 
  # posterior pi(theta, Sigma), rather than the conditional posterior pi(theta|Sigma).
  # Note: this is only used in the predictive mean, not variance. 
  if(include_sig_eps_prior) {
    stop("Inclusion of sig eps prior not yet implemented.") 
  }
  
  return(return_list)
  
} 


get_emulator_comparison_pred_df <- function(emulator_info_list, X_test, scale_inputs = FALSE,
                                            denormalize_predictions = TRUE, output_variables = NULL, 
                                            include_cov_mat = FALSE, transform_predictions = FALSE) {
  # A convenience function for computing different GP emulator predictions on the same set of validation points. 
  # This is mainly used for testing when comparing different emulators. This function loops over each emulator, 
  # returning the GP predictive means and variances at the test inputs `X_test`. All predictions are compiled 
  # in a single data.frame. GP predictions will be raw predictions if `transform_predictions` is FALSE; 
  # otherwise, the predictions will be transformed according to the specified GP predictive distribution 
  # (e.g. log-normal, truncated Gaussian, rectified Gaussian). In this case, the transformations are determined 
  # by `emulator_info_list`. 
  #
  # Args:
  #    emulator_info_list: Named list, where each element is an "emulator info" object. The emulator info objects 
  #                        are themselves lists with elements "gp_fits", "input_bounds", "output_stats", and "settings". 
  #    X_test: matrix, of dimension M x D, where M is the number of test inputs and D the dimension of the parameter space. 
  #            The test inputs at which to predict. 
  #    scale_inputs: logical(1), if TRUE, uses the "input_bounds" element of the emulator info lists to scale `X_test` prior to prediction. 
  #                  If FALSE, assumes `X_test` is already properly scaled. 
  #    denormalize_predictions: logical(1), if TRUE the GP predictions are returned on the original scale by inverting the normalization 
  #                             transformation. Otherwise, the GP predictions are returned on the normalized scale. 
  #    output_variables: character, vector of output variable names, which is passed to `predict_independent_GPs()`. If NULL, 
  #                      `predict_independent_GPs()` just assigns numbers to the outputs. 
  #    include_cov_mat: logical(1), if TRUE, calculates and returns the N_pred x N_pred predictive covariance matrix 
  #             over the set of input points. Otherwise, only calculates the pointwise predictive variances.
  #    transform_predictions: logical(1), GP predictions will be raw predictions if `transform_predictions` is FALSE; 
  #                           otherwise, the predictions will be transformed according to the specified GP predictive distribution 
  #                           (e.g. log-normal, truncated Gaussian, rectified Gaussian). In this case, the transformations are determined 
  #                           by `emulator_info_list`.
  #
  # Returns:
  #    list, with two elements "df" and "cov_list". The former is a 
  #    data.frame, with columns of the form "<test_label>_<pred_type>_<output_variable>". "<test_label>" identifies 
  #    which emulator the prediction pertains to; these names are taken from the names in the list `emulator_info_list`. 
  #    "<pred_type>" is either "mean", "var", or "var_comb" for the predictive mean, predictive variance, and combined 
  #    variance (i.e. including nugget), respectively. "<output_variable>" is the name of the output variable to 
  #    which the prediction  corresponds. Additional columns are added for each emulator test with the output variable 
  #    "total_SSR" (sum of squared residuals). These columns simply sum over all output variables to yield the predictive 
  #    quantities for the total SSR. Note that this is valid for the variances since the multi-output GPs consist 
  #    of independent GPs for each output. The element "cov_list" is the list (of lists) of predictive covariance  
  #    matrices returned by `predict_independent_GPs()`; see this function as well as `predict_GP()` for more 
  #    details on the predictive covariance list. 

  pred_list <- vector(mode = "list", length = length(emulator_info_list))
  gp_cov_list <- vector(mode = "list", length = length(emulator_info_list))
  test_labels <- names(emulator_info_list)
  names(gp_cov_list) <- test_labels
  
  if(!scale_inputs) X_test_scaled <- X_test

  for(j in seq_along(emulator_info_list)) {
    test_label <- test_labels[j]
    
    # Scale inputs. 
    if(scale_inputs) X_test_scaled <- scale_input_data(X_test, input_bounds = emulator_info_list[[j]]$input_bounds)
    
    # Determine transformation, if transforming predictions. 
    transformation_method <- ifelse(transform_predictions, emulator_info_list[[j]]$settings$transformation_method, 
                                                           NA_character_)

    # Compute GP predictive distribution. 
    gp_test_pred_list <- predict_independent_GPs(X_pred = X_test_scaled, 
                                                 gp_obj_list = emulator_info_list[[j]]$gp_fits, 
                                                 gp_lib = emulator_info_list[[j]]$settings$gp_lib, 
                                                 include_cov_mat = include_cov_mat, 
                                                 denormalize_predictions = denormalize_predictions,
                                                 output_stats = emulator_info_list[[j]]$output_stats, 
                                                 return_df = TRUE, 
                                                 output_variables = output_variables,  
                                                 transformation_method = transformation_method)
    gp_test_pred_df <- gp_test_pred_list$df
    gp_cov_list[[test_label]] <- gp_test_pred_list$cov_list
    
    # Sum across columns to compute predictive mean and variance of total sum of squared residuals 
    # (summed across all output variables). 
    mean_cols <- grep("mean", colnames(gp_test_pred_df))
    var_comb_cols <- grep("var_comb", colnames(gp_test_pred_df))
    var_cols <- setdiff(grep("var", colnames(gp_test_pred_df)), var_comb_cols)
    gp_test_pred_df[["mean_total_SSR"]] <- rowSums(gp_test_pred_df[, mean_cols, drop = FALSE])
    gp_test_pred_df[["var_comb_total_SSR"]] <- rowSums(gp_test_pred_df[, var_comb_cols, drop = FALSE])
    gp_test_pred_df[["var_total_SSR"]] <- rowSums(gp_test_pred_df[, var_cols, drop = FALSE])
    
    # Prepend test label to column names. 
    colnames(gp_test_pred_df) <- paste(test_label, colnames(gp_test_pred_df), sep = "_")
    
    # Combine with current data.frame. 
    if(j == 1) gp_pred_df <- gp_test_pred_df
    else gp_pred_df <- cbind(gp_pred_df, gp_test_pred_df)

  }
  
  return(list(df = gp_pred_df, cov_list = gp_cov_list))
  
}


# ------------------------------------------------------------------------------
# Evaluating GPs
# ------------------------------------------------------------------------------


get_gp_rmse <- function(gp_mean, y_true, ...) {
  # Computes the root mean squared error between GP mean predictions and a true  
  # baseline response.
  #
  # Args:
  #    gp_mean: numeric, vector of GP mean predictions. 
  #    y_true: numeric, vector of true response values. Must be equal length as `gp_mean`.
  #
  # Returns:
  #    numeric, the RMSE. 
  
  return(sqrt(sum((gp_mean - y_true)^2) / length(y_true)))
  
}


get_gp_srmse <- function(gp_mean, gp_var, y_true, ...) {
  # Computes the standardized root mean squared error between GP mean predictions and a true  
  # baseline response. The pointwise error is (mean_i - y_i)^2 / var_i; these pointwise 
  # errors are averaged, and the square root of the average is returned. 
  #
  # Args
  #    gp_mean: numeric, vector of GP mean predictions. 
  #    gp_var: numeric, vector of GP variance predictions.
  #    y_true: numeric, vector of true response values. Must be equal length as `gp_mean`.
  #
  # Returns:
  #    numeric, the standardized RMSE. 
  
  return(sqrt(sum((gp_mean - y_true)^2 / gp_var) / length(y_true)))
  
}


get_gp_nlpd_pointwise <- function(gp_mean, gp_var, y_true, transformation_method = NA_character_, ...) {
  # Returns the average GP pointwise negative log predictive density of GP predictions. Note that this 
  # does not take into account GP predictive covariance; the log predictive density is computed independently 
  # for each test point and then an average is taken across all test points. 
  # This is essentially a wrapper function for `get_GP_pointwise_predictive_density(..., log = TRUE)`, which returns  
  # the log predictive density evaluations. This function then averages these log density evaluations and negates 
  # the result so that it can be interpreted as a measure of error, rather than a score. 
  #
  # Args:
  #    gp_mean: numeric, vector of GP mean predictions. 
  #    gp_var: numeric, vector of GP variance predictions.
  #    y_true: numeric, vector of true response values. Must be equal length as `gp_mean`.
  #    transformation_method: character(1), string specifying a transformation method; currently supports
  #                           "LNP", "rectified" and "truncated". `NA_character_` will apply 
  #                           no transformation. This is required as the log predictive density depends on the 
  #                           predictive distribution, which may not always be Gaussian. 
  #  
  # Returns:
  #    numeric, the pointwise negative log predictive density error metric.
  
  lpd <- get_GP_pointwise_predictive_density(y_vals = y_true, 
                                             gp_mean = gp_mean, 
                                             gp_var = gp_var, 
                                             transformation_method = transformation_method, 
                                             log = TRUE)
  
  return(-mean(lpd))
  
}


get_gp_mah <- function(gp_mean, y_true, gp_cov = NULL, gp_L = NULL, eps = sqrt(.Machine$double.eps), ...) {
  # Computes the Mahalanobis distance between GP mean predictions and a true baseline response.
  #
  # Args:
  #    gp_mean: numeric, vector of GP mean predictions. 
  #    y_true: numeric, vector of true response values. Must be equal length as `gp_mean`.
  #    gp_cov: matrix, M x M GP predictive covariance matrix evaluated at the M test points. If NULL, 
  #            `gp_L` must be provided. 
  #    gp_L: matrix, M x M lower triangular matrix, the lower Cholesky factor of `gp_cov`. If NULL, 
  #          `gp_cov` must be provided. 
  #    eps: numeric(1), value to add to diagonal of predictive covariance matrix in order to ensure 
  #         positive definiteness. 
  #
  # Returns:
  #    numeric, the Mahalanobis distance. 
  
  if(is.null(gp_cov) && is.null(gp_L)) stop("Either predictive covariance or Cholesky factor must be provided.")
  
  if(is.null(gp_L)) {
    gp_L <- t(chol(gp_cov + diag(eps, nrow = nrow(gp_cov))))
  }
  
  return(sqrt(sum(forwardsolve(gp_L, y_true - gp_mean)^2)))
  
}


get_gp_crps <- function(gp_mean, gp_var, y_true, transformation_method = NA_character_, ...) {
  # Computes the continuous ranked probability score (CRPS) of a GP forecast based on 
  # on the observed data `y_true`. `gp_mean` and `gp_var` are assumed to be vectors of 
  # GP predictive means and variances at the test locations for the underlying/untransformed 
  # GP. The predictive distribution can be transformed into something non-Gaussian by specifying 
  # `transformation_method`. The current allowed transformations "LNP" (log-normal process)
  # "truncated" (zero-truncated Gaussian) and "rectified" (zero censored/rectified Gaussian)
  # all admit closed-form expressions for the CRPS. The CRPS is computed at each test location 
  # and the average over all test locations is returned. 
  #
  # Args:
  #    gp_mean: numeric, vector of GP mean predictions. 
  #    gp_var: numeric, vector of GP variance predictions.
  #    y_true: numeric, vector of true response values. Must be equal length as `gp_mean`.
  #    transformation_method: character(1), string specifying a transformation method; currently supports
  #                           "LNP", "rectified" and "truncated". `NA_character_` will apply 
  #                           no transformation. This is required as the log predictive density depends on the 
  #                           predictive distribution, which may not always be Gaussian.
  #
  # Returns:
  #    numeric, the average CRPS. 
  
  gp_sd <- sqrt(gp_var)
  
  if(is.na(transformation_method)) {
    crps_vals <- crps_norm(y_true, mean = gp_mean, sd = gp_sd)
  } else if(transformation_method == "LNP") {
    crps_vals <- crps_lnorm(y_true, meanlog = gp_mean, sdlog = gp_sd)
  } else if(transformation_method == "truncated") {
    crps_vals <- crps_tnorm(y_true, location = gp_mean, scale = gp_sd, lower = 0)
  } else if(transformation_method == "rectified") {
    crps_vals <- crps_cnorm(y_true, location = gp_mean, scale = gp_sd, lower = 0)
  } else {
    stop("Invalid transformation method: ", transformation_method)
  }
  
  return(mean(crps_vals))
  
}


get_independent_gp_metric <- function(gp_mean_df, Y_true, metric, gp_var_df = NULL, transformation_method = NA_character_, 
                                      gp_cov = NULL, gp_L = NULL) {
  # A convenience function to compute a GP metric for each output variable and return the 
  # results. The GPs modeling each output are assumed independent. The names of the metric 
  # functions are standardized as "get_gp_<metric>". For metrics that require the GP predictive covariance, 
  # it is assumed that `gp_var_df` contains the variances that are desired to be used as the diagonals 
  # of the respective covariance matrices (i.e. these may include or not include the nugget). Thus, 
  # the diagonal of the covariance matrix `gp_cov[[j]]` is replaced with the column `gp_var_df[,j]`. 
  # Note that `predict_GP()` never includes the nugget in the diagonal of the predictive covariance, hence 
  # the required adjustment here if necessary. 
  #
  # Args:
  #    gp_mean_df: data.frame or matrix, of dimension M x P, where M is the number of test inputs 
  #                and P the number of output variables. The GP predictive means at the test points. 
  #    Y_true: matrix, of dimension M x P containing the true responses at the test inputs. 
  #    gp_var_df: data.frame or matrix, of dimension M x P. The GP predictive variances at the test points.
  #               Not required for every metric, so the default is NULL. 
  #    transformation_method: character(1), string specifying a transformation method; currently supports
  #                           "LNP", "rectified" and "truncated". The default, `NA_character_` will apply 
  #                           no transformation. This argument is only required by some metrics, such as 
  #                           "nlpd_pointwise" and "nlpd". 
  #    gp_cov: list of matrices of length P. Each matrix is an M x M GP predictive covariance matrix evaluated  
  #            at the M test points. Only required  for some metrics, e.g. "mah". 
  #    gp_L: list of matrices matries of length P. Each matrix is an M x M lower triangular matrix, 
  #          the lower Cholesky factor of `gp_cov`. Either `gp_cov` or `gp_L` must be provided for
  #          metrics that utilize GP predictive covariance. 
  #
  # Returns: 
  #    numeric, vector of length P containing the prediction metric for each output/GP. 

  metric_func <- get(paste0("get_gp_", metric))
  
  apply_func <- function(j) {
    cov_mat <- gp_cov[[j]]
    if(!is.null(cov_mat)) diag(cov_mat) <- gp_var_df[,j]
    
    metric_func(gp_mean = gp_mean_df[,j], 
                gp_var = gp_var_df[, j], 
                y_true = Y_true[,j], 
                transformation_method = transformation_method, 
                gp_cov = cov_mat, 
                gp_L = gp_L)
  }
  
  metrics <- sapply(seq(1, ncol(gp_mean_df)), apply_func)
  
  return(metrics)
  
}


get_emulator_comparison_metrics_validation_list <- function(emulator_info_list, X_test_list, Y_test_list, metrics, scale_inputs, 
                                                            output_variables = NULL, include_nug = TRUE) {
  # A wrapper for `get_emulator_comparison_metrics()` that loops over a list of different validation/test data and combines 
  # all of the results in a single data.table. 
  
  f_apply <- function(idx) {
    metrics_dt <- get_emulator_comparison_metrics(emulator_info_list = emulator_info_list, 
                                                  X_test = X_test_list[[idx]], 
                                                  Y_test = Y_test_list[[idx]], 
                                                  metrics = metrics, 
                                                  scale_inputs = TRUE, 
                                                  output_variables = output_variables, 
                                                  include_nug = include_nug)$metrics
    metrics_dt[, validation_data := idx]
  }
  
  metrics_list <- lapply(seq_along(X_test_list), f_apply)
                                                                               
  return(rbindlist(metrics_list, use.names = TRUE))
  
}


get_emulator_comparison_metrics <- function(emulator_info_list, X_test, Y_test, metrics, scale_inputs = FALSE, output_variables = NULL,
                                            emulator_pred_list = NULL, include_nug = TRUE) {
  # This function computes error metrics (or scores) for a set of different GPs, and returns the results in a data.frame 
  # which can be used to compare the performance across the different GPs. The response data `Y_test` is assumed to be unnormalized. 
  # The test inputs `X_test` may or may not be scaled, which can be accounted for by the argument `scale_inputs`. Metrics are 
  # computed on the original (unnormalized) scale. 
  #
  # Args:
  #    emulator_info_list: Named list, where each element is an "emulator info" object. The emulator info objects 
  #                        are themselves lists with elements "gp_fits", "input_bounds", "output_stats", and "settings". 
  #    X_test: matrix, of dimension M x D, where M is the number of test inputs and D the dimension of the parameter space. 
  #            The test inputs at which to predict in order to compute prediction metrics. 
  #    Y_test: matrix, of dimension M x P, the baseline true responses corresponding to the inputs `X_test`. The pth column 
  #            contains the responses for the pth output variable. Must be on the original (unnormalize scale). 
  #    metrics: character(), vector of prediction metrics to include. Currently supports: "rmse" (root mean squared error), 
  #             "srmse" (Standardized RMSE; normalized by GP predictive variances), "crps" (continuous rank probability score), 
  #             "nlpd_pointwise" (negative log predictive density, summed over test points), 
  #             "nlpd" (negative log predictive density, taking into account GP predictive covariance), "mah" (Mahalanobis distance). 
  #             "mah" and "nlpd" require computation of the GP predictive covariance at the test points, while the other metrics do not. 
  #    scale_inputs: If TRUE, uses the "input_bounds" element of the emulator info lists to scale `X_test` prior to prediction. 
  #                  If FALSE, assumes `X_test` is already properly scaled. 
  #    output_variables: character, vector of output variable names, which is passed to `predict_independent_GPs()`. If NULL, 
  #                      `predict_independent_GPs()` just assigns numbers to the outputs. 
  #    emulator_pred_list: list, the return value to `get_emulator_comparison_pred_df()` if this has already available. Otherwise this 
  #                        function is run to compute the predictions. 
  #    transformation_method: character(1), string specifying a transformation method; currently supports
  #                           "LNP", "rectified" and "truncated". The default, `NA_character_` will apply 
  #                           no transformation. This argument is only required by some metrics, such as 
  #                           "nlpd_pointwise" and "nlpd". 
  #    include_nug: logical(1), if TRUE the predictive GP variance is the predictive variance of the latent function plug the nugget. 
  #                 Otherwise, the nugget variance is not added. 
  #
  # Returns:
  #    list, with elements "pred" (containing the data.frame of GP predictions) and "metrics" (containing the data.frame of GP metrics). 
  
  # TODO: ensure Y_test has named columns.
  # TODO: enusre output_variables is non-NULL before passing to `get_emulator_comparison_pred_df()`. 
  # TODO: for now not including any metrics that require predictive covariance; need to modify a lot of functions to allow this.
  # TODO: Note that when populating the metrics data.frame, the order of the output variables is crucial; should probably make this more 
  #       robust to avoid mistakes with this. 

  # Compute emulator predictions. Note that predictions are untransformed; i.e. these are GP predictions, not log-normal process
  # or any other transformed GP. The transformations are taken into account when the metrics are computed, if necessary. 
  if(is.null(emulator_pred_list)) {
    include_cov_mat <- any(metrics %in% c("mah", "nlpd"))
    emulator_pred_list <- get_emulator_comparison_pred_df(emulator_info_list = emulator_info_list, 
                                                          X_test = X_test, 
                                                          scale_inputs = scale_inputs, 
                                                          output_variables = output_variables, 
                                                          include_cov_mat = include_cov_mat, 
                                                          denormalize_predictions = TRUE, 
                                                          transform_predictions = TRUE)
  }
  emulator_pred_df <- emulator_pred_list$df
  
  # data.frame to store metrics. 
  test_labels <- names(emulator_info_list)
  df_pred_cols <- colnames(emulator_pred_df)
  Y_test <- Y_test[, output_variables, drop = FALSE]
  metrics_df <- as.data.table(expand.grid(test_labels, output_variables, metrics))
  colnames(metrics_df) <- c("test_label", "output_variable", "metric_name")
  metrics_df[, metric_value := NA_real_]
  metrics_df[, c("design_rep", "design") := NA_character_]
  
  # Compute metrics for each emulator.
  for(j in seq_along(emulator_info_list)) {
    
    # Extract columns pertaining to GP predictive means, variances, etc. 
    test_lab <- test_labels[j]
    test_lab_split <- strsplit(test_lab, split = "_")[[1]]
    metrics_df[test_label == test_lab, design := substr(test_lab_split[1], start = nchar(test_lab_split[1]), stop = nchar(test_lab_split[1]))]
    metrics_df[test_label == test_lab, design_rep := substr(test_lab_split[2], start = nchar(test_lab_split[2]), stop = nchar(test_lab_split[2]))]
    gp_mean_cols <- paste(test_lab, "mean", output_variables, sep = "_")
    gp_var_cols <- paste(test_lab, "var", output_variables, sep = "_")
    gp_var_comb_cols <- paste(test_lab, "var_comb", output_variables, sep = "_")

    if(include_nug) {
      gp_var_df <- emulator_pred_df[, gp_var_comb_cols, drop = FALSE]
    } else {
      gp_var_df <- emulator_pred_df[, gp_var_cols, drop = FALSE]
    }
    
    for(metric in metrics) {
      metric_vals <- get_independent_gp_metric(gp_mean_df = emulator_pred_df[, gp_mean_cols, drop = FALSE], 
                                               gp_var_df = gp_var_df,
                                               Y_true = Y_test, 
                                               metric = metric, 
                                               transformation_method = emulator_info_list[[j]]$settings$transformation_method, 
                                               gp_cov = emulator_pred_list$cov_list[[j]])
      metrics_df[(test_label == test_lab) & (metric_name == metric), metric_value := metric_vals]
    }
  }
  
  return(list(pred = emulator_pred_df, metrics = metrics_df))
  
}


compute_gp_metrics <- function(gp_mean_pred, output_true, metrics, gp_var_pred = NULL, transformation_method = NA_character_, 
                               gp_cov_pred = NULL, gp_L_pred = NULL) {
  # A convenience function to compute a set of (univariate) GP metrics, given GP predictions and the true baseline output values. 
  # The names of the metric functions are standardized as "get_gp_<metric>".
  #
  # Args:
  #    gp_mean_pred: numeric vector of length M. The GP mean predictions at the M test locations. 
  #    output_true: numeric vector of length M, the true response/output values at the M test locations. 
  #    metrics: character, vector of metric names; used to call function "get_gp_<metric>". 
  #    gp_var_pred: numeric vector of length M. The GP variance predictions at the M test locations. Not required for 
  #                 the computation of all metrics, so the default is NULL. 
  #    transformation_method: character(1), string specifying a transformation method; currently supports
  #                           "LNP", "rectified" and "truncated". The default, `NA_character_` will apply 
  #                           no transformation. This argument is only required by some metrics, such as 
  #                           "nlpd_pointwise" and "nlpd". 
  #    gp_cov_pred: matrix of dimension M x M, the GP predictive covariance matrix evaluated  
  #                 at the M test points. Only required  for some metrics, e.g. "mah". 
  #    gp_L_pred: matrix of dimension M x M, lower triangular matrix, the lower Cholesky factor
  #               of `gp_cov_pred`. Either `gp_cov_pred` or `gp_L_pred` must be provided for
  #               metrics that utilize GP predictive covariance. 
  #
  # Returns: 
  #    numeric, vector of length `length(metrics)`. The vector of computed metrics; the names attribute of the vector is 
  #    set to `metrics`. 
  
  apply_func <- function(metric) {
    
    metric_func <- get(paste0("get_gp_", metric))
    
    metric_func(gp_mean = gp_mean_pred, 
                gp_var = gp_var_pred, 
                y_true = output_true, 
                transformation_method = transformation_method, 
                gp_cov = gp_cov_pred, 
                gp_L = gp_L_pred)
  }
  
  metric_results <- sapply(metrics, apply_func)
  names(metric_results) <- metrics
  
  return(metric_results)
  
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
  # will be set to a log base 10 scale. 
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

# TODO: comment. 
plot_gp_fit_2d <- function(X_test, X_train = NULL, y_test = NULL, gp_mean_pred = NULL, gp_var_pred = NULL, post_samples = NULL,
                           true_theta = NULL, xlab = "", ylab = "", main_title = "", 
                           transformation_method = NA_character_, log_predictive_density = NULL, raster = FALSE) {
  # Set raster=TRUE when X_test is a grid of points (equally spaced). If not equally, spaced set to FALSE. 
  
  # Log predictive density evaluations at test points. 
  df_test <- as.data.frame(X_test)
  colnames(df_test) <- c("theta1", "theta2")
  
  if(is.null(log_predictive_density)) {
    log_predictive_density <- get_GP_pointwise_predictive_density(y_test, gp_mean_pred, gp_var_pred, 
                                                                  transformation_method = transformation_method, log = TRUE)
  }
  df_test <- cbind(df_test, density_vals = log_predictive_density)
  if(raster) {
    plt <- ggplot() + geom_tile(data = df_test, aes(x = theta1, y = theta2, fill = density_vals))
  } else {
    plt <- ggplot() + geom_point(data = df_test, aes(x = theta1, y = theta2, color = density_vals))
  }
  
  # Design points. 
  if(!is.null(X_train)) {
    df_train <- as.data.frame(X_train)
    colnames(df_train) <- c("theta1", "theta2")
    plt <- plt + geom_point(data = df_train, aes(x = theta1, y = theta2), color = "red")
  }
  
  # Contours of true posterior. 
  if(!is.null(post_samples)) {
    df <- as.data.frame(post_samples)
    colnames(df) <- c("theta1", "theta2")
    
    plt <- plt + geom_density_2d(data = df, aes(x = theta1, y = theta2))
  }
  
  # Marker for the true parameter values. 
  if(!is.null(true_theta)) {
    plt <- plt + geom_point(data = data.frame(theta1 = true_theta[1], theta2 = true_theta[2]), 
                            aes(x = theta1, y = theta2), color = "red", shape = 17)
  }
  
  # Labels and title. 
  plt <- plt + 
         xlab(xlab) + 
         ylab(ylab) + 
         ggtitle(main_title)
  
  
  return(plt)  
  
}


plot_emulator_comparison_violin_plots <- function(emulator_metrics_df, fill_col, metrics, output_variables, include_points = TRUE, nrow = 1) {
  # `fill_col` will typically either be "design" (when comparison across different designs is desired) or "emulator_setting" (if comparison 
  # across different emulator specifications is desired). 
  
  plts <- vector(mode = "list", length = length(metrics) * length(output_variables))
  N_metrics <- length(metrics)
  
  for(i in seq_along(output_variables)) {
    for(j in seq_along(metrics)) {
      idx <- N_metrics * (i - 1) + j
      plts[[idx]] <- ggplot(data = emulator_metrics_df[(metric_name == metrics[j]) & (output_variable == output_variables[i]),], 
                            aes(x = .data[[fill_col]], y = metric_value, fill = .data[[fill_col]])) + 
                      geom_violin() +
                      ylab(paste(metrics[j], output_variables[i], sep = ", "))
      if(include_points) plts[[idx]] <- plts[[idx]] + geom_jitter(color = "black", size = 0.8, alpha = 0.9, width = 0.2)
    }
  }
  
  return(arrangeGrob(grobs = plts, nrow = nrow))
  
}


get_GP_pointwise_predictive_density <- function(y_vals, gp_mean, gp_var, transformation_method = NA_character_, log = FALSE) {
  # Computes the pointwise GP predictive density at a set of input locations. The density evaluations
  # to not take into account GP predictive covariance, hence the "pointwise". The arguments 
  # `gp_mean` and `gp_var` must pertain to a GP, but the argument `transformation_method` can be used to 
  # transform to other predictive distributions. In this case, `y_vals` is assumed to already be on the 
  # desired scale (`y_vals` is never transformed). For example, the response variable y may have been 
  # modeled as a log-normal process (LNP) in which case log(y) was modeled as a GP. `gp_mean` and 
  # `gp_var` are the predictive moments for this GP. Setting `transformation_method` to "LNP" will 
  # then correctly treat the predictive density as a log-normal density. Note that the vectors 
  # `y_vals`, `gp_mean`, and `gp_var` are all assumed to be of equal length, with each respective entry 
  # corresponding to the same input point. 
  #
  # Args:
  #    y_vals: numeric(), vector of response values at which to evaluate the predictive density. 
  #    gp_mean: numeric(), vector of predictive mean evaluations of the underlying (untransformed) GP. 
  #    gp_var: numeric(), vector of predictive variance evaluations of the underlying (untransformed) GP.
  #    transformation_method: character(1), string specifying a transformation method; currently supports
  #                           "LNP", "rectified" and "truncated". The default, `NA_character_` will apply 
  #                           no transformation. 
  #    log: logical(1), if TRUE returns log density.
  #
  # Returns:
  #    numeric(), vector of (potentially transformed) GP pointwise density evaluations, of equal length 
  #    as `y_vals`. 

  gp_sd <- sqrt(gp_var)
  
  if(is.na(transformation_method)) {
    predictive_density <- dnorm(y_vals, gp_mean, gp_sd, log = log)
  } else if(transformation_method == "LNP") {
    predictive_density <- dlnorm(y_vals, gp_mean, gp_sd, log = log)
  } else if(transformation_method == "truncated") {
    predictive_density <- dtruncnorm(y_vals, a = 0, b = Inf, mean = gp_mean, sd = gp_sd)
    if(isTRUE(log)) predictive_density <- log(predictive_density)
  } else if(transformation_method == "rectified") {
    predictive_density <- density_rectified_norm(y_vals, mean_norm = gp_mean, sd_norm = gp_sd, log = log)
  } else {
    stop("Invalid transformation method: ", transformation_method)
  }
  
  return(predictive_density)
  
}


# TODO: allow SSR_true to be passed in to avoid re-computing every time. 
get_GP_SSR_log_pred_density <- function(computer_model_data, theta_vals, emulator_info) {
  # Computes the GP log predictive density for an independent GP model at a set of input points. 
  # In particular, computes the predictive distribution for each output individually and evaluates
  # the log predictive density for each output at the true SSR values. The overall log predictive 
  # density is then the sum of the log predictive densities for each output. 
  #
  # Args:
  #    computer_model_data: list, the standard computer model data list. 
  #    theta_vals: matrix, of shape M x d, where M is the number of input points and d is the dimension of 
  #                the parameter space. Each row is an input points. 
  #    emulator_info: list, the emulator info list, as passed into `mcmc_calibrate_ind_GP`. Note that the element 
  #                   "output_stats" of `emulator_info` is assumed to already be transformed correctly; for example, 
  #                   for a log-normal process, "output_stats" should be the log-transformed version. 
  #
  # Returns:
  #    numeric vector, of length equal to the number of rows in `theta_vals` containing the log predictive 
  #    density evaluations. 
  
  # Get true SSR values. 
  SSR_true <- get_computer_model_SSR(computer_model_data, theta_vals = theta_vals, na.rm = TRUE)
  
  # GP predictive mean and variance of the SSR values. 
  thetas_scaled <- scale_input_data(theta_vals, input_bounds = emulator_info$input_bounds)
  gp_pred_list <- predict_independent_GPs(X_pred = thetas_scaled, gp_obj_list = emulator_info$gp_fits, 
                                          gp_lib = emulator_info$settings$gp_lib, include_cov_mat = FALSE, denormalize_predictions = TRUE,
                                          output_stats = emulator_info$output_stats)
  
  # Evaluate GP predictive density at the true SSR values. 
  log_pred_density <- vector(mode = "numeric", length = nrow(theta_vals))
  
  for(j in seq(1, ncol(SSR_true))) {
    log_pred_density <- log_pred_density + get_GP_pointwise_predictive_density(SSR_true[,j], gp_pred_list[[j]]$mean, gp_pred_list[[j]]$var, 
                                                                               transformation_method = emulator_info$settings$transformation_method, log = TRUE)
  }
  
  return(log_pred_density)
  
}


estimate_integrated_neg_log_pred_density_err <- function(computer_model_data, theta_samples_post, emulator_info) {
  # Computes a Monte Carlo estimate of error for a GP-approximated posterior distribution, with respect to a true 
  # known posterior distribution. The error is defined as the expectation of the GP negative log predictive density
  # with respect to the true posterior. This is estimated with a Monte Carlo estimate based on samples `theta_samples_post`
  # from the true posterior. Note that the log predictive density is multiplied by -1 so that it can be interpreted as 
  # a measure of error, with smaller being better. 
  #
  # Args:
  #    computer_model_data: list, the standard computer model data list. 
  #    theta_samples_post: matrix, of shape M x d, where M is the number of input points sampled from the true posterior 
  #                        and d is the dimension of the parameter space. Each row is an input point sampled from the true posterior.
  #    emulator_info: list, the emulator info list, as passed into `mcmc_calibrate_ind_GP`. 
  #
  # Returns:
  #    list, the first element is the vector of negative log predictive density evaluations. The second is the estimate of the 
  #    error measure, which is just the average of the elements in the vector of the first element. 
  
  # Compute log predictive density at input points sampled from true posterior. 
  log_pred_density <- get_GP_SSR_log_pred_density(computer_model_data, theta_samples_post, emulator_info)
  
  # Compute sample average of negative log predictive density, which estimates the integral of the negative log predictive 
  # density with respect to the true posterior. 
  err_estimate <- -mean(log_pred_density)
  
  return(list(neg_log_pred_density = -log_pred_density, 
              err_estimate = err_estimate))
  
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


density_rectified_norm <- function(x, mean_norm = 0, sd_norm = 0, allow_inf = FALSE, log = FALSE) {
  # Computes the density of a rectified Gaussian distribution. The arguments 
  # `mean_norm` and `sd_norm` are the mean and standard deviation of the Gaussian 
  # that gives rise to the rectified Gaussian, not the moments of the rectified 
  # Gaussian itself. The rectified Gaussian density is infinite at the value 0, 
  # so the argument `allow_inf` allows the user to specify whether this should 
  # be considered an error or not. 
  #
  # Args:
  #    x: numeric(), vector of points at which to evaluate the density. 
  #    mean_norm: numeric(), vector of Gaussian means. 
  #    sd_norm: numeric(), vector of Gaussian standard deviations. 
  #    allow_inf: logical(1), if TRUE returns `inf` for values at `x` which are 0. 
  #               Otherwise, 0 values of `x` will invoke an error. Default is FALSE. 
  #    log: logical(1), if TRUE returns log density. 
  #
  # Returns:
  #    numeric() vector of length equal to length of `x` containing the density evaluations. 
  
  if(!allow_inf && any(x == 0)) {
    stop("Density of rectified Gaussian at x = 0 is infinite.")
  }
  
  x_densities <- vector(mode = "numeric", length = length(x))
  
  x_densities[x == 0] <- Inf
  x_densities[x != 0] <- dnorm(x, mean_norm, sd_norm) * as.numeric(x > 0)
  
  if(isTRUE(log)) return(log(x_densities))
  return(x_densities)
  
}


# ------------------------------------------------------------------------------
# Design Points
# ------------------------------------------------------------------------------

# TODO: update comments and think about a better way of handling log output stats. 
get_input_output_design <- function(N_points, computer_model_data, theta_prior_params, scale_inputs = TRUE, normalize_response = TRUE,
                                    param_ranges = NULL, output_stats = NULL, log_output_stats = NULL, transformation_method = NA_character_,
                                    design_method = "LHS", order_1d = TRUE, tail_prob_excluded = 0.01, na.rm = TRUE, 
                                    design_candidates = NULL, design_candidates_weights = NULL, design_seed = NULL) {
  # Generates input points in parameter space and runs the VSEM model at these points to obtain the corresponding outputs, which is the L2 error between 
  # the model outputs and observed data. Handles scaling of input data and normalization of response data. Also handles log-transformation of response data
  # in the case of the log-normal process. 
  #
  # Args:
  #    N_points: integer(1), the number of input points to generate. 
  #    computer_model_data: list, the standard computer model data list. 
  #    theta_prior_params: data.frame containing the prior distribution information of the input 
  #                        parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                        for the requirements of this data.frame. 
  #    scale_inputs: logical(1), if TRUE will linearly scale the samples to the unit hypercube. Will return the scaled 
  #             sample in addition to the un-scaled one, and will also return the information used for the scaling 
  #             so that the transformation can be inverted. 
  #    normalize_response: logical(1), if TRUE will compute a Z-score transformation of the response/output variable. Will return the normalized 
  #                        variable in addition to the unnormalized one and will return the information used for the normalization so that the 
  #                        transformation can be inverted.
  #    param_ranges: matrix, of shape 2 x d, where d is the dimension of the parameter space. The columns correspond 
  #                  to the respective rows in `theta_prior_params`. The first and second rows are the lower and
  #                  upper bounds on the sample ranges of each parameter, respectively. Default is NULL, in which case 
  #                  no additional bounds are imposed on the samples, other than those already provided by the prior 
  #                  distributions.
  #    output_stats: matrix of means/variances used for normalizing response, as returned by `prep_GP_training_data()`. If non-NULL and 
  #                  `normalize_response` if TRUE, then uses this information for the normalization. If NULL, then will instead compute 
  #                   the means/variances of the response data and use those computations for the normalization. 
  #    log_output_stats: the log-transformed analog of `output_stats`. This is only used when `transformation_method` is "LNP".  
  #    transformation_method: character(1), if `transformation_method` is "LNP", then will include the log-transformed response and output 
  #                           statistics in addition to the non-log transformed versions. Transformations "rectified" and "truncated" have 
  #                           no effect. 
  #    design_method: character(1), the algorithm used to generate the inputs. Currently supports "LHS" or "grid". 
  #    order_1d: logical(1), only relevant if the dimension of the input space (i.e. the number of 
  #              calibration parameters) is one-dimensional. In this case, if `order_1d` is TRUE then 
  #              the design points will be sorted in increasing order, with the training response 
  #              values ordered accordingly. This is convenient when plotting 1d GP plots. 
  #    tail_prob_excluded: numeric(), this is only relevant in certain cases, such as a Gaussian prior with "grid" design method. In this case 
  #                        the Gaussian has infinite support, but the grid method requires bounded support. If `tail_prob_excluded` is 0.01 then 
  #                        these bounds will be set to the .5% and 99.5% quantiles of the Gaussian. 
  #    na.rm: logical(1), whether or not to remove NA values from the sum of squares calculation, when computing L2 error between computer model 
  #           outputs and observed data. Default is FALSE. 
  #    design_candidates, design_candidates_weights: Used by `get_sample_candidates_design()`, see explanation in `get_input_design()`. 
  #    design_seed: integer, an optional random seed to set before generating the design. 
  #           
  # Returns:
  #    list, containing all elements returned by `get_input_design()`. In addition, will at least contain element "outputs", containing an N x p matrix 
  #    storing the squared L2 errors (N = number inputs, p = number output variables). Will also optionally contain elements "outputs_normalized", "output_stats", 
  #    "log_outputs", "log_outputs_normalized", and "log_output_stats". 
  
  if(!is.null(design_seed)) set.seed(design_seed)
  
  # Input points. 
  design_list <- get_input_design(N_points = N_points, 
                                  theta_prior_params = theta_prior_params, 
                                  design_method = design_method, 
                                  scale_inputs = scale_inputs, 
                                  param_ranges = param_ranges, 
                                  order_1d = order_1d, 
                                  tail_prob_excluded = tail_prob_excluded, 
                                  design_candidates = design_candidates, 
                                  design_candidates_weights = design_candidates_weights)

  # Run model at inputs to produce outputs. 
  design_list[["outputs"]] <- get_computer_model_SSR(computer_model_data, theta_vals = design_list$inputs, na.rm = na.rm)
  
  # Normalize outputs. 
  if(normalize_response) {
    if(is.null(output_stats)) {
      design_list[c("outputs_normalized", "output_stats")] <- prep_GP_training_data(Y = design_list$outputs, normalize_Y = TRUE)[c("Y", "output_stats")]
    } else {
      design_list[["outputs_normalized"]] <- normalize_output_data(Y = design_list$outputs, output_stats) 
    }
  }
  
  # Include log-transformed data for log-normal process. 
  if(isTRUE(transformation_method == "LNP")) {
    design_list[["log_outputs"]] <- log(design_list$outputs)
    
    if(is.null(log_output_stats)) {
      design_list[c("log_outputs_normalized", "log_output_stats")] <- prep_GP_training_data(Y = design_list$log_outputs, normalize_Y = TRUE)[c("Y", "output_stats")]
    } else {
      design_list[["log_outputs_normalized"]] <- normalize_output_data(Y = design_list$log_outputs, log_output_stats)
    }
    
  }
  
  if(!is.null(output_stats) && isTRUE(transformation_method == "LNP") && is.null(log_output_stats)) {
    message("output_stats is non-NULL but log_output_stats is NULL; may be a mistake.")
  }
  
  return(design_list)
  
}


get_input_design <- function(N_points, theta_prior_params, design_method, scale_inputs, param_ranges = NULL, 
                             order_1d = FALSE, tail_prob_excluded = 0.01, design_candidates = NULL, 
                             design_candidates_weights = NULL) {
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
  #    design_method: character(1), the algorithm used to generate the inputs. Currently supports "LHS", "grid", or "sample_candidates".  
  #    scale_inputs: logical(1), if TRUE will linearly scale the samples to the unit hypercube. Will return the scaled 
  #             sample in addition to the un-scaled one, and will also return the information used for the scaling 
  #             so that the transformation can be inverted. 
  #    param_ranges: matrix, of shape 2 x d, where d is the dimension of the parameter space. The columns correspond 
  #                  to the respective rows in `theta_prior_params`. The first and second rows are the lower and
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
  #    design_candidates: matrix, with d columns. Only relevant for design method "sample_candidates". 
  #                       A set of input points (with each row being a point) from which the design will be sampled. The points will 
  #                       be sampled without replacement. By default the sampling will be uniform, but `design_candidates_weights` can 
  #                       be specified for non-uniform sampling. 
  #    design_candidates_weights: numeric, of length equal to the number of rows in `design_candidates`. Only relevant for design method "sample_candidates".
  #                               The values must be non-negative, and will be normalized and used as weights for the respective points in `design_candidates`
  #                               when sampling to select the design points. These weights are passed to the base R function `sample(..., replace = FALSE)`
  #                               so see documentation of this function for more details. 
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
  } else if(design_method == "sample_candidates") {
    X_design <- get_sample_candidates_design(N_points, candidates = design_candidates, 
                                             weights = design_candidates_weights, param_ranges = param_ranges)
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
  #    param_ranges: matrix, of shape 2 x d, where d is the dimension of the parameter space. The columns correspond 
  #                  to the respective rows in `theta_prior_params`. The first and second rows are the lower and
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
  d <- nrow(theta_prior_params)
  
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
  
  colnames(X_LHS) <- rownames(theta_prior_params)
  
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
  #    param_ranges: matrix, of shape 2 x d, where d is the dimension of the parameter space. The columns correspond 
  #                  to the respective rows in `theta_prior_params`. The first and second rows are the lower and
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
  d <- nrow(theta_prior_params)
  
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


get_sample_candidates_design <- function(N_points, candidates, weights = NULL, param_ranges = NULL) {
  # Samples a subset of the collection of points `candidates` without replacement, potentially 
  # weighted by `weights` and ensuring bound constraints `param_ranges`. 
  #
  # Args:
  #    N_points: integer(1), the number of points to sample. Must be less than the number of rows in `candidates`. 
  #    candidates: matrix, each row is considered a point. This is the collection of candidate points from which the 
  #                sampling is done. 
  #    weights: numeric, of length equal to the number of rows of `candidates`. Non-negative weights used for the 
  #             sampling, passed to the `prob` argument of the base R function `sample()`. If NULL, uses uniform weights.
  #    param_ranges: matrix, of shape 2 x d, where d is the dimension of the parameter space. The columns correspond 
  #                  to the respective rows in `theta_prior_params`. The first and second rows are the lower and
  #                  upper bounds on the sample ranges of each parameter, respectively. Default is NULL, in which case 
  #                  no additional bounds are imposed on the samples, other than those already provided by the prior 
  #                  distributions.
  #
  # Returns:
  #    matrix, of dimension N_points x d; a matrix consisting of the sub-collection of rows of `candidates`. 

  N_candidates <- nrow(candidates)
                 
  # If `param_ranges` is specified, set weights to zero for candidates outside of the parameter range. 
  if(!is.null(param_ranges)) {
    if(is.null(weights)) weights <- rep(1, N_candidates)
    
    outside_range <- rep(FALSE, N_candidates)
    for(d in 1:ncol(candidates)) {
      outside_range <- outside_range | (candidates[,d] < param_ranges[1,d]) # Check lower bound. 
      outside_range <- outside_range | (candidates[,d] > param_ranges[2,d]) # Check upper bound.
    }
    
    weights[outside_range] <- 0
    
  }
  
  # Sample rows from `candidates`.                
  sample_idx <- sample(1:N_candidates, size = N_points, replace = FALSE, prob = weights)
  X_design <- candidates[sample_idx,,drop = FALSE]
  
  return(X_design)
  
}


get_design_list <- function(design_settings, computer_model_data, theta_prior_params, reps = 1, include_log = FALSE, ...) {
  # Creates a list of different designs (sets of emulator training points). Each row of 
  # `design_settings` specifies one specific design specification. By settings `reps` > 1, 
  # each design specification will be generated `reps` number of times. This is useful when 
  # comparing emulators, in order to average over the randomness in the design algorithms. 
  #
  # Args:
  #    design_settings: data.frame, with column names "N_design", "design_method", "scale_inputs", 
  #                     and "normalize_response". Currently supported design methods are 
  #                     "grid", "LHS", and "sample_candidates". "N_design" is the number of design points. 
  #                     The scaling and normalizing columns are logicals, indicating whether 
  #                     design points should be scaled to lie within a unit hypercube and 
  #                     whether the design outputs/response values should be undergo a Z-score
  #                     normalization. 
  #    computer_model_data: list, the standard computer model data list. 
  #    theta_prior_params: data.frame containing the prior distribution information of the input 
  #                        parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                        for the requirements of this data.frame. 
  #    reps: integer(1), the number of reps of each design specification to create. Default is 1. 
  #    include_log: logical(1), whether or not to include the log-transformed response in the designs. 
  #                 Only relevant for the log-normal process emulator.
  #    ...: other arguments passed to `get_input_output_design()`. For example, should include "design_candidates" 
  #         and "design_candidates_weights" if design method "sample_candidates" is being used. 
  #
  # Returns: 
  #    Returns a list of length `nrow(design_settings) * reps`. The names of the list are of the form 
  #    "design<i>_rep<j>" where <i> is the corresponding row number in `design_settings` and <j> indexes 
  #    the reps. 
  
  design_list <- list()
  transformation_method <- ifelse(include_log, "LNP", NA_character_)
  
  for(i in seq(1, nrow(design_settings))) {
    for(j in seq_len(reps)) {
      label <- paste0("design", i, "_rep", j)
      design_list[[label]] <- get_input_output_design(N_points = design_settings[i, "N_design"], 
                                                      computer_model_data = computer_model_data, 
                                                      theta_prior_params = theta_prior_params, 
                                                      scale_inputs = design_settings[i, "scale_inputs"], 
                                                      normalize_response = design_settings[i, "normalize_response"], 
                                                      transformation_method = transformation_method,
                                                      design_method = design_settings[i, "design_method"], 
                                                      order_1d = TRUE,
                                                      na.rm = TRUE, ...)
    }
  }
  
  return(design_list)
  
}


get_design_list_test_data <- function(N_test_points, N_test_sets, design_method, design_list, computer_model_data, 
                                      theta_prior_params, transformation_method = NA_character_, ...) {
  # Creates validation/test sets given a list of design//training sets. The main impetus for this function 
  # is that typically test points should only be generated within the bounds of the design set to avoid 
  # GP extrapolation. Therefore, when testing across different designs, the validation sets must satisfy 
  # the bounds for all designs. This function ensures all of these constraints are met so that each 
  # validation set can be used to test every design without extrapolation. Note that the design list 
  # is assumed to contain designs for a single emulator specification, meaning that if for example 
  # a log-normal process is utilized, then all validation sets will contain the log-transformed 
  # responses. It should also be noted that when transformed GP predictions back to the original, 
  # un-normalized scale, the transformation applied will be different for each design. Therefore, 
  # only the untransformed responses are returned here, rather than having to store a different 
  # transformed validation set for each design set. Similarly, the Therefore, the validation sets must be 
  # transformed when comparing to GP predictions that are on the original scale. Finally, note that 
  # the validation sets created here are based on the prior distributions on the calibration parameters. 
  # For certain validation exercises, it may be more desirable to test with respect to samples from the 
  # posterior distribution. 
  #
  # Args:
  #    N_test_points: integer(1), the number of validation points in the validation sets. 
  #    N_test_sets: integer(1), the number of replicate validation sets to create. Determines the length of the 
  #                 list returned by this function. 
  #    design_method: character(1), determines the method used to sample the validation points. Currently supports
  #                   "LHS", "grid", or "sample_candidates". 
  #    design_list: list, as returned by `get_design_list()`. Only utilized here to set the parameter ranges to 
  #                 be used when sampling the validation points. 
  #    computer_model_data: list, the computer model data list. 
  #    theta_prior_params: data.frame containing the prior distribution information of the input 
  #                        parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                        for the requirements of this data.frame. 
  #    transformation_method: character(1), if `transformation_method` is "LNP", then will include the log-transformed response and output 
  #                           statistics in addition to the non-log transformed versions. Transformations "rectified" and "truncated" have 
  #                           no effect. 
  #    ...: other arguments passed to `get_input_output_design()`. For example, should inlcude "design_candidates" 
  #         and "design_candidates_weights" if design method "sample_candidates" is being used.
  #
  # Returns:
  #    list, of length `N_test_sets`. Each element is itself a list returned by `get_input_output_design()`. 
  
  validation_list <- vector(mode = "list", length = N_test_sets)
  
  # Set the param_ranges so that validation points are within the ranges for all designs. 
  min_bounds_mat <- do.call("rbind", lapply(design_list, function(l) l$input_bounds[1,]))
  max_bounds_mat <- do.call("rbind", lapply(design_list, function(l) l$input_bounds[2,]))
  
  min_bounds <- apply(min_bounds_mat, 2, max)
  max_bounds <- apply(max_bounds_mat, 2, min)
  param_ranges <- rbind(min_bounds, max_bounds)

  for(j in seq_along(validation_list)) {
    validation_list[[j]] <- get_input_output_design(N_points = N_test_points, 
                                                    computer_model_data = computer_model_data, 
                                                    theta_prior_params = theta_prior_params, 
                                                    scale_inputs = FALSE, 
                                                    normalize_response = FALSE, 
                                                    param_ranges = param_ranges, 
                                                    transformation_method = transformation_method, 
                                                    design_method = design_method, 
                                                    order_1d = TRUE, na.rm = TRUE, ...)
    
  }
  
  return(validation_list)
  
}


# ------------------------------------------------------------------------------
# Sampling GPs
# ------------------------------------------------------------------------------


sample_GP_pointwise <- function(gp_means, gp_vars, transformation_method = NA_character_, idx_selector = NULL) {
  # Draws a single n-dimensional sample from a GP (or transformation of a GP) at a set of n input locations. 
  # This sample does not consider correlation across input points, hence the "pointwise" in the function name. 
  #
  # Args:
  #    gp_mean: numeric(), vector of means of Gaussian distributions (not the means of the transformed GPs). 
  #    gp_var: numeric(), vector of variances of the Gaussian distributions. 
  #    transformation_method: character(1), if "truncated", will convert the GP distribution to truncated Gaussian distribution. 
  #                           If "rectified", will instead transform to rectified Gaussian. If "LNP" will exponentiate the GP, resulting in a 
  #                           log-normal process. The default is not to perform any transformation. 
  #    idx_selector: integer(), vector of indices of which to select in `gp_means` and `gp_vars`. Default selects all indices. 
  #
  # Returns:
  #    numeric(), vector of length equal to length(gp_means) = length(gp_vars). The GP (or transformed GP) samples. 
  
  if(length(gp_means) != length(gp_vars)) {
    stop("gp_means and gp_vars must have equal length.")
  }
  
  if(!is.null(idx_selector)) {
    gp_means <- gp_means[idx_selector]
    gp_vars <- gp_vars[idx_selector]
  }
  n <- length(gp_means)
  
  if(is.na(transformation_method)) {
    sample <- gp_means + sqrt(gp_vars) * rnorm(n)
  } else if(transformation_method == "truncated") {
    sample <- rtruncnorm(1, a = 0, b = Inf, mean = gp_means, sd = sqrt(gp_vars))
  } else if(transformation_method == "rectified") {
    sample <- pmax(0, gp_means + sqrt(gp_vars) * rnorm(n))
  } else if(transformation_method == "LNP") {
    sample <- exp(gp_means + sqrt(gp_vars) * rnorm(n))
  } else {
    stop("Invalid transformation method: ", transformation_method)
  }
  
  if(any(is.na(sample)) || is.null(sample)) {
    stop("GP sample is NA or NULL.")
  }
  
  return(sample)
  
}


sample_independent_GPs_pointwise <- function(gp_pred_list, transformation_methods = NA_character_, idx_selector = NULL, include_nugget = TRUE) {
  # A generalization of `sample_GP_pointwise()` that allows samples to be drawn from multiple independent GPs. 
  # Allows for potentially different output transformations for each GP. 
  #
  # Args:
  #    gp_pred_list: list, with one element per GP. Each element must itself be a list with elements "mean" and "var"
  #                  storing numeric vectors for the GP means and variances, respectively. These are always the means and variances
  #                  of the GP, not the transformed GP. 
  #    transformation_methods: either a numeric vector of equal length as `gp_pred_list` storing the transformation to apply to each 
  #                            GP ("truncated", "rectified", "LNP", or NA_character_ for no transformation). If a single element, will 
  #                            use the same method for all GPs. 
  #    idx_selector: integer(), vector of indices of which to select in the mean/variance numeric vectors. Default selects all indices. Note that 
  #                  the same indices will be selected across all GPs.
  #    include_nugget: logical(1), if TRUE then the variance used in the sampling is the observation variance, which is the variance of the latent
  #                    function values plug the "nugget" term. Otherwise, the nugget will not be added to the variance. 
  # 
  # Returns: 
  #    matrix, of dimension N x length(gp_pred_list) where N is the number of input points. The matrix stores the samples for each GP in the 
  #    columns of the matrix. This will work even if the number of input points differs per GP, but note that this will result in NAs for the 
  #    GPs with fewer input points. 
  
  N_GPs <- length(gp_pred_list)
  if(length(transformation_methods) == 1) {
    transformation_methods <- rep(transformation_methods, N_GPs)
  } 
  
  if(!is.null(idx_selector)) {
    N_row <- length(idx_selector)
  } else {
    N_inputs <- sapply(gp_pred_list, function(x) length(x$mean))
    N_row <- max(N_inputs)
    if(!all(N_inputs == N_row)) {
      message("Returned matrix with samples will have NAs due to different numbers of input locations across GPs.")
    }
  }
  
  gp_samples <- matrix(nrow = N_row, ncol = N_GPs)
  for(j in seq_len(N_GPs)) {
    if(include_nugget) vars <- gp_pred_list[[j]]$var_comb
    else vars <- gp_pred_list[[j]]$var
    
    gp_samples[,j] <- sample_GP_pointwise(gp_pred_list[[j]]$mean, vars, transformation_methods[j], idx_selector)
  }
  
  return(gp_samples)
  
}


sample_GP_lpost_theta <- function(theta_vals_scaled = NULL, theta_vals_unscaled = NULL, emulator_info_list,  
                                  computer_model_data, theta_prior_params, sig2_eps, N_samples = 1, gp_pred_list = NULL, 
                                  include_sig_eps_prior = FALSE) {
  # In the loss emulation setting, where the squared error maps (SSR) are modeled as GPs, these GPs induce a 
  # random field approximation on the unnormalized log (conditional) posterior 
  # log pi(theta|Sigma) := log p(Y|theta, Sigma) + log pi_0(theta). This random field approximation is also 
  # a Gaussian process over the input space of calibration parameters `theta`. This function sampled from this 
  # random field representation of the unnormalized log (conditional) posterior at a discrete set of input 
  # locations `theta_vals_unscaled`; `N_samples` are drawn at each input location. Currently, the sampling does 
  # not take into account the covariance structure of the random field; the samples are drawn independently at each
  # input location. If `include_sig_eps_prior` is TRUE, then the samples are instead drawn from the joint unnormalized 
  # log posterior density log pi(theta, Sigma) := log p(Y|theta, Sigma) + log pi_0(theta) + log pi_0(Sigma). 
  #
  # Args:
  #    theta_vals_scaled: matrix of dimension M x D of input locations (scaled to lie in unit hypercube) at which to sample
  #                       the log density values. Each row is a location. 
  #    theta_vals_scaled: matrix, of dimension M x D; the same as `theta_vals_scaled` but the inputs are unscaled. Both the 
  #                       scaled and unscaled inputs are required here, but only one of them need be passed, as the scaling 
  #                       can be performed using the information in `emulator_info_list$input_bounds`. 
  #    emulator_info_list: list, the emulator info list. 
  #    computer_model_data: list, the computer model data list. 
  #    theta_prior_params: data.frame containing the prior distribution information of the input 
  #                        parameters, with each row corresponding to a parameter. See `calc_lprior_theta()`
  #                        for the requirements of this data.frame. 
  #    sig2_eps: numeric, vector of length P containing the observation variance parameters. 
  #    N_samples: integer, the number of samples to draw at each input location. 
  #    gp_pred_list: A list, as returned by `predict_independent_GPs()`. This allows the predictive means and 
  #                  variances of the underlying GP emulators evaluated at the input locations at the to be passed 
  #                  if they have already been computed. If NULL, they are computed here. 
  #    include_sig_eps_prior: If TRUE, sample from the unnormalized log joint posterior, instead of the conditional; 
  #                           see description above for clarification. 
  #
  # Returns:
  #    matrix, of dimension `N_samples` x M. The ith row contains the sampled values at the M input locations. 

  if(is.null(theta_vals_scaled) && is.null(theta_vals_unscaled)) {
    stop("Either `theta_vals_scaled` or `theta_vals_unscaled` must be non-NULL.")
  } else if(is.null(theta_vals_scaled)) {
    theta_vals_scaled <- scale_input_data(theta_vals_unscaled, input_bounds = emulator_info_list$input_bounds)
  } else if(is.null(theta_vals_unscaled)) {
    theta_vals_unscaled <- scale_input_data(theta_vals_scaled, input_bounds = emulator_info_list$input_bounds, inverse = TRUE)
  }
  
  # Log prior evaluations. 
  lprior_vals <- calc_lprior_theta(theta_vals_unscaled, theta_prior_params)
  
  # Compute GP predictive means and variances. 
  if(is.null(gp_pred_list)) {
    gp_pred_list <- predict_independent_GPs(X_pred = theta_vals_scaled, 
                                            gp_obj_list = emulator_info_list$gp_fits,  
                                            gp_lib = emulator_info_list$settings$gp_lib, 
                                            denormalize_predictions = TRUE, 
                                            output_stats = emulator_info_list$output_stats)
  }
  
  # Draw samples from random field approximation of (unnormalized) log posterior density at inputs `theta_vals_unscaled`. 
  lpost_samp <- matrix(NA, nrow = N_samples, ncol = length(gp_pred_list[[1]]$mean))
  
  for(t in seq_len(N_samples)) {
    SSR_samp <- sample_independent_GPs_pointwise(gp_pred_list = gp_pred_list, 
                                                 transformation_methods = emulator_info_list$settings$transformation_method,
                                                 include_nugget = TRUE)
    
    lpost_samp[t,] <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data, 
                                                   SSR = SSR_samp, 
                                                   vars_obs = sig2_eps, 
                                                   na.rm = TRUE, 
                                                   lprior_vals = lprior_vals,
                                                   return_list = FALSE)
  }
  
  # Optionally include prior on likelihood variance parameters, in which case `lpost` refers to the point unnormalized 
  # posterior pi(theta, Sigma), rather than the conditional posterior pi(theta|Sigma). 
  if(include_sig_eps_prior) {
    stop("Inclusion of sig eps prior not yet implemented.") 
  }
  
  if(N_samples == 1) return(lpost_samp[1,])
  
  return(lpost_samp)
  
}
                    

# ------------------------------------------------------------------------------
# lpost emulator functions. 
#    - These functions directly work with the random field approximation to the 
#      unnormalized log posterior density (which we refer to as `lpost`). 
#    - This lpost emulator is induced by the underlying GP emulators (e.g. the 
#      GPs approximating the loss functions). 
# ------------------------------------------------------------------------------

get_lpost_emulator_obj <- function(emulator_info_list, design_info_list, computer_model_data, sig2_eps, theta_prior_params, 
                                   center_output = TRUE, scale_output = TRUE) {
  # Returns a list that represents a fit emulator object, namely the random field approximation to 
  # the unnormalized log posterior density induced by the underlying GPs. The observed input locations 
  # are thus the same as the underlying GPs, while the observed response vector are the unnormalized  
  # log posterior density evaluations computed by `calc_lpost_theta_product_lik()`. Note that the lpost 
  # emulator is also a function of the variance parameters `sig2_eps`. This the lpost emulator object is 
  # primarily used for sequential design purposes, especially batch sequential design heuristics such as 
  # kriging believer and constant liar. For these heuristics, we often want to condition the lpost emulator 
  # on new data, which is different from conditioning the underlying GPs on new data. 
  #
  # Args:
  #    emulator_info_list: list, the emulator information list for the underlying GPs. 
  #    design_info_list: list, the design information list for the underlying GPs. 
  #    computer_model_data: list, the computer model data list. 
  #    sig2_eps: numeric vector of the observation variance parameters. 
  #    theta_prior_params: data.frame, containing prior distributions for the calibration parameters. 
  #    include_nugget: If TRUE, includes nugget variances from underlying GPs when computing lpost inverse kernel matrix. 
  #
  # Returns:
  #    list, the lpost emulator object. The list has named elements "emulator_info_list" (pertaining to 
  #    the underlying GP emulators), "design_info_list" (also pertaining to the underlying GP emulators), 
  #    "outputs_lpost" (the observed response vector for the lpost emulator). The data in 
  #    `lpost_emulator_obj$design_info_list` and `lpost_emulator_obj$inputs_lpost` are allowed to differ, 
  #    as the former pertains to the underlying GPs, while the latter is for the lpost emulator, which 
  #    may be conditioned on new data. 
  
  # Create lpost emulator object and store information relating to underlying GPs. 
  lpost_emulator_obj <- list()
  lpost_emulator_obj$emulator_info_list <- emulator_info_list
  lpost_emulator_obj$design_info_list <- design_info_list
  
  # Store lpost emulator design. 
  lpost_emulator_obj$inputs_lpost <- design_info_list[c("inputs", "inputs_scaled")]
  
  # Compute lpost response vector. 
  lpost_emulator_obj$outputs_lpost <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data, 
                                                                   vars_obs = sig2_eps, na.rm = TRUE, 
                                                                   theta_vals = design_info_list$inputs, 
                                                                   theta_prior_params = theta_prior_params, 
                                                                   return_list = FALSE)
  
  # Optionally normalize response (i.e. unnormalized log posterior density values). 
  lpost_emulator_obj$center_output <- center_output
  lpost_emulator_obj$scale_output <- scale_output
  if(center_output) lpost_emulator_obj$mean_center <- mean(lpost_emulator_obj$outputs_lpost)
  if(scale_output) lpost_emulator_obj$sd_scale <- sd(lpost_emulator_obj$outputs_lpost)
  lpost_emulator_obj$outputs_lpost <- normalize_lpost_outputs(lpost_emulator_obj$outputs_lpost, lpost_emulator_obj)
  
  # lpost emulator is a function of the likelihood parameters (observation variances), prior on the calibration 
  # parameters, and the number of observations (n_obs). 
  lpost_emulator_obj$sig2_eps <- sig2_eps
  lpost_emulator_obj$theta_prior_params <- theta_prior_params
  lpost_emulator_obj$n_obs <- computer_model_data$n_obs

  # Compute lpost inverse kernel matrix. 
  lpost_emulator_obj$K_inv <- chol2inv(chol(calc_lpost_kernel(lpost_emulator_obj, inputs_scaled_1 = design_info_list$inputs_scaled, include_nugget = TRUE)))
  
  # Store log prior density evaluations at design points. 
  lpost_emulator_obj$lprior_design <- calc_lprior_theta(lpost_emulator_obj$inputs_lpost$inputs, theta_prior_params)
  
  # Store prior mean function evaluations at design points. 
  lpost_emulator_obj$mu0_design <- calc_lpost_mean(lpost_emulator_obj, inputs_scaled = lpost_emulator_obj$inputs_lpost$inputs_scaled, 
                                                   lprior_vals = lpost_emulator_obj$lprior_design)
  
  # Store the prior variance. Currently, it is assumed that the prior variance is constant so only a scalar is stored here. This is 
  # in contrast to `mu0_design` above, which is a vector due to the fact that the prior mean function can vary across inputs.
  # `var_prior` excludes nuggets of underlying GPs while `var_comb_prior` includes the nuggets. 
  # TODO: this assumes both the loss emulation approach and hetGP package; should generalize later. 
  lpost_emulator_obj$var_prior <- drop(calc_lpost_kernel(lpost_emulator_obj, inputs_scaled_1=design_info_list$inputs_scaled[1,,drop=FALSE], include_nugget = FALSE))
  lpost_emulator_obj$var_comb_prior <- drop(calc_lpost_kernel(lpost_emulator_obj, inputs_scaled_1=design_info_list$inputs_scaled[1,,drop=FALSE], include_nugget = TRUE))
  
  return(lpost_emulator_obj)
  
}


normalize_lpost_outputs <- function(outputs, lpost_emulator, inverse = FALSE) {
  
  if(inverse) {
    if(lpost_emulator$scale_output) outputs <- outputs * lpost_emulator$sd_scale
    if(lpost_emulator$center_output) outputs <- outputs + lpost_emulator$mean_center
  } else {
    if(lpost_emulator$center_output) outputs <- outputs - lpost_emulator$mean_center
    if(lpost_emulator$scale_output) outputs <- outputs / lpost_emulator$sd_scale
  }

  return(outputs)
  
}


calc_lpost_kernel <- function(lpost_emulator_obj, inputs_scaled_1, inputs_scaled_2 = NULL, include_nugget = TRUE) {
  # Computes the kernel matrix k(`inputs_scaled_1`, `inputs_scaled_2`) for the unnormalized log posterior density emulator,  
  # based on the kernel function induced on lpost by the underlying GPs. If `inputs_scaled_2` is NULL, computes 
  # k(`inputs_scaled_1`, `inputs_scaled_1`). 
  # TODO: this function should be generalized to work with other packages other than "hetGP". 
  #
  # Args:
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    inputs_scaled_1: matrix, of dimension M1 x D where M1 is the number of new inputs and D the dimension of the 
  #                     input parameter space. These inputs should already be properly scaled - no scaling is done 
  #                     in this function.
  #    inputs_scaled_2: matrix, of dimension M2 x D, the second set of scaled inputs. Default is NULL. 
  #    sig2_eps: numeric vector of length P = the number of output variables. The observation variances. 
  #    include_nugget: If TRUE and `inputs_scaled_2` is NULL, includes nugget variances on the diagonal of the underlying 
  #                    GP kernel matrices used to compute the lpost kernel. 
  #    
  # Returns:
  #    The M1 x M2 kernel matrix k(`inputs_scaled_1`, `inputs_scaled_2`). i.e. the covariance matrix (or cross covariance 
  #    matrix) of the lpost emulator between the two sets of inputs. 
  
  if(lpost_emulator_obj$emulator_info_list$settings$emulator_target == "SSR") {
    gp_fits <- lpost_emulator_obj$emulator_info_list$gp_fits
    N_GPs <- length(gp_fits)
    
    GP_scaled_ker_mats <- vector(mode = "list", length = N_GPs)
    for(j in 1:N_GPs) {
      C <- hetGP::cov_gen(X1 = inputs_scaled_1, X2 = inputs_scaled_2, theta = gp_fits[[j]]$theta, type = gp_fits[[j]]$covtype)
      
      if(is.null(inputs_scaled_2)) {
        nug <- rep(gp_fits[[j]]$eps, nrow(C))
        if(include_nugget) nug <- nug + gp_fits[[j]]$g
        C <- hetGP:::add_diag(C, nug)
      }
      
      # Map back to original unnormalized scale and scale by observation variances. 
      output_stats <- lpost_emulator_obj$design_info_list$output_stats
      GP_scaled_ker_mats[[j]] <- output_stats["var_Y", j] * gp_fits[[j]]$nu_hat * C / lpost_emulator_obj$sig2_eps[j]^2                               
    }
    
    K_lpost <- 0.25 * Reduce("+", GP_scaled_ker_mats)
  } else {
    stop("Other emulator targets not yet implemented.")
  }
  
  # Scale according to lpost normalization. 
  if(lpost_emulator_obj$scale_output) K_lpost <- K_lpost / lpost_emulator_obj$sd_scale^2
  
  return(K_lpost)
  
}


calc_lpost_mean <- function(lpost_emulator, inputs_scaled = NULL, inputs_unscaled = NULL, lprior_vals = NULL) {
  # Computes the prior mean function for the unnormalized log posterior density emulator,  
  # based on the mean functions induced on lpost by the underlying GPs.
  # TODO: this function should be generalized to work with other packages other than "hetGP". 
  # TODO: generalize to allow non-constant mean function. 
  #
  # Args:
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    inputs_scaled: matrix, of dimension M x D where M1 is the number of new inputs and D the dimension of the 
  #                   input parameter space. These inputs should already be properly scaled - no scaling is done 
  #                   in this function. Currently, this argument has no impact on the function as the prior mean 
  #                   is constant, but it could play a role if more general mean functions are specified. 
  #   
  # Returns:
  #    The M1 x M2 kernel matrix k(`inputs_scaled_1`, `inputs_scaled_2`). i.e. the covariance matrix (or cross covariance 
  #    matrix) of the lpost emulator between the two sets of inputs. 
  
  # Requires scaled and unscaled inputs. 
  if(is.null(inputs_scaled) && is.null(inputs_unscaled)) {
    stop("Either scaled or unscaled inputs must be provided.")
  } else if(is.null(inputs_scaled)) {
    inputs_scaled <- scale_input_data(inputs_unscaled, input_bounds = lpost_emulator$design_info_list$input_bounds)
  } 
  
  # Underlying GP means, de-normalize to return to original scale. 
  GP_means <- sapply(lpost_emulator$emulator_info_list$gp_fits, function(gp) gp$beta0)
  output_stats <- lpost_emulator$design_info_list$output_stats
  GP_means <- drop(normalize_output_data(matrix(GP_means, nrow=1), output_stats, inverse = TRUE))
  
  # Prior on calibration parameters. 
  if(is.null(lprior_vals)) {
    if(is.null(inputs_unscaled)) inputs_unscaled <- scale_input_data(inputs_scaled, input_bounds = lpost_emulator$design_info_list$input_bounds, inverse = TRUE)
    lprior_vals <- calc_lprior_theta(theta_vals = inputs_unscaled, theta_prior_params = lpost_emulator$theta_prior_params)
  }
    
  # Prior mean function induced on lpost emulator. Note that since the underlying GPs are (for now) assumed to have constant prior 
  # mean functions, that the only potential variation in the lpost prior mean function comes from the log prior evaluations. 
  # Also scale mean function according to the lpost centering transformation. 
  means <- -0.5 * (sum(lpost_emulator$n_obs * log(2*pi*lpost_emulator$sig2_eps)) + sum(GP_means / lpost_emulator$sig2_eps)) + lprior_vals
  
  # Scale according to lpost normalization. 
  means <- normalize_lpost_outputs(means, lpost_emulator)
  
  return(matrix(means, ncol = 1))
  
}


update_lpost_emulator <- function(lpost_emulator, inputs_new_scaled, outputs_lpost_new = NULL, 
                                  inputs_new_unscaled = NULL, outputs_normalized = FALSE, verbose = TRUE) {
  # Updates the random field approximation to the unnormalized log posterior density induced 
  # by the underlying GPs by conditioning on newly observed data {`input_new`, `output_lpost_new`}. 
  # In the primary use case of this function, `output_lpost_new` is "pseudo-data" used for heuristic 
  # batch sequential design updates such as kriging believer and constant liar. If `output_lpost_new` 
  # is NULL, then it will be set to the predictive mean of the lpost emulator at inputs `input_new`.
  # This has the effect of updating the lpost emulator predictive variance without affecting the
  # predictive mean. No hyperparameter re-estimation is performed here, only conditioning on data. 
  # Also, no updating of the underlying GPs is performed here. This is because this function is 
  # mainly used for conditioning on pseudo-data for batch heuristic strategies, rather than 
  # actually updating the emulator on true observed data. The updates on `lpost_emulator` include 
  # appending new data to `lpost_emulator$inputs_lpost` and `lpost_emulator`$outputs_lpost` as 
  # well as updating the inverse kernel matrix via partitioned matrix inverse equations. 
  # 
  # Args:
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    input_new_scaled: matrix, of dimension M x D where M is the number of new inputs and D the dimension of  
  #               the input parameter space. These inputs should already be properly scaled - no scaling is done 
  #               in this function. 
  #    output_lpost_new: numeric, the vector of length M of lpost response values associated with the inputs 
  #                      `input_new`. If NULL, this is set to the expected value of the lpost emulator at 
  #                      as inputs `input_new`. 
  #    inputs_new_unscaled: Optionally provide unscaled version of `inputs_new_scaled` as well, which is required to compute  
  #                         predictive means (as the predictive mean depends on the prior distribution on the calibration 
  #                         parameters). 
  #    outputs_normalized: logical(1), indicates whether or not `outputs_lpost_new` has already been centered (according to the 
  #                        centering transformation defined by `lpost_emulator$mean_center`).
  #
  # Returns:
  #    The updated lpost emulator object. 
  
  # Remove any repeated observations. If the only new observation is removed, return `lpost_emulator` unchanged. 
  # TODO: think better about how to handle this; does it even matter? 
  # repeated_obs_sel <- inputs_new_scaled %in% lpost_emulator$inputs_lpost$inputs_scaled
  # if(any(repeated_obs_sel)) {
  #   if(verbose) message("Removing repeated input(s) prior to updating emulator; indices: ", seq_along(inputs_new_scaled)[repeated_obs_sel])
  # }
  # inputs_new_scaled <- inputs_new_scaled[!repeated_obs_sel,, drop = FALSE]
  # if(length(inputs_new_scaled) == 0) return(lpost_emulator)
  
  if(is.null(inputs_new_unscaled)) {
    inputs_new_unscaled <- scale_input_data(inputs_new_scaled, input_bounds = lpost_emulator$emulator_info_list$input_bounds, inverse = TRUE)
  } else {
    # TODO: update this
    # inputs_new_unscaled <- inputs_new_unscaled[!repeated_obs_sel,, drop = FALSE] 
  }
  
  # If no new responses are provided, set to GP expectation. 
  if(is.null(outputs_lpost_new)) {
    outputs_lpost_new <- predict_lpost_emulator(inputs_new_scaled, lpost_emulator = lpost_emulator, return_vals = "mean",
                                                inputs_new_unscaled = inputs_new_unscaled, unscale = FALSE, uncenter = FALSE)$mean
  } else {
    # TODO: update this
    # outputs_lpost_new <- outputs_lpost_new[!repeated_obs_sel]
    if(!outputs_normalized) outputs_lpost_new <- normalize_lpost_outputs(outputs_lpost_new, lpost_emulator)
  }
   
  # Add log prior density evaluations at new points. 
  lprior_vals_new <- calc_lprior_theta(inputs_new_unscaled, lpost_emulator$theta_prior_params)
  lpost_emulator$lprior_design <- c(lpost_emulator$lprior_design, lprior_vals_new)

  # Add prior mean function evaluations at new points. Note that prior variance of lpost emulator is assumed constant 
  # so no update is needed for `lpost_emulator$var_prior` and `lpost_emulator$var_comb_prior`. 
  lpost_emulator$mu0_design <- rbind(lpost_emulator$mu0_design, calc_lpost_mean(lpost_emulator, inputs_scaled = inputs_new_scaled, lprior_vals = lprior_vals_new))
                                                                                
  # Update inverse kernel matrix and design data. 
  for(i in 1:nrow(inputs_new_scaled)) {
    
    # Inverse kernel matrix. 
    lpost_emulator$K_inv <- update_lpost_inverse_kernel_matrix(lpost_emulator, input_new_scaled=inputs_new_scaled[i,,drop=FALSE])
                                                               
    # Design data. 
    lpost_emulator$inputs_lpost$inputs_scaled <- rbind(lpost_emulator$inputs_lpost$inputs_scaled, inputs_new_scaled[i,,drop=FALSE])
    lpost_emulator$inputs_lpost$inputs <- rbind(lpost_emulator$inputs_lpost$inputs, inputs_new_unscaled[i,,drop=FALSE])
    lpost_emulator$outputs_lpost <- c(lpost_emulator$outputs_lpost, outputs_lpost_new[i])
  }
  
  return(lpost_emulator)
  
}


predict_lpost_emulator <- function(inputs_new_scaled, lpost_emulator, return_vals = c("mean", "var"), 
                                   inputs_new_scaled_2 = NULL, include_nugget = TRUE, 
                                   inputs_new_unscaled = NULL, prior_mean_vals_new = NULL, 
                                   verbose = TRUE, unscale = TRUE, uncenter = TRUE) {
  # Compute predictive mean and variance equations for the `lpost_emulator`. 
  # Importantly, note that this function will potentially yield different predictions than `predict_lpost_GP_approx()`. 
  # The latter always generates predictions based on the underlying GPs, encapsulated in `emulator_info_list`. On the 
  # other hand, this function predicts based on the `lpost_emulator` object, which may have been conditioned on new 
  # data without affecting the underlying GPs. The reason for this discrepancy is to facilitate conditioning the lpost 
  # emulator on pseudo data, as is required in batch sequential design heuristics such as kriging believer and constant 
  # liar. Even with the same GP hyperparameters, this function will produce different predictions than `predict_lpost_GP_approx()`. 
  # This is due to the fact that this function directly operates on the lpost prior induced by the underlying GPs, and then 
  # conditions on the observed log posterior values. On the other hand, `predict_lpost_GP_approx()` conditions each GP 
  # individually on the observed SSR values and then plugs the resulting predictive GP distributions in the log posterior 
  # expression. Note that `predict_lpost_emulator()` assumes the prior GP distribution is covariance stationary; in particular, 
  # the prior variance is constant. 
  #
  # Args:
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    input_new_scaled: matrix, of dimension M x D where M is the number of inputs at which to predict and 
  #                      D the dimension of the input parameter space. These inputs should already be properly scaled -  
  #                      no scaling is done in this function.
  #                      Alternatively, a numeric vector of length D, which will be converted to a 1 x D matrix. 
  #    include_nugget: If TRUE, includes nugget variances on underlying GP kernel computations used to compute 
  #                    the lpost kernel. This implies the returned predictive variance is the variance of the 
  #                    latent function plus the nugget variance. Otherwise, the returned variance excludes the  
  #                    nugget variance. 
  #    return_vals: character, either "mean", "var", or c("mean", "var") depending on whether the predictive mean
  #                 or variance is desired, or both.
  #    inputs_new_unscaled: Optionally provide unscaled version of `input_new_scaled` as well, which is required to compute  
  #                         predictive means (as the predictive mean depends on the prior distribution on the calibration 
  #                         parameters). 
  #    prior_mean_vals_new: numeric, vector of evaluations of the lpost prior mean function at inputs `inputs_new_scaled`. 
  #                         If predicting at the same set of points repeatedly, pre-computing these mean function 
  #                         evaluations can result in very large speedups. Default is NULL, in which case the mean function 
  #                         evaluations are computed in the function. 
  #    denormalize: if TRUE, undoes the centering transformation defined by `lpost_emulator$mean_center`. 
  #
  # Returns:
  #    list, with names "mean" and "var", containing respectively the predictive mean and variance values computed 
  #    at the test inputs `inputs_new_scaled`. The associated value with be NULL if it is not included in 
  #    `return_vals`. 
  
  return_list <- list(mean = NULL, var = NULL, cov = NULL)
  
  # Cross covariances. 
  kn <- calc_lpost_kernel(lpost_emulator, inputs_scaled_1 = inputs_new_scaled, inputs_scaled_2 = lpost_emulator$inputs_lpost$inputs_scaled)
           
  # Predictive mean. 
  if("mean" %in% return_vals) {
    # Prior mean function evaluations at new points. 
    if(is.null(prior_mean_vals_new)) prior_mean_vals_new <- calc_lpost_mean(lpost_emulator, inputs_scaled = inputs_new_scaled, inputs_unscaled = inputs_new_unscaled)

    return_list$mean <- drop(prior_mean_vals_new + kn %*% (lpost_emulator$K_inv %*% (matrix(lpost_emulator$outputs_lpost, ncol=1) - lpost_emulator$mu0_design)))
  }
  
  # Predictive variance. 
  if("var" %in% return_vals) {
    prior_var <- ifelse(include_nugget, lpost_emulator$var_comb_prior, lpost_emulator$var_prior)
    pred_vars <- as.vector(prior_var - hetGP:::fast_diag(kn, tcrossprod(lpost_emulator$K_inv, kn)))

    # Ordinary kriging variance correction. Assumes underlying GPs either all use ordinary kriging or all use simple kriging.
    if(lpost_emulator$emulator_info_list$gp_fits[[1]]$trendtype == "OK") {
      pred_vars <- pred_vars + (1 - tcrossprod(rowSums(lpost_emulator$K_inv), kn))^2/sum(lpost_emulator$K_inv)
    }

    # Set negative variances due to numerical rounding errors to 0.
    neg_var_sel <- (pred_vars < 0)
    if(any(neg_var_sel)) {
      if(verbose) message("Warning: `predict_lpost_emulator()` setting negative variances to 0.")
      pred_vars[neg_var_sel] <- 0
    }
    
    return_list$var <- drop(pred_vars)
    
  }
  
  # Predictive Covariance. 
  if("cov" %in% return_vals) {
    if(is.null(inputs_new_scaled_2)) stop("`inputs_new_scaled` is required to compute predictive covariance.")
    
    # Cross covariance between design points and second set of inputs. 
    kn2 <- calc_lpost_kernel(lpost_emulator, inputs_scaled_1 = lpost_emulator$inputs_lpost$inputs_scaled, inputs_scaled_2 = inputs_new_scaled_2)
    
    # Predictive covariance matrix. 
    if(nrow(inputs_new_scaled) < nrow(inputs_new_scaled_2)) {
      pred_cov <- calc_lpost_kernel(lpost_emulator, inputs_new_scaled, inputs_new_scaled_2) - (kn %*% lpost_emulator$K_inv) %*% kn2
    } else {
      pred_cov <- calc_lpost_kernel(lpost_emulator, inputs_new_scaled, inputs_new_scaled_2) - kn %*% (lpost_emulator$K_inv %*% kn2)
    }
    
    # Ordinary kriging variance correction. Assumes underlying GPs either all use ordinary kriging or all use simple kriging.
    if(lpost_emulator$emulator_info_list$gp_fits[[1]]$trendtype == "OK") {
      pred_cov <- pred_cov + crossprod(1 - tcrossprod(rowSums(lpost_emulator$K_inv), kn), 1 - rowSums(lpost_emulator$K_inv) %*% kn2)/sum(lpost_emulator$K_inv)
    }
    
    return_list$cov <- pred_cov
    
  }
  
  # Optionally invert lpost normalization. 
  if(unscale || uncenter) {
    
    if(unscale && lpost_emulator$scale_output) {
      if("var" %in% return_vals) return_list$var <- lpost_emulator$sd_scale^2 * return_list$var
      if("mean" %in% return_vals) return_list$mean <- lpost_emulator$sd_scale * return_list$mean
      if("cov" %in% return_vals) return_list$cov <- lpost_emulator$sd_scale^2 * return_list$cov
    }
    
    if(uncenter && lpost_emulator$center_output) return_list$mean <- lpost_emulator$mean_center + return_list$mean

  }
  
  return(return_list)
  
}


update_lpost_inverse_kernel_matrix <- function(lpost_emulator, input_new_scaled) {
  # Uses partitioned matrix inverse equations to perform a fast update of the inverse kernel matrix for the lpost 
  # emulator given a new set of input point, `input_new_scaled`. Currently only considers the addition of a single input 
  # point as this is the primary use case. Inverse kernel matrix update corresponding to conditioning on multiple inputs 
  # can be performed by calling this function one point at a time. Credit to `hetGP` function `update_Ki()` for 
  # providing a guide for this implementation.
  #
  # Args:
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    input_new_scaled: matrix, of dimension 1 x D where  D the dimension of the input parameter space. 
  #                      This input should already be properly scaled - no scaling is done in this function.
  #                      Alternatively, a numeric vector of length D, which will be converted to a 1 x D matrix. 
  #
  # Returns: 
  #    The updated inverse kernel matrix. 
  
  # Convert to matrix, if necessary. 
  if(!is.matrix(input_new_scaled)) input_new_scaled <- matrix(input_new_scaled, nrow = 1)
  
  # Construct matrix inverse via partitioned inverse equations. 
  k_new_old <- calc_lpost_kernel(lpost_emulator, inputs_scaled_1=input_new_scaled, inputs_scaled_2=lpost_emulator$inputs_lpost$inputs_scaled)
  k_new <- lpost_emulator$var_comb_prior # Since prior lpost variance is assumed stationary. 
  
  K_inv_k_new_old <- tcrossprod(lpost_emulator$K_inv, k_new_old)
  nu_inv <- 1 / drop(k_new - k_new_old %*% K_inv_k_new_old)
  z <- -nu_inv * K_inv_k_new_old
  K_inv_top_left <- lpost_emulator$K_inv + (1/nu_inv) * tcrossprod(z)
  
  return(rbind(cbind(K_inv_top_left, z), c(z, nu_inv)))

}


sample_lpost_emulator <- function(inputs_new_scaled, lpost_emulator, N_samples = 1, inputs_new_unscaled = NULL, 
                                  include_nugget = TRUE, prior_mean_vals = NULL, verbose = TRUE, denormalize = TRUE) {
  # Samples from the GP predictive distribution of the lpost emulator (which is assumed to be a GP). 
  #
  # Args:
  #    theta_vals_scaled: matrix of dimension M x D of input locations (scaled to lie in unit hypercube) at which to sample
  #                       the log density values. Each row is a location. 
  #    lpost_emulator: list, the lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    N_samples: integer, the number of samples to draw at each input location. 
  #    theta_vals_unscaled: matrix, of dimension M x D; the same as `theta_vals_scaled` but the inputs are unscaled. Both the 
  #                       scaled and unscaled inputs are required here, but only one of them need be passed, as the scaling 
  #                       can be performed using the information in `emulator_info_list$input_bounds`. 
  #    include_nugget: If TRUE, includes nugget variances on underlying GP kernel computations used to compute 
  #                    the lpost kernel. This implies the returned predictive variance is the variance of the 
  #                    latent function plus the nugget variance. Otherwise, the returned variance excludes the  
  #                    nugget variance. 
  #    prior_mean_vals: numeric, vector of evaluations of the lpost prior mean function at inputs `inputs_new_scaled`. 
  #                     See `predict_lpost_emulator()`. Default is NULL, in which case the mean function 
  #                     evaluations are computed in the function. 
  #    denormalize: if TRUE, undoes the centering transformation defined by `lpost_emulator$mean_center`. 
  #
  # Returns:
  #    matrix, of dimension `N_samples` x M. The ith row contains the sampled values at the M input locations. 
  
  # Compute predictive mean and variance. 
  lpost_pred <- predict_lpost_emulator(inputs_new_scaled = inputs_new_scaled, lpost_emulator = lpost_emulator, return_vals = c("mean", "var"), 
                                       include_nugget = include_nugget, inputs_new_unscaled = inputs_new_unscaled, 
                                       prior_mean_vals_new = prior_mean_vals, verbose = verbose)
                                     
  # Draw pointwise samples from marginal predictive distributions. 
  samp <- matrix(nrow = N_samples, ncol = nrow(inputs_new_scaled))
  for(t in seq_len(N_samples)) {
    samp[t,] <- sample_GP_pointwise(lpost_pred$mean, lpost_pred$var) 
  }
  
  if(denormalize) samp <- normalize_lpost_outputs(samp, lpost_emulator, inverse = TRUE)
  
  if(N_samples == 1) return(samp[1,])
  return(samp)
  
}


convert_to_post_emulator_log_moments <- function(lpost_emulator, lpost_means, lpost_vars, return_vals = c("log_mean", "log_var")) {
  
  return_list <- list(log_mean = NULL, log_var = NULL)
  
  # Log Mean. 
  if("log_mean" %in% return_vals) return_list$log_mean <- lpost_means + 0.5 * lpost_vars
  
  # Log Variance. 
  if("log_var" %in% return_vals) {
    
    # Numerically stable computation of log exponential term. 
    idx_approx_sel <- (lpost_vars >= 100)
    log_exp_term <- vector(mode = "numeric", length = length(lpost_vars))
    log_exp_term[!idx_approx_sel] <- log(exp(lpost_vars[!idx_approx_sel]) - 1)
    log_exp_term[idx_approx_sel] <- lpost_vars[idx_approx_sel] # Apply approximation. 
    
    return_list$log_var <- 2 * lpost_means + lpost_vars + log_exp_term
    
  }
  
  return(return_list)
  
}


get_lpost_emulator_metric_comparison <- function(lpost_emulator_list, lpost_emulator_validation_list, metrics, include_nugget = TRUE) {
  # Computes emulator metrics for a list of different lpost emulators on a single validation set. Note that the list 
  # `lpost_emulator_validation_list` is in one-to-one correspondence with `lpost_emulator_list`, due to the fact that 
  # the input scaling and likelihood parameters may differ by emulator, resulting in different scaled inputs and different 
  # outputs in the validation sets. 
  #
  # Args:
  #    lpost_emulator_list: list, of lpost emulator objects, as returned by `get_lpost_emulator_obj()`. 
  #    lpost_emulator_validation_list: list, of equal length as `lpost_emulator_list` containing validation data information 
  #                                    for each emulator. A list of the required form is returned by 
  #                                    `format_lpost_emulator_validation_data()`.
  #    metrics: character, vector of metric names; used to call function "get_gp_<metric>". 
  #    include_nugget: logical, whether or not to include nugget in computing lpost emulator metrics. 
  #
  # Returns:
  #    data.frame, with column names "emulator_name", "metric_name", and "metric_value". The emulator names are taken from the 
  #    names attribute of the list `lpost_emulator_list`, is non-NULL. 
  
  # Labels for each emulator. 
  emulator_names <- names(lpost_emulator_list)
  if(is.null(emulator_names)) emulator_names <- paste0("emulator", 1:length(lpost_emulator_list))
  
  # data.frame to store computed metrics. 
  metrics_out <- data.frame(emulator_name = character(), metric_name = character(), metric_value = numeric())
  
  for(j in seq_along(lpost_emulator_list)) {
    
    # Compute validation metrics. 
    metrics_results <- get_lpost_emulator_metrics(lpost_emulator = lpost_emulator_list[[j]], 
                                                  lpost_validation_inputs = lpost_emulator_validation_list[[j]]$inputs, 
                                                  lpost_validation_inputs_scaled = lpost_emulator_validation_list[[j]]$inputs_scaled,
                                                  lpost_validation_outputs = lpost_emulator_validation_list[[j]]$outputs,
                                                  metrics = metrics,  
                                                  include_nugget = include_nugget)
    
    # Append to data.frame. 
    metrics_out <- rbind(metrics_out, data.frame(emulator_name = emulator_names[j], metric_name = names(metrics_results), metric_value = metrics_results))
    
  }
  
  return(metrics_out)

}
                                                   

get_lpost_emulator_metrics <- function(lpost_emulator, lpost_validation_inputs_scaled, lpost_validation_outputs, metrics, 
                                       lpost_validation_inputs = NULL, include_nugget = TRUE) {
  # TODO: generalize this to allow computing metrics that require predictive covariance matrix. 
  # Computes lpost emulator metrics by first computing lpost emulator mean and variance predictions at a set of validation inputs, 
  # and then comparing these predictive quantities to the true observed outputs `lpost_validation_outputs`. 
  #
  # Args:
  #    lpost_emulator: list, lpost emulator object, as returned by `get_lpost_emulator_obj()`. 
  #    lpost_validation_inputs_scaled: matrix, of dimension M x D where M is the number of scaled validation inputs and D the dimension of 
  #                                    the input space. 
  #    lpost_validation_outputs: numeric, vector of length M. The observed unnormalized log posterior evaluations used as the ground 
  #                              truth when computing emulator metrics. 
  #    metrics: character, vector of metric names; used to call function "get_gp_<metric>". 
  #    lpost_validation_inputs: matrix, of the same dimension as `lpost_validation_inputs_scaled`. The unscaled validation inputs. 
  #    include_nugget: logical, whether or not to include nugget in computing lpost emulator metrics.
  #
  # Returns:
  #    numeric, named vector of length `length(metrics)`. The names are set to `names(metrics)` and the values to the corresponding 
  #    computed matric values. 
  
  # Compute lpost emulator predictions at validation inputs. 
  pred <- predict_lpost_emulator(lpost_validation_inputs_scaled, lpost_emulator, return_vals = c("mean", "var"), include_nugget = include_nugget, 
                                 inputs_new_unscaled = lpost_validation_inputs)
  
  # Compute metrics using emulator predictions and true outputs. 
  metrics_results <- compute_gp_metrics(gp_mean_pred = pred$mean, gp_var_pred = pred$var, output_true = lpost_validation_outputs, metrics = metrics)
  
  return(metrics_results)
  
}


format_lpost_emulator_validation_data <- function(lpost_emulator_list, inputs_unscaled, computer_model_data) {
  
  # Create one validation data info object for each lpost emulator. 
  validation_info_list <- vector(mode = "list", length = length(lpost_emulator_list))
  names(validation_info_list) <- names(lpost_emulator_list)

  # Compute L2 error corresponded to inputs `inputs_unscaled`. This error will be the same for all emulators. 
  # However the true lpost values may differ by emulator due to different values of likelihood parameters. 
  SSR_outputs <- get_computer_model_SSR(computer_model_data, theta_vals = inputs_unscaled, na.rm = TRUE)
  
  for(j in seq_along(lpost_emulator_list)) {
    
    # Store both scaled and unscaled inputs. The scaling transformation may differ from emulator to emulator, hence the need 
    # to have a different validation data info object for each.
    validation_info_list[[j]][["inputs"]] <- inputs_unscaled
    validation_info_list[[j]][["inputs_scaled"]] <- scale_input_data(inputs_unscaled, input_bounds = lpost_emulator_list[[j]]$emulator_info_list$input_bounds)
    
    # Run full forward model to obtain true lpost values corresponding to `inputs_unscaled`.
    validation_info_list[[j]][["outputs"]] <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data,
                                                                           SSR = SSR_outputs,  
                                                                           vars_obs = lpost_emulator_list[[j]]$sig2_eps, na.rm = TRUE, 
                                                                           theta_vals = inputs_unscaled, 
                                                                           theta_prior_params = lpost_emulator_list[[j]]$theta_prior_params, 
                                                                           return_list = FALSE)
  }
  
  return(validation_info_list)
  
}





