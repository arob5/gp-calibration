#
# mcmc_helper_functions.r
# Primarily functions to post-process, organize, and plot MCMC samples. 
#
# Andrew Roberts
# 

library(HDInterval)
library(kde1d)

# ------------------------------------------------------------------------------
# MCMC Formatting Functions. 
# ------------------------------------------------------------------------------

format_mcmc_output <- function(samp_list, test_label) {
  # This function re-formats the output returned by a single MCMC run. This MCMC run is identified with 
  # the label `test_label`. The function combines matrices containing MCMC samples of different parameter 
  # types (see `samp_list` description below) into a single data.table. The returned data.table is in a 
  # long format in that it has columns "test_label", "param_type", "param_name", "itr", and "sample". 
  #
  # Args:
  #    samp_list: named list, each element must have name set to the relevant parameter type (e.g. "theta" or "sig_eps"). 
  #               Each element is a matrix where rows contain MCMC samples of that parameter type. The column names of 
  #               these matrices must be set to the parameter names. Given that the samples are all from the same MCMC run, 
  #               the matrices will all typically have the same number of rows. However, this is not required. 
  #    test_label: character, a string providing a label for the MCMC run. 
  #
  # Returns: 
  #    data.table with column names "test_label", "param_type", "param_name", "itr", and "sample". The column "test_label" will be 
  #    constant with value set to `test_label`. The values in column "param_type" are taken from the names of `samp_list`, 
  #    while the values in column "param_name" are taken from the column names in the matrices within `samp_list`. The 
  #    column `itr` contains the integer MCMC iteration. The column `sample` contains the MCMC numeric sample values. 
  
  samp_dt <- data.table(param_type=character(), itr=integer(), param_name=character(), sample=numeric())
  
  for(j in seq_along(samp_list)) {
    # Format samples for current variable. 
    samp_param_dt <- as.data.table(samp_list[[j]])
    
    if(nrow(samp_param_dt) > 0) {
      samp_param_dt[, param_type := names(samp_list)[j]]
      samp_param_dt[, itr := 1:.N]
      samp_param_dt <- melt.data.table(data=samp_param_dt, id.vars=c("param_type", "itr"), 
                                       variable.name="param_name", value.name="sample", na.rm=TRUE, 
                                       variable.factor=FALSE)
      
      # Append to samples for existing variables. 
      samp_dt <- rbindlist(list(samp_dt, samp_param_dt), use.names=TRUE)  
    }
  }
  
  # Add test label. 
  samp_dt[, test_label := test_label]
  
  
  return(samp_dt)
  
}


append_mcmc_output <- function(mcmc_samp_dt, samp_list, test_label) {
  # Appends a new `samp_list` to an existing MCMC samp data.table. First 
  # formats the list via `format_mcmc_output()` then appends it. 
  # Ensures that the `test_label` is not already present in the data.table, 
  # but this could be relaxed if needed in the future. 
  
  assert_that(is.data.table(mcmc_samp_dt))
  assert_that(setequal(colnames(mcmc_samp_dt), 
                       c("param_type", "itr", "param_name", "sample", "test_label")))
  assert_that(!(test_label %in% mcmc_samp_dt$test_label))
  
  mcmc_samp_dt_new <- format_mcmc_output(samp_list, test_label=test_label)
  mcmc_samp_dt <- rbindlist(list(mcmc_samp_dt, mcmc_samp_dt_new), use.names=TRUE)
  
  return(mcmc_samp_dt)
  
}


select_mcmc_samp <- function(samp_dt, burn_in_start=NULL, test_labels=NULL, param_types=NULL, param_names=NULL,
                             return_burnin=FALSE) {
  # Operates on the long data.table format, as returned by `format_mcmc_output()`. Selects rows corresponding to 
  # valid combinations of `test_labels`, `param_types`, `param_names`. Also removes iterations specified as "burn-in". 
  # See `burn_in_start` for details. 
  #
  # Args:
  #    samp_dt: data.table of MCMC samples, in long format as returned by format_mcmc_output()`. 
  #    burn_in_start: If NULL, selects all MCMC iterations. If integer of length 1, this is interpreted as the starting 
  #                   iteration for all parameters - all earlier iterations are dropped. If vector of length > 1, must 
  #                   be a named vector with names set to valid test label values. This allows application of a different 
  #                   burn-in start iteration for different test labels. 
  #    test_labels, param_types, param_names: vectors of values to include in selection corresponding to columns 
  #                                           "test_label", "param_type", and "param_name" in `samp_dt`. A NULL  
  #                                            value includes all values found in `samp_dt`. 
  #    return_burnin: If TRUE, returns only the burn-in iterations (i.e. all iterations strictly less
  #                   than the burn_in_start values). Otherwise, the default behavior is to drop 
  #                   the burn-in, returning all iterations greater than or equal to the 
  #                   burn_in_start values. 
  #
  # Returns:
  #    data.table, containing subset of rows from `samp_dt`. 
  
  samp_dt_subset <- copy(samp_dt)
  
  # If not provided, select all. 
  if(is.null(test_labels)) test_labels <- samp_dt_subset[, unique(test_label)]
  if(is.null(param_types)) param_types <- samp_dt_subset[, unique(param_type)]
  if(is.null(param_names)) param_names <- samp_dt_subset[, unique(param_name)]
  
  # Select rows corresponding to label-type-name combinations. 
  samp_dt_subset <- samp_dt[(test_label %in% test_labels) & 
                              (param_type %in% param_types) & 
                              (param_name %in% param_names)]
  
  # Remove (or select) burn-in iterations; burn-in can differ by test label. 
  samp_dt_subset <- remove_mcmc_samp_burnin(samp_dt_subset, burn_in_start, return_burnin)
  
  return(samp_dt_subset)
  
}


select_mcmc_samp_mat <- function(samp_dt, test_label, param_type, param_names=NULL, burn_in_start=NULL, 
                                 return_burnin=FALSE) {
  # Converts MCMC samples from long to wide format. Wide format means that a 
  # matrix is returned where each row is a sample. A `param_type` (e.g. "theta") must be selected so 
  # that each row of the matrix is a valid parameter vector (e.g. each row is a sampled parameter 
  # calibration vector). Alternatively, a subset of the parameters for a specific parameter type 
  # can be returned by passing the `param_names` argument. It is recommended to pass this argument 
  # even if all parameters of a specific type are to be returned, as passing `param_names` will 
  # order to the columns according to `param_names`, which ensures that the ordering is correct. 
  #
  # Args:
  #    samp_dt: data.table, must be of the format described in `format_mcmc_output()`. 
  #    test_label: character(1), the single test label to select. 
  #    param_type: character(1), the single parameter type to select. Must be a type within the  
  #                specified `test_label`. 
  #    param_names: character, if provided then selects the specific parameter names provided, and 
  #                also orders the columns of the returned matrix in the order they are provided
  #                in this argument. The parameters specified here must be of the type specified 
  #                by `param_type` and contained within the test label specified by `test_label`. 
  #                If NULL, selects all parameters within the param type, and no explicit ordering 
  #                is performed. 
  #    burn_in_start: integer(1), the starting iteration to use for the samples from the test label. 
  #                   If NULL, does not drop any burn-in. 
  #    return_burnin: If TRUE, returns only the burn-in iterations (i.e. all iterations strictly less
  #                   than the burn_in_start values). Otherwise, the default behavior is to drop 
  #                   the burn-in, returning all iterations greater than or equal to the 
  #                   burn_in_start values. 
  #
  # Returns:
  #    matrix, where each row is a sample from the selected parameters, with burn-in dropped or 
  #    returned if specified. Rownames are set to the iteration numbers. 
  
  if(length(test_label)>1 || length(param_type)>1) stop("Must select single test label and param type.")
  
  samp <- select_mcmc_samp(samp_dt, burn_in_start, test_label, param_type, param_names, return_burnin)
  if(nrow(samp)==0) stop("Cant convert zero length data.table to wide matrix format.")
  
  samp <- dcast(data=samp, formula=itr~param_name, value.var="sample")
  itrs <- samp$itr
  samp_mat <- as.matrix(samp[, .SD, .SDcols=!"itr"])
  rownames(samp_mat) <- itrs
  
  if(!is.null(param_names)) samp_mat <- samp_mat[,param_names, drop=FALSE]
  
  return(samp_mat)
  
}


remove_mcmc_samp_burnin <- function(samp_dt, burn_in_start, return_burnin=FALSE) {
  # Acts on a data.table in long format containing MCMC samples and either drops 
  # burn-in iterations or returns only the burn-in iterations. Burn-in can either 
  # differ by test label or a single burn-in cutoff can be applied to all test 
  # labels. 
  #
  # Args:
  #    samp_dt: data.table, must be of the format described in `format_mcmc_output()`. 
  #    burn_in_start: If NULL, selects all MCMC iterations. If integer of length 1 AND unnamed, this is 
  #                   interpreted as the starting iteration for all parameters - all earlier iterations
  #                   are dropped. If a named vector, must  be a named vector with names set
  #                   to valid test label values. This allows application of a different 
  #                   burn-in start iteration for different test labels. If only some test labels are 
  #                   specified by the vector, then the others will not be affected. 
  #    return_burnin: If TRUE, returns only the burn-in iterations (i.e. all iterations strictly less
  #                   than the burn_in_start values). Otherwise, the default behavior is to drop 
  #                   the burn-in, returning all iterations greater than or equal to the 
  #                   burn_in_start values. 
  #
  # Returns:
  #    data.table, a copy of `samp_dt` which is row subsetted to either drop the burn-in or 
  #    non-burn-in iterations, depending on the value of `return_burnin`. 
  
  samp_dt_subset <- copy(samp_dt)
  
  # If `burn_in_start` is NULL, either returns entire data.table or empty data.table. 
  if(is.null(burn_in_start)) {
    if(return_burnin) return(samp_dt_subset[0])
    return(samp_dt_subset)
  }
  
  # If `burn_in_start` unnamed vector of length 1, apply to all test labels. 
  if((length(burn_in_start) == 1) && is.null(names(burn_in_start))) {
    if(return_burnin) return(samp_dt_subset[itr < burn_in_start]) 
    else return(samp_dt_subset[itr >= burn_in_start]) 
  }
  
  # If `burn_in_start` is named vector, apply burn-in values separately for each test label. 
  test_labels <- samp_dt_subset[, unique(test_label)]
  extra_labels <- setdiff(names(burn_in_start), test_labels)
  missing_labels <- setdiff(test_labels, names(burn_in_start))
  if(length(extra_labels) > 0) message("Extra labels not used: ", paste(extra_labels, collapse=", "))
  if(length(missing_labels) > 0) message("Labels not affected: ", paste(missing_labels, collapse=", "))
  test_labels <- setdiff(test_labels, missing_labels)
  
  for(lbl in test_labels) {
    if(return_burnin) samp_dt_subset <- samp_dt_subset[(test_label != lbl) | (itr < burn_in_start[lbl])]
    else samp_dt_subset <- samp_dt_subset[(test_label != lbl) | (itr >= burn_in_start[lbl])]
  }
  
  return(samp_dt_subset)
}


# ------------------------------------------------------------------------------
# Computing statistics and errors from MCMC samples. 
# ------------------------------------------------------------------------------

compute_mcmc_param_stats <- function(samp_dt, burn_in_start=NULL, test_labels=NULL, param_types=NULL,
                                      param_names=NULL, subset_samp=TRUE, format_long=FALSE) {
  # Currently just computes sample means and variances for the selected parameters/variables in `samp_dt`. 
  #
  # Args:
  #   samp_dt: data.table, must be of the format described in `format_mcmc_output()`
  #   burn_in_start, param_types, param_names: passed to `select_mcmc_samp()` to subset `samp_dt` to determine which
  #                                            parameters and samples will be included in the metric computations.
  #   subset_samp: If FALSE, indicates that `samp_dt` is already in the desired form so do not call `select_mcmc_samp()`. 
  #   format_long: If FALSE (the default), then one column is created for each statistic computed with column names 
  #                of the form `stat_<stat_name>`. Otherwise, two columns are added: "stat_name" storing the statistic name,
  #                and "stat_value" storing the associated value. 
  #
  # Returns:
  #    data.table, with columns "samp_mean", "samp_var". 
  
  # Select rows and columns `samp_dt` required for computing metrics. 
  if(subset_samp) {
    samp_dt <- select_mcmc_samp(samp_dt, burn_in_start=burn_in_start, test_labels=test_labels, 
                                param_types=param_types, param_names=param_names)
  }
  
  # Compute statistics. 
  mcmc_param_stats <- samp_dt[, .(mean=mean(sample), var=var(sample)), by=.(test_label, param_type, param_name)]
  
  # Convert to long format, if requested. 
  if(format_long) {
    mcmc_param_stats <- melt.data.table(mcmc_param_stats, id.vars=c("test_label", "param_type", "param_name"), 
                                        variable.name="stat_name", value.name="stat_value")
  } else {
    setnames(mcmc_param_stats, c("mean", "var"), c("stat_mean", "stat_var"))
  }
  
  return(mcmc_param_stats)
  
}


compute_mcmc_comparison_metrics <- function(samp_dt, test_label_1, test_label_2, metrics, burn_in_start=NULL,
                                            param_types=NULL, param_names=NULL) {
  # Computes metrics that quantify the difference between two distributions, given samples from the two distributions. 
  # The samples from the distributions are assumed to both be stored in `samp_dt`, with the column "test_label" 
  # identifying the two distributions. The distribution associated with the label `test_label_1` is assumed to be 
  # the "reference" distribution, so relative metrics are computed with respect to this distribution. This function 
  # computes both 1.) "individual metrics", which compare samples associated with each parameter one at a time, and 
  # 2.) "aggregate metrics", which take into account information across different parameters; e.g. computing a 
  # covariance matrix across a set of parameters. 
  #
  # Args:
  #    samp_dt: data.table, must be of the format described in `format_mcmc_output()`. 
  #    test_label_1: character(1), the label identifying the reference distribution. 
  #    test_label_2: character(1), the label identifying the second distribution. 
  #    metrics: character, vector of metrics. Currently only supports "mean" and "cov". 
  #    burn_in_start, param_types, param_names: passed to `select_mcmc_samp` to subset `samp_dt` to determine which
  #                                             parameters and samples will be included in the metric computations. 
  #
  # Returns: 
  #    list, with names "metrics_individual" and "metrics_agg". Each element is a data.table storing the metrics. 
  
  # Select rows and columns `samp_dt` required for computing metrics. 
  test_labels <- c(test_label_1, test_label_2)
  samp_dt_subset <- select_mcmc_samp(samp_dt, burn_in_start=burn_in_start, test_labels=test_labels, 
                                     param_types=param_types, param_names=param_names)
  
  # Compute univariate MCMC means and variance estimates.  
  mcmc_param_stats <- compute_mcmc_param_stats(samp_dt_subset, burn_in_start=burn_in_start, test_labels=test_labels, 
                                               param_types=param_types, param_names=param_names, subset_samp=FALSE)
  
  
  # Create data.table for storing aggregate metrics (i.e. metrics involving one or more parameters). 
  metrics_agg <- data.table(param_type = character(), 
                            metric = character(), 
                            test_lab_1 = character(), 
                            test_lab_2 = character(), 
                            diff = character(), 
                            diff_rel = character())
  setnames(metrics_agg, c("test_lab_1", "test_lab_2"), test_labels)
  
  # L2 distance.
  L2 <- function(x) sqrt(sum(x^2))
  
  if("mean" %in% metrics) {
    
    # Compare sample means parameter-by-parameter.
    means <- dcast(mcmc_param_stats, param_type+param_name ~ test_label, value.var = "samp_mean")
    means$diff <- abs(means[[test_label_1]] - means[[test_label_2]])
    means$mean_diff_rel <- means$diff / abs(means[[test_label_1]])
    means$metric <- "mean_abs_diff"
    
    # L2 distance between mean vectors (over all calibration parameters).
    mean_vec_diffs <- means[ , lapply(.SD, L2), .SDcols = c("diff", test_label_1, test_label_2), by=param_type]
    print(mean_vec_diffs)
    mean_vec_diffs[, diff_rel := diff / get(test_label_1)]
    mean_vec_diffs[, metric := "mean_L2"]
    metrics_agg <- rbindlist(list(metrics_agg, mean_vec_diffs), use.names = TRUE)
    
  }
  
  if("cov" %in% metrics) {
    
    for(pt in samp_dt_subset[, unique(param_type)]) {
      samp_wide_test_1 <- dcast(samp_dt_subset[(test_label == test_label_1) & (param_type == pt)], 
                                formula = itr~param_name, value.var = "sample")[, .SD, .SDcols = !"itr"]
      samp_wide_test_2 <- dcast(samp_dt_subset[(test_label == test_label_2) & (param_type == pt)], 
                                formula = itr~param_name, value.var = "sample")[, .SD, .SDcols = !"itr"]
      C1 <- cov(samp_wide_test_1)
      C2 <- cov(samp_wide_test_2)
      C1_frobenius <- L2(C1)
      
      cov_diffs <- data.table(param_type = pt, metric = "cov_frobenius")
      cov_diffs[, (test_label_1) := C1_frobenius]
      cov_diffs[, (test_label_2) := L2(C2)]
      cov_diffs[, diff := L2(C1 - C2)]
      cov_diffs[, diff_rel := diff / C1_frobenius]
      metrics_agg <- rbindlist(list(metrics_agg, cov_diffs), use.names=TRUE)
    }
    
  }
  
  return(list(metrics_individual=means, metrics_agg=metrics_agg))
  
}


compute_mcmc_running_err_multivariate <- function(samp_dt, mean_true, cov_true, 
                                                  param_type, test_labels=NULL,
                                                  param_names=NULL, burn_in_start=NULL, 
                                                  init_running_err_using_burnin=TRUE) {
  # A wrapper around `compute_samp_running_err_multivariate()` which computes the running mean
  # and covariance error (with respect to the baseline provided in `true_mean_list` and 
  # `true_cov_list`) across different test labels (for a single param_type). One call to
  # `compute_samp_running_err_multivariate()` is made for each test label. By default, 
  # the error measures are computed using all parameters within the specified parameter
  # type; however, `param_names` can also be specified to select a subset of parameters
  # within the parameter type; the multivariate error measures will then be computed for 
  # this subset. 
  #
  # Details for `init_running_err_using_burnin`: 
  # If no burn-in is specified (`burn_in_start` is NULL) then the running mean/covariance 
  # caclulations will start from the third iteration, with the first two iterations used 
  # to initialize the empirical mean and covariance estimates. If a burn-in is specified, 
  # then there are two possible behaviors. 1.) If `init_running_err_using_burnin` is 
  # TRUE, then the empirical mean and covariance will be initialized by computing the
  # mean and covariance estimates from the burn-in. The running mean and covariance 
  # estimates will be iteratively updated for the non-burn-in iterations. 2.) if 
  # init_running_err_using_burnin` is FALSE, then the burn-in will be dropped and 
  # the remaining iterations will be treated as in the NULL `burn_in_start` case 
  # (the mean and covariance estimates are initialized using the first two non-burn-in
  # iterations). However, the original iteration labels will still be used in the 
  # returned data.table. 
  #
  # Args:
  #    samp_dt: data.table, must be of the format described in `format_mcmc_output()`.    
  #    mean_true: numeric or matrix with one row, the true mean. If numeric vector, 
  #               must have  names set to parameter names; if matrix, must have column 
  #               names set to parameter names. If `param_names` is provided, `mean_true`
  #               will be subsetted to select only the relevant parameters. 
  #    cov_true: matrix of dimension (ncol(samp), ncol(samp)). The true covariance matrix. 
  #              Must have row and column names set to parameter names. If `param_names` is 
  #             provided, `cov_true` will be subsetted to select only the relevant parameters. 
  #    param_type: character(1), the parameter type to select. 
  #    test_labels: character, vector of test labels. If NULL uses all test labels in `samp_dt`. 
  #    param_names: character, if non-NULL, selects a subset of parameters within the param 
  #                 type. Also used to order parameters. 
  #    burn_in_start: integer, the burn-in vector; see remove_mcmc_samp_burnin() for details. 
  #                   See above for details on how the burn-in specification affects the 
  #                   running error calculations. 
  #    init_running_err_using_burnin: logical(1), whether or not to use the burn-in iterations 
  #                                   to initialize the running mean and covariance estimates. 
  #                                   See above description for details. 
  #
  # Returns:
  #    data.table, with columns "test_label", "param_type", "itr", "mean_err", "cov_err". See
  #    `compute_samp_running_err_multivariate()` for details on the error measures. Note that
  #    these error measures are computed with respect to the parameter vectors determined by 
  #    `param_names`. 
  
  err_dt <- data.table(test_label=character(), 
                       itr=integer(),
                       mean_err=numeric(),
                       cov_err=numeric())
  
  if(is.null(test_labels)) test_labels <- unique(samp_dt$test_label)
  
  mean_true <- drop(mean_true)
  if(!is.null(param_names)) {
    mean_true <- mean_true[param_names]
    cov_true <- cov_true[param_names, param_names]
  }
  
  # Don't exclude burn-in yet. 
  samp_dt <- select_mcmc_samp(samp_dt, test_labels=test_labels, param_types=param_type, 
                              param_names=param_names)
  
  for(lbl in test_labels) {
    samp_mat <- select_mcmc_samp_mat(samp_dt, test_label=lbl, param_type=param_type,
                                     burn_in_start=burn_in_start)
    
    burn_in_required <- isTRUE(burn_in_start[lbl] > 1) || 
      ((length(burn_in_start)==1) && burn_in_start > 1)
    if(init_running_err_using_burnin && burn_in_required) {
      samp_burnin <- select_mcmc_samp_mat(samp_dt, test_label=lbl, param_type=param_type,
                                          burn_in_start=burn_in_start, return_burnin=TRUE)
      mean_init <- colMeans(samp_burnin)
      cov_init <- cov(samp_burnin)
    } else {
      mean_init <- NULL
      cov_init <- NULL
    }
    
    err_list <- compute_samp_running_err_multivariate(samp_mat, mean_true, cov_true, mean_init, cov_init)
    dt_curr <- data.table(test_label=lbl, itr=err_list$itr, mean_err=err_list$mean, cov_err=err_list$cov)
    err_dt <- rbindlist(list(err_dt, dt_curr), use.names=TRUE)               
  }
  
  err_dt[, param_type := param_type]
  return(err_dt)
  
}


compute_samp_running_err_multivariate <- function(samp, mean_true, cov_true, mean_curr=NULL, cov_curr=NULL) {
  # Computes 1.) running L2 norm between running sample mean and a given true mean and 
  # 2.) running Frobenius norm between running empirical covariance and a given true 
  # covariance matrix. 
  #
  # Args:
  #    samp: matrix, with rows equal to samples, i.e. must be in the wide 
  #          format as returned by `select_mcmc_samp_mat()`. Must have column names set to the 
  #          parameter names and rownames set to the iteration numbers. 
  #    mean_true: numeric or matrix with one row, the true mean. Dimension must agree
  #               with the number of columns of `samp`. If numeric vector, must have 
  #               names set to parameter names; if matrix, must have column names set to
  #               parameter names. 
  #    cov_true: matrix of dimension (ncol(samp), ncol(samp)). The true covariance matrix. 
  #              Must have row and column names set to parameter names. 
  #    mean_curr: numeric or matrix with one row. Typically used when `samp` is excluding 
  #               some earlier samples and one wants to intialize the running mean to the 
  #               sample mean of these excluded samples. 
  #    cov_curr: matrix of dimension (ncol(samp), ncol(samp)), analoglous to `mean_curr`. 
  #              This is especially useful when early values of `samp` are excluded  
  #              in order to prevent a singular covariance matrix (e.g. when the early 
  #              iterations of MCMC are stuck at a single value). 
  #
  # Returns:
  #    list, with named elements "mean", "cov", and "itr". Each element is a vector and 
  #    all are of the same length. The first two are numeric vectors storing the mean 
  #    and covariance errors. The third is an integer vector storing the iteration number
  #    associated with each respective error. If `mean_curr` and `cov_curr` are provided 
  #    then these vectors will have length equal to `nrow(samp)`. If these arguments 
  #    are NULL, then the vectors will have length equal to `nrow(samp)-2`, as the first 
  #    2 iterations will be used to initialize the empirical covariance estimate. 
  
  N_samp <- nrow(samp)
  mean_err <- vector(mode="numeric", length=N_samp) 
  cov_err <- vector(mode="numeric", length=N_samp)
  mean_true <- drop(mean_true)
  if(is.null(names(mean_true))) stop("`mean_true` lacking parameter names in names attribute.")
  if(is.null(colnames(cov_true))) stop("`cov_true` lacking parameter names in colnames attribute.")
  if(xor(is.null(mean_curr), is.null(cov_curr))) stop("`mean_curr` and `cov_curr` must be both NULL or both non-NULL.")
  
  # Ensure proper ordering. 
  params_ordered <- names(mean_true)
  cov_true <- cov_true[params_ordered, params_ordered, drop=FALSE]
  samp <- samp[, params_ordered, drop=FALSE]
  
  # If not provided, initialize mean and cov using first two samples. 
  if(is.null(cov_curr)) {
    mean_curr <- colMeans(samp[1:2,,drop=FALSE])
    cov_curr <- cov(samp[1:2,,drop=FALSE])
    idx_start <- 3
    itrs <- seq(1, N_samp)
    itr <- itrs[idx_start]
  } else {
    mean_curr <- drop(mean_curr)[params_ordered]
    cov_curr <- cov_curr[params_ordered, params_ordered]
    idx_start <- 1
    itrs <- as.integer(rownames(samp))
    itr <- itrs[1]
  }
  
  for(i in seq(idx_start, nrow(samp))) {
    # Update mean and covariance. 
    cov_curr <- tcrossprod(samp[i,]-mean_curr)/itr + (itr-2)/(itr-1)*cov_curr
    mean_curr <- samp[i,]/itr + ((itr-1)/itr)*mean_curr
    
    # Compute error measures. 
    mean_err[i] <- sqrt(sum((mean_curr - mean_true)^2))
    cov_err[i] <- sqrt(sum((cov_curr - cov_true)^2))
    
    # Update iteration. 
    itr <- itr + 1
  }    
  
  return(list(mean=mean_err[idx_start:N_samp], cov=cov_err[idx_start:N_samp], 
              itr=as.integer(rownames(samp))[idx_start:N_samp]))
  
}


# ------------------------------------------------------------------------------
# MCMC Plotting Functions. 
# ------------------------------------------------------------------------------


get_hist_plot <- function(samples_list, col_sel=1, bins=30, vertical_line=NULL, xlab="samples", ylab="density", 
                          main_title="Histogram", data_names=NULL) {
  # Generates a single plot with one or more histograms. The input data `samples_list` must be a list of matrices. 
  # For example, the first element of the list could be a matrix of samples from a reference distribution, and the 
  # second element a matrix of samples from an approximate distribution. This function will select one column from 
  # each matrix based on the index `col_sel`. Each selected column will be turned into a histogram. 
  # Note that the matrices may have different lengths (numbers of samples).
  #
  # Args:
  #    samples_list: list of matrices, each of dimension (num samples, num parameters). 
  #    col_sel: integer(1), the column index to select from each matrix. i.e. this selects a particular parameter, whose
  #             samples will be transformed into a histogram. 
  #    bins: integer(1), number of bins to use in histogram. 
  #    vertical_line: If not NULL, the x-intercept to be used for a plotted vertical line. 
  #    xlab, ylab, main_title: x and y axis labels and plot title. 
  #    data_names: character(), vector equal to the length of the list `samples_list`. These are the names used in 
  #                the legend to identify the different histograms. 
  #
  # Returns:
  #    ggplot2 object. 
  
  # Pad with NAs in the case that the number of samples is different across different parameters. 
  N_samp <- sapply(samples_list, nrow)
  if(length(unique(N_samp)) > 1) {
    N_max <- max(N_samp)
    for(j in seq_along(samples_list)) {
      if(nrow(samples_list[[j]]) < N_max) {
        samples_list[[j]] <- rbind(samples_list[[j]], matrix(NA, nrow=N_max - nrow(samples_list[[j]]), ncol=ncol(samples_list[[j]])))
      }
    }
  }
  
  dt <- as.data.table(lapply(samples_list, function(mat) mat[,col_sel]))
  if(!is.null(data_names)) setnames(dt, colnames(dt), data_names)
  dt <- melt(dt, measure.vars = colnames(dt), na.rm=TRUE)
  
  plt <- ggplot(data = dt, aes(x=value, color=variable)) + 
    geom_histogram(aes(y=..density..), bins=bins, fill="white", alpha=0.2, position="identity") + 
    xlab(xlab) + 
    ylab(ylab) + 
    ggtitle(main_title)
  
  if(!is.null(vertical_line)) {
    plt <- plt + geom_vline(xintercept=vertical_line, color="red")
  }
  
  return(plt)
  
}


get_trace_plots <- function(samp_dt, burn_in_start=NULL, test_labels=NULL, param_types=NULL, 
                            param_names=NULL, save_dir=NULL) {
  # Operates on a data.table of MCMC samples in the long format, as returned by `format_mcmc_output()`. Generates one 
  # MCMC trace plot per valid `param_name`-`param_type`-`test_label` combination. 
  #
  # Args:
  #    samp_dt: data.table of MCMC samples, in long format as returned by format_mcmc_output()`. 
  #    burn_in_start: If NULL, selects all MCMC iterations. If integer of length 1, this is interpreted as the starting 
  #                   iteration for all parameters - all earlier iterations are dropped. If vector of length > 1, must 
  #                   be a named vector with names set to valid test label values. This allows application of a different 
  #                   burn-in start iteration for different test labels. 
  #    test_labels, param_types, param_names: vectors of values to include in selection corresponding to columns 
  #                                           "test_label", "param_type", and "param_name" in `samp_dt`. A NULL  
  #                                            value includes all values found in `samp_dt`. 
  #    save_dir: character(1), if not NULL, a file path to save the plots to. 
  #
  # Returns:
  #    list of ggplot objects, one for each `param_name`-`param_type`-`test_label` combination. 
  
  # Determine which plots to create by subsetting rows of `samp_dt`. 
  samp_dt_subset <- select_mcmc_samp(samp_dt, burn_in_start=burn_in_start, test_labels=test_labels, 
                                     param_types=param_types, param_names=param_names)
  plt_id_vars <- unique(samp_dt_subset[, .(test_label, param_type, param_name)])
  
  # Generate plots. 
  plts <- list()
  for(j in 1:nrow(plt_id_vars)) {
    test_label_curr <- plt_id_vars[j, test_label]
    param_type_curr <- plt_id_vars[j, param_type]
    param_name_curr <- plt_id_vars[j, param_name]
    plt_label <- paste(test_label_curr, param_type_curr, param_name_curr, sep="_")
    
    plts[[plt_label]] <- ggplot(data = samp_dt_subset[(test_label == test_label_curr) & 
                                                        (param_type == param_type_curr) & 
                                                        (param_name == param_name_curr)], aes(x=itr, y=sample)) + 
      geom_line() + 
      ggtitle(paste0("Trace Plot: ", plt_label)) + 
      xlab("Iteration")
  }
  
  if(!is.null(save_dir)) save_plots(plts, "trace", save_dir)
  
  return(plts)
  
}


get_1d_kde_plot_comparisons <- function(samp_dt, N_kde_pts=100, burn_in_start=NULL, test_labels=NULL, 
                                        param_types=NULL, param_names=NULL, test_label_baseline=NULL,
                                        xlab="parameter", ylab="kde", save_dir=NULL) {
  # Produces one plot per unique param type-param name combination. Each plot contains one line
  # per `test_label` that has samples for that specific parameter. Each line is produced from 
  # a 1d kernel density estimate (KDE) constructed fromn the MCMC samples. N_kde_pts` is the
  # `number of points at which the KDE approximation is evaluated for each univariate variable.

  # Determine which plots to create by subsetting rows of `samp_dt`. 
  if(!is.null(test_label_baseline) && !is.null(test_labels) && !(test_label_baseline %in% test_labels)) {
    test_labels <- c(test_labels, test_label_baseline)
  }
  samp_dt_subset <- select_mcmc_samp(samp_dt, burn_in_start=burn_in_start, test_labels=test_labels, 
                                     param_types=param_types, param_names=param_names)
  
  # One plot is generated for each unique `param_type`-`param_name` combination. 
  plt_id_vars <- unique(samp_dt_subset[, .(param_type, param_name)])
  
  # Separate out data to be used as the baseline for comparison in each plot, if provided. 
  if(!is.null(test_label_baseline)) {
    samp_dt_baseline <- samp_dt_subset[test_label==test_label_baseline]
    samp_dt_subset <- samp_dt_subset[test_label != test_label_baseline]
  }
  
  # Generate plots. 
  plts <- list()
  for(j in 1:nrow(plt_id_vars)) {
    
    # Define the parameter being plotted. One KDE will be computed for each test label
    # that has samples for this parameter. 
    param_type_curr <- plt_id_vars[j, param_type]
    param_name_curr <- plt_id_vars[j, param_name]
    plt_label <- paste(param_type_curr, param_name_curr, sep="_")
    samp_dt_param <- samp_dt_subset[(param_type == param_type_curr) & 
                                    (param_name == param_name_curr), .(test_label, sample)] 
    test_labels_curr <- unique(samp_dt_param$test_label)
    
    # Determine the grid of points at which the KDE will be evaluated. 
    bound_lower <- quantile(samp_dt_param$sample, 0.025)
    bound_upper <- quantile(samp_dt_param$sample, 0.975)
    kde_pts <- seq(bound_lower, bound_upper, length.out=N_kde_pts)
    
    # Loop over test labels, constructing KDE for each label. 
    kde_mat <- matrix(nrow=N_kde_pts, ncol=length(test_labels_curr), 
                      dimnames=list(NULL, test_labels_curr))
    for(lbl in test_labels_curr) {
      kde_fit <- kde1d(samp_dt_subset[test_label==lbl, sample])
      kde_mat[,lbl] <- dkde1d(kde_pts, kde_fit)
    }

    # Construct KDE for the baseline label, if provided. 
    kde_baseline <- NULL
    if(!is.null(test_label_baseline)) {
      samp_baseline_param <- samp_dt_baseline[(param_type == param_type_curr) & 
                                              (param_name == param_name_curr), sample] 
      kde_fit <- kde1d(samp_baseline_param)
      kde_baseline <- dkde1d(kde_pts, kde_fit)
    }
    
    # Construct KDE comparison plot for the current parameter. 
    plt_curr <- plot_curves_1d_helper(kde_pts, kde_mat, y_new=kde_baseline,
                                      plot_title=plt_label, xlab=plt_label, ylab="kde")
    plts[[plt_label]] <- plt_curr
  }
  
  if(!is.null(save_dir)) save_plots(plts, "kde1d", save_dir)
  
  return(plts)
  
}


get_mcmc_moments_scatter_plot_comparisons <- function(samp_dt, test_label_baseline, burn_in_start=NULL,
                                                      test_labels=NULL, param_types=NULL, param_names=NULL,
                                                      xlab="observed", ylab="predicted", save_dir=NULL) {
  # This function currently produces one plot per unique (param type, stat name) combination. 
  # Currently the stats are hard-coded to be mean and standard deviation, though this can easily 
  # be generlized to allow the user to pass in any number of statistics of interest. 
  # The plots are scatter plots summarizing the mean/standard deviation MCMC estimates to some 
  # "baseline" estimates. Currently, the baseline estimates are also computed from samples, using 
  # the test label `test_label_baseline`, which must be present in `samp_dt` (this can be generalized
  # to allow the baseline estimates to be explicitly passed if desired; e.g., if they are known exactly 
  # as is the case in linear Gaussian problems). Each plot contains a point for each (param name, test label)
  # combination, where the param names are those that belong to the param type being considered in 
  # that plot. The param names are differentiated by shape and the test labels are differentiated by color. 
  
  # Determine which plots to create by subsetting rows of `samp_dt`. 
  if(!is.null(test_labels) && !(test_label_baseline %in% test_labels)) {
    test_labels <- c(test_labels, test_label_baseline)
  }
  samp_dt_subset <- select_mcmc_samp(samp_dt, burn_in_start=burn_in_start, test_labels=test_labels, 
                                     param_types=param_types, param_names=param_names)
  plt_id_vars <- unique(samp_dt_subset[, .(test_label, param_type, param_name)])
  
  # Compute statistics from posterior samples. 
  samp_stats <- compute_mcmc_param_stats(samp_dt_subset, subset_samp=FALSE, format_long=TRUE) 
  samp_stats[stat_name=="var", `:=`(stat_name="sd", stat_value=sqrt(stat_value))]

  # Separate out data to be used as the baseline for comparison in each plot. 
  samp_stats_baseline <- samp_stats[test_label==test_label_baseline]
  samp_stats <- samp_stats[test_label != test_label_baseline]
  samp_stats <- merge(samp_stats, samp_stats_baseline, all.x=TRUE, 
                      by=c("param_type", "param_name", "stat_name"), suffixes=c("", "_baseline"))
                      
  # Produce one plot per unique stat_name, param_type, or (stat_name, param_type, param_name) combination. 
  plt_id_vars <- unique(samp_stats[, .(param_type, stat_name)])
  plt_list <- list()
  for(i in 1:nrow(plt_id_vars)) {
    param_type_curr <- plt_id_vars[i, param_type]
    stat_name_curr <- plt_id_vars[i, stat_name]
    plt_lbl <- paste(param_type_curr, stat_name_curr, sep="-")
    samp_stats_curr <- samp_stats[(param_type==param_type_curr) & (stat_name==stat_name_curr)]
    
    plt <- ggplot(samp_stats_curr) + 
            geom_point(aes(x=stat_value_baseline, y=stat_value, color=test_label, shape=param_name)) + 
            geom_abline(slope=1, intercept=0, color="red") + 
            ggtitle(plt_lbl) + xlab(paste0(plt_lbl, ", ", test_label_baseline)) + 
            ylab(plt_lbl)
    plt_list[[plt_lbl]] <- plt
  }
  
  if(!is.null(save_dir)) save_plots(plts, "kde1d", save_dir)
  
  return(plt_list)
  
}


get_hist_plot_comparisons <- function(samp_dt, burn_in_start=NULL, test_labels=NULL, param_types=NULL, param_names=NULL,
                                      test_label_baseline=NULL, xlab="samples", ylab="density", bins=30, save_dir=NULL) {
  # Operates on a data.table of MCMC samples in the long format, as returned by `format_mcmc_output()`. Generates one 
  # MCMC marginal histogram plot per valid `param_name`-`param_type`-`test_label` combination. If `test_label_baseline`
  # is non-NULL, then each plot will include a second histogram corresponding to the specified test label; this is 
  # useful if wanting to compare approximate samples against some sort of baseline, for example. 
  #
  # Args:
  #    samp_dt: data.table of MCMC samples, in long format as returned by format_mcmc_output()`. 
  #    burn_in_start: If NULL, selects all MCMC iterations. If integer of length 1, this is interpreted as the starting 
  #                   iteration for all parameters - all earlier iterations are dropped. If vector of length > 1, must 
  #                   be a named vector with names set to valid test label values. This allows application of a different 
  #                   burn-in start iteration for different test labels. 
  #    test_labels, param_types, param_names: vectors of values to include in selection corresponding to columns 
  #                                           "test_label", "param_type", and "param_name" in `samp_dt`. A NULL  
  #                                            value includes all values found in `samp_dt`. 
  #    test_label_baseline: character(1) or NULL. If non-NULL, then must be a valid test label with associated 
  #                         samples in `samp_dt`. In this case, the histograms corresponding to this baseline 
  #                         test label will be overlaid on the plots for all other test labels. 
  #    xlab, ylab, bins: ggplot arguments, all passed to `get_hist_plot()`. 
  #    save_dir: character(1), if not NULL, a file path to save the plots to. 
  #
  # Returns: 
  #    list, each element being a ggplot object. 
  
  # Determine which plots to create by subsetting rows of `samp_dt`. 
  if(!is.null(test_label_baseline) && !is.null(test_labels) && !(test_label_baseline %in% test_labels)) {
    test_labels <- c(test_labels, test_label_baseline)
  }
  samp_dt_subset <- select_mcmc_samp(samp_dt, burn_in_start=burn_in_start, test_labels=test_labels, 
                                     param_types=param_types, param_names=param_names)
  plt_id_vars <- unique(samp_dt_subset[, .(test_label, param_type, param_name)])
  
  # Separate out data to be used as the baseline for comparison in each plot, if provided. 
  if(!is.null(test_label_baseline)) {
    samp_dt_baseline <- samp_dt_subset[test_label==test_label_baseline]
    plt_id_vars <- plt_id_vars[test_label != test_label_baseline]
  }
  
  # Generate plots. 
  plts <- list()
  for(j in 1:nrow(plt_id_vars)) {
    test_label_curr <- plt_id_vars[j, test_label]
    param_type_curr <- plt_id_vars[j, param_type]
    param_name_curr <- plt_id_vars[j, param_name]
    plt_label <- paste(test_label_curr, param_type_curr, param_name_curr, sep="_")
    samp <- samp_dt_subset[(test_label == test_label_curr) & 
                             (param_type == param_type_curr) & 
                             (param_name == param_name_curr), sample]
    samp_list <- list()
    samp_list[[1]] <- matrix(samp, ncol=1)
    data_names <- test_label_curr
    
    if(!is.null(test_label_baseline)) {
      samp_baseline <- samp_dt_baseline[(param_type == param_type_curr) & 
                                          (param_name == param_name_curr), sample] 
      
      samp_list[[2]] <- matrix(samp_baseline, ncol=1)
      data_names <- c(data_names, test_label_baseline)
    }
    
    plts[[plt_label]] <- get_hist_plot(samp_list, bins=bins, xlab=param_name_curr, ylab="density", 
                                       main_title=test_label_curr, data_names=data_names) 
    
  }
  
  if(!is.null(save_dir)) save_plots(plts, "hist", save_dir)
  return(plts)
  
}


get_1d_coverage_plots <- function(samp_dt, test_label_baseline, burn_in_start=NULL, test_labels=NULL,
                                  param_types=NULL, param_names=NULL, xlab="observed", ylab="predicted", 
                                  save_dir=NULL, probs=seq(0.5, 1.0, .1)) {
  # Produces one plot per unique (param type, param name) combination. Each plot will contain one line 
  # per test label that has samples associated with the specific parameter being plotted. These lines 
  # summarize the coverage of the distributions (i.e., the samples for each test label) with respect 
  # to some "baseline" test label, specified by `test_label_baseline` (this typically represents 
  # some notion of the "true" distribution). For a specific plot, let us consider how a single 
  # point is computed: 
  #    (1) A highest posterior density (HPD) interval of probability p is estimated from samples 
  #        (using the HDInterval package) for each test label. 
  #    (2) Let [a,b] be the bounds of this estimated interval for a specific test label. Next we 
  #        estimate the probability of the baseline distribution falling within the interval [a,b]. 
  #        This is accomplished via a 1d kernel density estimate (KDE) of the baseline distribution 
  #        using the package kde1d. Suppose we estimate this probability to be q. 
  #    (3) The point (p,q) is then plotted. This is repeated for each test label, and also repeated
  #        for a set of different probabilities p, specified by the arguments `probs`. The points  
  #        corresponding to the same test label are connected with interpolating lines. 
                                  
  # Determine which plots to create by subsetting rows of `samp_dt`. 
  if(!is.null(test_labels) && !(test_label_baseline %in% test_labels)) {
    test_labels <- c(test_labels, test_label_baseline)
  }
  samp_dt_subset <- select_mcmc_samp(samp_dt, burn_in_start=burn_in_start, test_labels=test_labels, 
                                     param_types=param_types, param_names=param_names)
  plt_id_vars <- unique(samp_dt_subset[, .(param_type, param_name)])
  
  # Separate out data to be used as the baseline for comparison in each plot.
  samp_dt_baseline <- samp_dt_subset[test_label==test_label_baseline]
  samp_dt_subset <- samp_dt_subset[test_label != test_label_baseline]
  
  # Create one plot per unique (param_type, param_name) combination. 
  plt_list <- list()
  for(i in 1:nrow(plt_id_vars)) {
    # (param_type, param_name) combination for current plot. 
    param_type_curr <- plt_id_vars[i, param_type]
    param_name_curr <- plt_id_vars[i, param_name]
    samp_dt_param <- samp_dt_subset[(param_type==param_type_curr) & (param_name==param_name_curr)]
    lbl_curr <- paste(param_type_curr, param_name_curr, sep="-")
    test_labels_curr <- unique(samp_dt_param$test_label)
    
    # Estimate highest posterior density interval for each test label and probability. 
    dt_plt <- data.table(test_label=character(), prob=numeric(), hdi_interval_lower=numeric(),
                         hdi_interval_upper=numeric())
    for(lbl in test_labels_curr) {
      compute_hdi <- function(p) HDInterval::hdi(samp_dt_param[test_label==lbl, sample], p)
      hdi_bounds <- mapply(compute_hdi, probs)
      dt_plt_lbl <- data.table(test_label=lbl, prob=probs, hdi_interval_lower=hdi_bounds["lower",], 
                               hdi_interval_upper=hdi_bounds["upper",])
      dt_plt <- rbindlist(list(dt_plt, dt_plt_lbl), use.names=TRUE)
    }
    
    # Compute KDE for the baseline distribution. 
    kde_fit <- kde1d(samp_dt_baseline[(param_type==param_type_curr) & (param_name==param_name_curr), sample])

    # Use KDE to estimate probability of each of the HPD intervals computed above. 
    dt_plt[, coverage_prob := pkde1d(hdi_interval_upper, kde_fit) - pkde1d(hdi_interval_lower, kde_fit)]
    
    # Generate plot. 
    plt <- ggplot(dt_plt) + 
            geom_point(aes(x=prob, y=coverage_prob, color=test_label)) + 
            geom_line(aes(x=prob, y=coverage_prob, color=test_label)) + 
            geom_abline(slope=1, intercept=0, color="red", linetype="dashed") + 
            ggtitle(paste0("Coverage: ", lbl_curr)) + xlab("Prob") + ylab("Coverage Prob")
        
    plt_list[[lbl_curr]] <- plt
  }
  
  # Optionally save plots to file. Return plot list. 
  if(!is.null(save_dir)) save_plots(plt_list, "coverage", save_dir)
  return(list(plots=plt_list, plot_data=dt_plt))

}


get_2d_density_contour_plot <- function(samples_list, col_sel=c(1,2), xlab="theta1", ylab="theta2", main_titles=NULL) {
  # Plots the contours of a 2D kernel density estimate. If `samples_list` contains multiple elements, then one plot will be 
  # returned per element. Each element of `samples_list` is matrix of dimension N_samples x N_params. The vector `col_sel`
  # determines which 2 columns will be used for the 2D KDE plot in each matrix. 
  #
  # Args:
  #    samples_list: list of matrices, each of dimension (num samples, 2). 
  #    col_sel: integer(1), the two column indees to select from each matrix. i.e. this selects two particular parameters, whose
  #             samples will be used to compute the KDE. 
  #    xlab, ylab, main_title: x and y axis labels and plot title. 
  #
  # Returns:
  #    list, of equal length to `samples_list`. Each element is a ggplot object. 
  
  if(is.null(main_titles)) main_titles <- paste0("2D KDE Countours: ", seq_along(samples_list))
  plts <- vector(mode = "list", length = length(samples_list))
  
  for(j in seq_along(samples_list)) {
    df <- data.frame(samples_list[[j]])
    x <- colnames(df)[col_sel[1]]
    y <- colnames(df)[col_sel[2]]
    plts[[j]] <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) + 
      geom_density_2d_filled() + 
      xlab(xlab) + 
      ylab(ylab) + 
      ggtitle(main_titles[j])
  }
  
  return(plts)
  
}


get_overlaid_2d_density_contour_plot <- function(samp_baseline, samp_overlay, col_sel=c(1,2), main_title=NULL) {
  # Plots the contours of 2D kernel density estimates for samples from two different distributions over 2D input spaces 
  # so that they may compared. The samples `samp_baseline` will be plotted with a "filled" KDE heatmap, while the
  # the contours for `samp_overlay` will be overlaid on top and not filled. The vector `col_sel`
  # determines which 2 columns of both matrices will be used for the 2D KDE plot. 
  #
  # Args:
  #    samp_baseline: matrix, of dimension (num samples, 2). 
  #    samp_overlay: matrix, of dimension (num samples, 2).
  #    col_sel: integer or character. If integer, the two column indices to select from each matrix. i.e. this selects two particular parameters, whose
  #             samples will be used to compute the KDE. If character, the column names to select. The column names of teh selected columns 
  #             from the two matrices must align. 
  #    main_title: plot title. 
  #
  # Returns:
  #    ggplot object.  
  
  # Select columns. 
  df_baseline <- data.frame(samp_baseline[, col_sel])
  df_overlay <- data.frame(samp_overlay[, col_sel])
  if(!all(colnames(df_baseline) == colnames(df_overlay))) stop("Column names of selected cols in `samp_baseline` and `samp_overlay` do not match.")
  colnames(df_baseline) <- colnames(df_overlay) <-  c("theta1", "theta2")
  
  # Plot title. 
  if(is.null(main_title)) main_title <- "2D KDE Countours"
  
  # Generate plot.                      
  plt <- ggplot(df_baseline, mapping = aes(x = theta1, y = theta2)) + 
    geom_density_2d_filled() + 
    geom_density_2d(df_overlay, mapping = aes(x = theta1, y = theta2), color = "white") + 
    xlab(colnames(df_baseline)[1]) + 
    ylab(colnames(df_baseline)[2]) + 
    ggtitle(main_title)
  
  return(plt)
  
}


get_2d_heatmap_plot <- function(X, y, param_names, samples_kde=NULL, samples_points=NULL,  
                                raster=FALSE, point_coords=NULL, main_title="Heatmap", 
                                bigger_is_better=TRUE, legend_label="y", log_scale=FALSE, 
                                point_coords_shape=8, point_coords_col="black", 
                                samples_points_size=1, point_coords_size=3, 
                                samples_kde_lab="KDE", samples_points_lab="samples_points", KDE_opacity=1.0) {
  # Plots a 2d heatmap of a scalar quantity `y`. Optionally overlays contours of a 2d 
  # kernel density estimate from `samples`. The input locations are given by the 
  # M x 2 matrix `X`. If these input locations correspond to an evenly-spaced grid, 
  # then `raster = TRUE` may be set to create a classic heatmap produced over a dense 
  # grid. Altenatively, if `X` consists of more sparsely sampled or non-evenly-spaced 
  # locations, `raster = FALSE` should be set and the resulting plot will plot the 
  # individual points, which will still be colored in heatmap fashion. 
  #
  # Args:
  #    X: matrix, of dimension M x D with colnames specified. `param_names` will be used to 
  #       select the two input dimensions of `X`. These define the input locations used in the 
  #       heatmap. 
  #    y: numeric(M), the scalar output value used to determine the colors in the heatmap. 
  #    samples_kde: matrix, with colnames specified. `param_names` will be used to 
  #                 select the two input dimensions of `samples_kde`. These input points 
  #                 will not be part of the heatmap. Instead, they will be used to construct 
  #                 a 2D KDE estimate and the contours of this KDE will be overlaid on the 
  #                 heatmap. 
  #    samples_points: matrix, with colnames specified. `param_names` will be used to 
  #                    select the two input dimensions of `samples_points`. These input points 
  #                    will be directly plotted as points on the plot and colored red. This 
  #                    argument is typically used to plot design points. 
  #    raster: logical(1), see above description. Set to TRUE when `X` is a dense grid of evenly-spaced  
  #            points FALSE when the points in `X` are not evenly-spaced or are sparse. 
  #    point_coords: numeric(2), coordinates to plot a single point as a red triangle. This typically 
  #                  corresponds to the location of the true parameter value. 
  #    main_title: character(1), the title of the plot. 
  #    bigger_is_better: logical(1), if TRUE interprets `y` such that larger values are "better" (i.e. higher 
  #                      posterior density). Set to FALSE for values like SSR, where smaller is better. 
  #    legend_label: character(1), the title for the legend which indicates the color scale. 
  
  
  if(length(param_names) != 2) stop("<param_names> must have length 2.")
  color_direction <- ifelse(bigger_is_better, 1, -1)
  color_breaks <- c()
  color_values <- c()
  
  df <- as.data.frame(cbind(X[, param_names], y))
  colnames(df) <- c("theta1", "theta2", "y")
  
  plt_transformation <- ifelse(log_scale, "log10", "identity")
  if(log_scale) {
    main_title <- paste0(main_title, ", log10 scale")
  }
  
  #
  # Heatmap. 
  #
  # NOTE: it is essential that in neither case a "color scale" is added, since this causes difficulties when trying to 
  # manually set a color scale later. To make this work in the non-raster case, I set "`hape = 21`, which is a point 
  # that has both fill and color attributes. We can then use the fill attribute as the mapping, while removing 
  # the border of these points with `stroke = NA`. 
  if(raster) {
    plt <- ggplot() + 
      geom_tile(data = df, aes(x = theta1, y = theta2, fill = y)) + 
      scale_fill_viridis(discrete=FALSE, direction = color_direction, trans = plt_transformation) + 
      labs(fill = legend_label)
  } else {
    plt <- ggplot() + 
      geom_point(data = df, aes(x = theta1, y = theta2, fill = y), shape = 21, stroke = NA) + # changed color to fill
      scale_fill_viridis(discrete=FALSE, direction = color_direction, trans = plt_transformation) + # changed scale_color_viridis to scale_fill_viridis
      labs(fill = legend_label) # changed color to fill
  }
  
  # Title and axis labels. 
  plt <- plt + ggtitle(main_title) + xlab(param_names[1]) + ylab(param_names[2])
  
  # Density contours from samples. 
  if(!is.null(samples_kde)) {
    if(!all(param_names %in% colnames(samples_kde))) stop("<param_names> must be column names of <samples_kde>.")
    
    samples_kde <- as.data.frame(samples_kde[, param_names])
    colnames(samples_kde) <- c("theta1", "theta2")
    
    plt <- plt + geom_density_2d(data = samples_kde, mapping = aes(x = theta1, y = theta2, color = samples_kde_lab), alpha = KDE_opacity)
    color_breaks <- c(color_breaks, samples_kde_lab)
    color_values <- c(color_values, setNames("blue", samples_kde_lab))
  }
  
  # Plot points. 
  if(!is.null(samples_points)) {
    if(!all(param_names %in% colnames(samples_points))) stop("<param_names> must be column names of <samples_points>.")
    samples_points <- as.data.frame(samples_points[, param_names])
    colnames(samples_points) <- c("theta1", "theta2")
    
    plt <- plt + geom_point(data = samples_points, mapping = aes(x = theta1, y = theta2, color = samples_points_lab), size = samples_points_size)
    color_breaks <- c(color_breaks, samples_points_lab)
    color_values <- c(color_values, setNames("red", samples_points_lab))
  }
  
  # Mark specific point in plot. 
  if(!is.null(point_coords)) {
    plt <- plt + geom_point(data = data.frame(theta1 = point_coords[1], theta2 = point_coords[2]), 
                            aes(x = theta1, y = theta2), color = point_coords_col, 
                            shape = point_coords_shape, size = point_coords_size)
  }
  
  # Legend. 
  if(length(color_breaks) > 0) {
    plt <- plt + scale_colour_manual(aesthetics = "colour", name = "", breaks = color_breaks, values = color_values)
  }
  
  return(plt)
  
}


get_2d_heatmap_plots <- function(X, Y, param_names, samples_kde=NULL, samples_points=NULL,  
                                 raster=FALSE, point_coords=NULL, main_title=NULL, 
                                 base_main_title="Heatmap", bigger_is_better=TRUE, legend_label="y") {
  # A wrapper function around `get_2d_headmap_plot()` that allows `Y` to be multivariate. 
  # Each column of `Y` is interpreted as the response variable for a different heatmap plot. 
  # Each column is fed to `get_2d_headmap_plot()` and a list of plots is returned. 
  #
  # Args:
  #    The same as `get_2d_headmap_plot()` except that `Y` is now a matrix, with each 
  #    column a variable to use as the response in a heatmap plot. `Y` should have column names, 
  #    which will be used to generate plot titles. The column names will be appended to 
  #    `base_main_title` to create the full titles. If `main_title` is specified, this 
  #    will be used as the plot titles, overwriting `base_main_title`. 
  #
  # Returns: 
  #    list, of length equal to the number of columns of `Y`. The names of the list will be 
  #    set to the corresponding column names of `Y`. Each element is the corresponding plot 
  #    returned from `get_2d_heatmap_plot()`. 
  
  if(is.null(colnames(Y))) output_variables <- paste("output", seq(1, ncol(Y)))
  else output_variables <- colnames(Y)
  
  plts <- vector(mode = "list", length = ncol(Y))
  for(j in seq_along(plts)) {
    
    plt_title <- ifelse(is.null(main_title), paste(base_main_title, output_variables[j], sep = ": "), main_title)
    
    plts[[j]] <- get_2d_heatmap_plot(X = X, 
                                     y = Y[,j], 
                                     param_names = param_names, 
                                     samples_kde = samples_kde, 
                                     samples_points = samples_points,  
                                     raster = raster, 
                                     point_coords = point_coords, 
                                     main_title = plt_title, 
                                     bigger_is_better = bigger_is_better, 
                                     legend_label = legend_label)
  }
  
  names(plts) <- output_variables
  
  return(plts)
  
}


get_2d_Bayes_opt_heatmap_plot <- function(theta_vals, computer_model_data, param_names, samples_kde = NULL, init_design_points = NULL,  
                                          sequential_design_points = NULL, raster = FALSE, point_coords = NULL,  
                                          main_title = "Heatmap: Sequential Design", bigger_is_better = TRUE, 
                                          legend_label = "y", log_scale = FALSE, SSR_vals = NULL, llik_vals = NULL, lprior_vals = NULL, 
                                          lpost_vals = NULL, sig2_eps = NULL, init_design_points_size = 1, sequential_design_points_size = 1) {
  # TODO: currently main title and legend label arguments do nothing. 
  
  # Log Posterior response surface plot. 
  plt <- get_2d_response_surface_plot_posterior(theta_vals = theta_vals, 
                                                computer_model_data = computer_model_data, 
                                                param_names = param_names, 
                                                output_variables = output_variables, 
                                                raster = raster, 
                                                point_coords = point_coords, 
                                                samples_kde = samples_kde, 
                                                samples_points = init_design_points, 
                                                SSR_vals = SSR_vals, 
                                                llik_vals = llik_vals,  
                                                lprior_vals = lprior_vals, 
                                                lpost_vals = lpost_vals, 
                                                theta_prior_params = theta_prior_params, 
                                                sig2_eps = sig2_eps, 
                                                main_title = main_title, 
                                                samples_points_size = init_design_points_size)
  
  # Add sequentially chosen design points. 
  if(!is.null(sequential_design_points)) {
    sequential_design_points <- as.data.frame(sequential_design_points[, param_names, drop=FALSE])
    colnames(sequential_design_points) <- c("theta1", "theta2")
    sequential_design_points[, "ID"] <- as.character(seq(1, nrow(sequential_design_points)))
    
    plt <- plt + geom_text(data = sequential_design_points, 
                           mapping = aes(x = theta1, y = theta2, label = ID), color = "red", size = sequential_design_points_size)
    
  }  
  
  return(plt)
  
}


get_2d_response_surface_plot <- function(computer_model_data, theta_vals, param_names, response_surface, 
                                         theta_prior_params = NULL, output_variables = NULL, 
                                         combine_outputs = TRUE, raster = FALSE, point_coords = NULL, 
                                         samples_kde = NULL, samples_points = NULL, scale_inputs = FALSE, 
                                         input_bounds = NULL, SSR_vals = NULL, llik_vals = NULL, 
                                         lprior_vals = NULL, lpost_vals = NULL, ...) {
  # A wrapper for `get_2d_heatmap_plot()` that specializes in plotting quantities of interest 
  # related to the calibration problem. This function can be used when the dimension of the 
  # calibration parameter space is 2, or to plot 2-dimensional projections when the input 
  # dimension is larger than 2. This function can plot a heatmap 
  # of 1.) the SSR (sum of squared residual) surface, 2.) the likelihood surface,  
  # 3.) the posterior surface, or 4.) the prior surface. It also allows the user to control 
  # whether individual plots are generated for each output variable, or if all the output variables 
  # should  be combined to plot the overall likelihood, posterior, etc. It also allows points 
  # to be overlaid on the plot, which typically represent design points or samples 
  # from a ground truth distribution. Additionally a single point can be added to 
  # mark the true value of the parameters. 
  #
  # Args:
  #    response_surface: character(1), either "SSR", "likelihood", "prior", or "posterior".
  #    combine_outputs: logical(1), currently only relevant if `response_surface` is "likelihood". 
  #                     If TRUE, returns a single plot of the overall log-likelihood. Otherwise, 
  #                     returns one plot per output variable, where each plot is a heatmap for the 
  #                     log likelihood of that output. 
  #    ...: Other arguments passed to ggplot(). 
  #
  #
  
  # If requested, scale input parameters based on `input_bounds`. Must save an unscaled version as well so that 
  # SSR/llik/lpost values can be computed by running the forward model on the original (unscaled) inputs. Also, 
  # the lprior values also require the unscaled inputs. 
  theta_vals_unscaled <- as.matrix(theta_vals)
  if(scale_inputs) {
    theta_vals <- scale_input_data(theta_vals, input_bounds)
    if(!is.null(point_coords)) point_coords <- scale_input_data(matrix(point_coords, nrow=1), input_bounds)
    if(!is.null(samples_kde)) samples_kde <- scale_input_data(samples_kde, input_bounds)
    if(!is.null(samples_points)) samples_points <- scale_input_data(samples_points, input_bounds)
  } 
  
  if(response_surface == "SSR") {
    
    if(is.null(SSR_vals)) SSR_vals <- get_computer_model_SSR(computer_model_data, theta_vals = theta_vals_unscaled)
    
    plts <- get_2d_response_surface_plot_SSR(theta_vals = theta_vals, 
                                             SSR_vals = SSR_vals, 
                                             param_names = param_names, 
                                             output_variables = output_variables, 
                                             raster = raster, 
                                             point_coords = point_coords, 
                                             samples_kde = samples_kde, 
                                             samples_points = samples_points)
  } else if(response_surface == "prior") {
    if(is.null(lprior_vals)) lprior_vals <- calc_lprior_theta(theta_vals_unscaled, theta_prior_params)
    
    plts <- get_2d_response_surface_plot_prior(theta_vals = theta_vals, 
                                               param_names = param_names, 
                                               raster = raster, 
                                               point_coords = point_coords, 
                                               samples_kde = samples_kde, 
                                               samples_points = samples_points,
                                               lprior_vals = lprior_vals, 
                                               theta_prior_params = theta_prior_params)
    
  } else if(response_surface == "likelihood") {
    if(is.null(SSR_vals)) SSR_vals <- get_computer_model_SSR(computer_model_data, theta_vals = theta_vals_unscaled)
    
    plts <- get_2d_response_surface_plot_likelihood(theta_vals = theta_vals,
                                                    computer_model_data = computer_model_data, 
                                                    param_names = param_names, 
                                                    output_variables = output_variables, 
                                                    combine_outputs = combine_outputs, 
                                                    raster = raster, 
                                                    point_coords = point_coords, 
                                                    samples_kde = samples_kde, 
                                                    samples_points = samples_points, 
                                                    SSR_vals = SSR_vals, 
                                                    llik_vals = llik_vals)
    
  }  else if(response_surface == "posterior") {
    
    if(is.null(SSR_vals) && is.null(llik_vals) && is.null(lpost_vals)) SSR_vals <- get_computer_model_SSR(computer_model_data, theta_vals = theta_vals_unscaled)
    if(is.null(lpost_vals) && is.null(lprior_vals)) lprior_vals <- calc_lprior_theta(theta_vals_unscaled, theta_prior_params)
    
    plts <- get_2d_response_surface_plot_posterior(theta_vals = theta_vals,
                                                   computer_model_data = computer_model_data, 
                                                   param_names = param_names, 
                                                   output_variables = output_variables, 
                                                   raster = raster, 
                                                   point_coords = point_coords,
                                                   samples_kde = samples_kde, 
                                                   samples_points = samples_points,
                                                   SSR_vals = SSR_vals, 
                                                   llik_vals = llik_vals, 
                                                   lprior_vals = lprior_vals, 
                                                   lpost_vals = lpost_vals, 
                                                   theta_prior_params = theta_prior_params)
  } else {
    stop("Invalid response surface: ", response_surface)
  }
  
  return(plts)
  
}


get_2d_response_surface_plot_SSR <- function(theta_vals, SSR_vals, param_names, output_variables = NULL, 
                                             raster = FALSE, point_coords = NULL, samples_kde = NULL,
                                             samples_points = NULL) {
  
  if(!is.null(output_variables)) {
    SSR_vals <- SSR_vals[, output_variables]
  }
  
  plts <- get_2d_heatmap_plots(X = theta_vals, 
                               Y = SSR_vals, 
                               param_names = param_names,
                               samples_kde = samples_kde, 
                               samples_points = samples_points, 
                               raster = raster, 
                               base_main_title = "SSR",
                               point_coords = point_coords, 
                               bigger_is_better = FALSE, 
                               legend_label = "SSR")
  
  return(plts)
  
}


get_2d_response_surface_plot_prior <- function(theta_vals, param_names, raster = FALSE, 
                                               point_coords = NULL, samples_kde = NULL, samples_points = NULL,
                                               lprior_vals = NULL, theta_prior_params = NULL) {
  # Even though this function is guaranteed to return a single plot, it returns a one element list to be 
  # consistent with the other response surface plotting functions. 
  
  # Log prior evaluations. 
  if(is.null(lprior_vals)) lprior_vals <- calc_lprior_theta(theta_vals, theta_prior_params)
  lprior_vals <- matrix(lprior_vals, ncol=1)
  
  # Heatmap plot. 
  plts <- get_2d_heatmap_plots(X = theta_vals, 
                               Y = lprior_vals, 
                               param_names = param_names,
                               samples_kde = samples_kde, 
                               samples_points = samples_points, 
                               raster = raster, 
                               main_title = "Log Prior", 
                               point_coords = point_coords,
                               bigger_is_better = TRUE, 
                               legend_label = "Log Prior")
  
  return(plts)
  
}


get_2d_response_surface_plot_likelihood <- function(theta_vals, computer_model_data, param_names, output_variables,  
                                                    combine_outputs = TRUE, raster = FALSE, point_coords = NULL, 
                                                    samples_kde = NULL, samples_points = NULL, SSR_vals = NULL, llik_vals = NULL) {
  
  # Log likelihood evaluations. 
  if(is.null(llik_vals)) {
    llik_vals <- llik_product_Gaussian(computer_model_data = computer_model_data, 
                                       vars_obs = diag(computer_model_data$Sig_eps), 
                                       theta_vals = theta_vals, 
                                       SSR = SSR_vals, 
                                       na.rm = TRUE, 
                                       sum_output_lliks = combine_outputs)
  }
  
  if(combine_outputs) {
    llik_vals <- matrix(llik_vals, ncol = 1)
    main_title <- "Log Likelihood"
  } else {
    main_title <- NULL
  }
  
  # Heatmap plot(s). 
  plts <- get_2d_heatmap_plots(X = theta_vals, 
                               Y = llik_vals, 
                               param_names = param_names,
                               samples_kde = samples_kde, 
                               samples_points = samples_points, 
                               raster = raster, 
                               base_main_title = "Log Likelihood", 
                               main_title = main_title, 
                               point_coords = point_coords,
                               bigger_is_better = TRUE, 
                               legend_label = "Log Likelihood")
  
  return(plts)
  
}


get_2d_response_surface_plot_posterior <- function(theta_vals, computer_model_data, param_names, output_variables, 
                                                   raster = FALSE, point_coords = NULL, samples_kde = NULL, 
                                                   samples_points = NULL, SSR_vals = NULL, llik_vals = NULL,  
                                                   lprior_vals = NULL, lpost_vals = NULL, theta_prior_params = NULL, 
                                                   sig2_eps = NULL, main_title = NULL, samples_points_size = 1) {
  
  # Log posterior evaluations. 
  if(is.null(lpost_vals)) {
    if(is.null(sig2_eps)) sig2_eps <- diag(computer_model_data$Sig_eps)
    
    lpost_vals <- calc_lpost_theta_product_lik(computer_model_data = computer_model_data, 
                                               theta_vals = theta_vals, 
                                               SSR = SSR_vals, 
                                               vars_obs = sig2_eps, 
                                               na.rm = TRUE, 
                                               theta_prior_params = theta_prior_params, 
                                               lprior_vals = lprior_vals,
                                               return_list = FALSE)
    
  }
  
  # Heatmap plot. 
  if(is.null(main_title)) main_title <- "Log Posterior"
  plt <- get_2d_heatmap_plot(X = theta_vals, 
                             y = lpost_vals, 
                             param_names = param_names,
                             samples_kde = samples_kde, 
                             samples_points = samples_points, 
                             raster = raster, 
                             main_title = main_title, 
                             point_coords = point_coords,
                             bigger_is_better = TRUE, 
                             legend_label = "Log Posterior", 
                             samples_points_size = samples_points_size)
  
  return(plt)
  
}


save_plots <- function(plts, type, save_dir) {
  
  for(i in seq_along(plts)) {
    plt_name <- paste0(type, "_", names(plts)[i], ".png")
    ggsave(filename=file.path(save_dir, plt_name), plot=plts[[i]])
  }
  
}


get_trace_plots_wide <- function(samp_df, burn_in_start = 1, ...) {
  # See `select_mcmc_samp_cols()` for details on how to specify columns which to plot. Operates on the 
  # wide MCMC data.frame format. 
  
  # Select columns to plot.  
  n <- nrow(samp_df)
  cols_sel <- select_mcmc_samp_cols(...)
  samp_df_plot <- samp_df %>% select(matches(cols_sel))
  col_names <- colnames(samp_df_plot)
  
  # Get starting iteration number for each parameter. 
  start_itrs <- get_mcmc_burn_in_start_itrs(burn_in_start, colnames(samp_df_plot))
  
  # Generate one trace plot per column. 
  plts <- vector(mode = "list", length = length(col_names))
  
  for(j in seq_along(plts)) {
    y <- col_names[j]
    itr_start <- start_itrs[j]
    df <- samp_df_plot[itr_start:n, y, drop = FALSE]
    colnames(df) <- "param"
    df$itr <- itr_start:n
    
    plts[[j]] <- ggplot(data = df, aes(x = itr, y = param)) + 
      geom_line() + 
      ggtitle(paste0("Trace Plot: ", y)) + 
      xlab("Iteration")
  }
  
  return(plts)
  
}


get_mcmc_marginal_hist_plot <- function(samp_df, param_names, burn_in_start=1, bins=30, vertical_lines=NULL, ...) {
  # Generates one plot per parameter name. 
  
  samp_df_plot <- select_mcmc_samp(mcmc_samp_df=samp_df, param_names=param_names, burn_in_start=burn_in_start, ...)
  
  # Produce one plot per parameter name. 
  plts <- vector(mode = "list", length = length(param_names))
  for(j in seq_along(plts)) {
    param_name <- param_names[j]
    
    plts[[j]] <- ggplot(data = samp_df_plot, aes(x = value, color = variable)) + 
      geom_histogram(aes(y = ..density..), bins = bins, fill = "white", alpha = 0.2, position = "identity") + 
      xlab(param_names[j]) + 
      ylab("Frequency") + 
      ggtitle(paste0("Marginal Distribution: ", param_names[j]))
    
  }
  
  
  # Select columns to plot.  
  n <- nrow(samp_df)
  cols_sel <- select_mcmc_samp_cols(param_names = param_names, ...)
  samp_df_plot <- samp_df %>% select(matches(cols_sel))
  col_names <- colnames(samp_df_plot)
  
  # Get starting iteration number for each parameter. 
  start_itrs <- get_mcmc_burn_in_start_itrs(burn_in_start, colnames(samp_df_plot))
  
  # Convert data.frame to long format. 
  samp_df_plot <- melt(as.data.table(samp_df_plot), measure.vars = colnames(samp_df_plot), na.rm = TRUE)
  
  # Produce one plot per parameter name. 
  plts <- vector(mode = "list", length = length(param_names))
  for(j in seq_along(plts)) {
    df_param_cols <- grep(param_names[j], col_names, value = TRUE)
    plts[[j]] <- ggplot(data = samp_df_plot[variable %in% df_param_cols], aes(x = value, color = variable)) + 
      geom_histogram(aes(y = ..density..), bins = bins, fill = "white", alpha = 0.2, position = "identity") + 
      xlab(param_names[j]) + 
      ylab("Frequency") + 
      ggtitle(paste0("Marginal Distribution: ", param_names[j]))
    
    if(!is.null(vertical_lines)) {
      plts[[j]] <- plts[[j]] + geom_vline(xintercept = vertical_lines[param_names[j]], color = "red")
    }
    
  }
  
  return(plts)
  
}














