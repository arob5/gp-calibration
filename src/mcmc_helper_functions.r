#
# mcmc_helper_functions.r
# Primarily functions to post-process, organize, and plot MCMC samples. 
#
# Andrew Roberts
#
# Dependencies:
#    general_helper_functions.r, plotting_helper_functions.r

library(HDInterval)
library(kde1d)

# ------------------------------------------------------------------------------
# Description of required data.table formatting:
#
# Bespoke functions for storing and manipulating MCMC output geared towards the  
# use case of comparing different MCMC algorithms. The primary object that 
# these functions operate on is a data.table that by convention is named 
# `samp_dt`. This data.table, storing the MCMC samples, is required to have the 
# following columns: test_label, chain_idx, param_type, param_name, itr, sample.
#    test_label: a character uniquely identifying an MCMC "test". This is 
#                typically used to differentiate different algorithms, but it 
#                could also allow for other groupings (e.g., the same algorithm
#                run with different settings).
#    chain_idx: Integer used to index different MCMC chains/trajectories within 
#               a single MCMC test. Must be unique within a single `test_label`.
#    param_type: character, used to define groups of parameters. For example, 
#                output from an MCMC run might consist of samples from a 
#                coefficient vector and a variance parameter. Two `param_type`
#                labels could be defined to group these two types of parameters.
#                The `param_type` labels must be unique within a unique 
#                (test_label, chain) combination.
#    param_name: character, the parameter names. Must be unique within a unique
#                (test_label, chain, param_type) combination.
#    itr: integer, the MCMC iteration index. Index starts from 1. Must be unique
#         within a unique (test_label, chain, param_type, param_name) combination.
#    sample: numeric, the sample values.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Functions for creating `samp_dt` objects.
# ------------------------------------------------------------------------------

assert_is_samp_dt <- function(samp_dt) {
  # Checks that the argument `samp_dt` satisfies the requirements for the 
  # data.table storing MCMC output.
  # TODO: add checks for the uniqueness requirements.
  #
  # Args:
  #    samp_dt: object to check.
  #
  # Returns:
  #  None, will throw exception if `samp_dt` fails any of the tests.
  
  required_cols <- c("test_label"="character", "chain_idx"="integer", 
                     "param_type"="character", "param_name"="character", 
                     "itr"="integer", "sample"="numeric")
  
  assert_that(is.data.table(samp_dt))
  assert_that(setequal(colnames(samp_dt), names(required_cols)))
  assert_that(all(sapply(samp_dt, class)[names(required_cols)] == required_cols))
}

format_mcmc_output <- function(samp_list, test_label, chain_idx=1L) {
  # Formats output from a single MCMC run (i.e., a single test label and chain)
  # into the required data.table format. 
  #
  # Args:
  #    samp_list: named list, with one element per parameter type. Each element 
  #               is a matrix with rows containing MCMC samples of the respective
  #               parameter type, and columns corresponding to different 
  #               parameters falling within that type. The names of `samp_list`
  #               must be set to the parameter type names. Given that the samples 
  #               are all from the same MCMC run, the matrices will all 
  #               typically have the same number of rows. However, this is not 
  #               required.
  #    test_label: character, the test label uniquely identifying the MCMC run.
  #    chain_idx: integer, the index used to uniquely identify different MCMC chains
  #               within a single MCMC run. This function assumes only a single 
  #               chain has been run; for mutliple chains, see 
  #               `format_mcmc_output_multi_chain()`.
  #
  # Returns:
  # data.table with column names:
  # "test_label", "chain_idx", param_type", "param_name", "itr", and "sample". 
  # The columns "test_label" and "chain_idx" will be constant with values set to 
  # the `test_label` and `chain` function arguments, respectively. The values 
  # in column "param_type" are taken from the names of `samp_list`, while the 
  # values in column "param_name" are taken from the column names in the 
  # matrices within `samp_list`. The column `itr` contains the integer MCMC 
  # iteration. The column `sample` contains the MCMC numeric sample values. 
  
  samp_dt <- data.table(param_type=character(), itr=integer(), 
                        param_name=character(), sample=numeric())
  
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
  
  # Add test label and chain index.
  samp_dt[, test_label := test_label]
  samp_dt[, chain_idx := as.integer(chain_idx)]
  
  return(samp_dt)
}


append_mcmc_output <- function(samp_dt, samp_list, test_label, chain_idx=1L) {
  # Appends a new MCMC test to an existing `samp_dt` object. The new test
  # label `test_label` must not already be present in `samp_dt`. This is 
  # a convenience wrapper that essentially calls 
  # `format_mcmc_output(samp_list, test_label)` and then appends the result
  # to `samp_dt`. This function only appends a single chain; for multiple
  # chains see `append_mcmc_output_multi_chain()`.
  # 
  # Args:
  #    samp_dt: existing data.table object storing samples.
  #    samp_list, test_label, chain_idx: all passed to `format_mcmc_output()`.
  #
  # Returns:
  #    The data.table `samp_dt` with the new MCMC test appended.
  
  assert_is_samp_dt(samp_dt)

  # Ensure test label is not already in list.
  assert_that(!(test_label %in% samp_dt$test_label))
  
  samp_dt_new <- format_mcmc_output(samp_list, test_label=test_label, 
                                    chain_idx=chain_idx)
  samp_dt <- rbindlist(list(samp_dt, samp_dt_new), use.names=TRUE)
  
  return(samp_dt)
}

append_chain <- function(samp_dt, samp_list, test_label, chain_idx=NULL) {
  # Appends sample output to `samp_dt` corresponding to a new chain of an MCMC 
  # test that is already present in `samp_dt`.
  #
  # Args:
  #    samp_dt: existing data.table object storing samples.
  #    samp_list: list satisfying the requirements specified in 
  #               `format_mcmc_output()`. The sample output from a new chain/
  #               trajectory of the MCMC test identified by `test_label`.
  #    test_label: character, a test label already present in `samp_dt`.
  #    chain_idx: integer, the index for the new chain. If NULL, looks for the 
  #               current maximum chain index for the specified test label, and 
  #               increments by 1.
  #
  # Returns:
  #  The data.table `samp_dt` with the output from the new chain appended.
  
  assert_is_samp_dt(samp_dt)
  
  # Ensure the test is already present.
  assert_that(test_label %in% samp_dt$test_label)

  # If chain index not specified, increment current maximum index by 1.
  if(is.null(chain_idx)) {
    chain_idx <- 1L + samp_dt[test_label==test_label, max(chain_idx)]
  } else {
    assert_that(is.integer(chain_idx))
    assert_that(!(chain_idx %in% samp_dt[test_label==test_label, unique(chain_idx)]))
  }
  
  # Format the MCMC output for the new chain.
  samp_dt_new <- format_mcmc_output(samp_list, test_label, chain_idx=chain_idx)
  
  # Append new output to current data.table.
  samp_dt <- rbindlist(list(samp_dt, samp_dt_new), use.names=TRUE)
  
  return(samp_dt)
}


format_mcmc_output_multi_chain <- function(chain_samp_list, test_label) {
  # A wrapper around `format_mcmc_output()` that allows for appending multiple 
  # chains. See `format_mcmc_output()` for more details.
  #
  # Args:
  #    chain_samp_list: a list, where each element is a `samp_list`, as 
  #                     described in `format_mcmc_output()`; i.e., each element
  #                     of the outer list corresponds to a different chain.
  #    test_label: character, the unique label identifying the MCMC test.
  #
  # Returns:
  #  A `samp_dt` data.table. The chain indices are set to 1,2,...,n_chains.
  
  assert_that(is.list(chain_samp_list))
  n_chains <- length(chain_samp_list)
  
  # First create the data.table for only the first chain.
  samp_dt <- format_mcmc_output(chain_samp_list[[1]], test_label, chain_idx=1L)
  if(n_chains==1L) return(samp_dt)
  
  # Now loop over remaining chains and append to data.table.
  for(idx in 2:n_chains) {
    samp_dt <- append_chain(samp_dt, chain_samp_list[[idx]], test_label, chain_idx=idx)
  }
  
  return(samp_dt)
}


append_mcmc_output_multi_chain <- function(samp_dt, chain_samp_list, test_label) {
  # A wrapper around `append_mcmc_output()` that allows for appending a new MCMC
  # test containing multiple chains. See `format_mcappend_mcmc_output()` for more 
  # details.
  # 
  # Args:
  #    samp_dt: existing data.table object storing samples.
  #    chain_samp_list: a list, where each element is a `samp_list`, as 
  #                     described in `format_mcmc_output()`; i.e., each element
  #                     of the outer list corresponds to a different chain.
  #    test_label: character, the label uniquely identifying the new MCMC test.
  #                Must not already be present in `samp_dt`.
  #
  # Returns:
  #    The data.table `samp_dt` with the new MCMC test appended.
  
  assert_is_samp_dt(samp_dt)
  
  # Ensure test label is not already in list.
  assert_that(!(test_label %in% samp_dt$test_label))
  
  samp_dt_new <- format_mcmc_output_multi_chain(chain_samp_list, test_label)
  samp_dt <- rbindlist(list(samp_dt, samp_dt_new), use.names=TRUE)
  
  return(samp_dt)
}


# ------------------------------------------------------------------------------
# Functions for manipulating/subsetting `samp_dt` objects.
# ------------------------------------------------------------------------------

select_mcmc_samp <- function(samp_dt, test_labels=NULL, param_types=NULL,
                             param_names=NULL, itr_start=1L, itr_stop=NULL, 
                             chain_idcs=NULL) {
  # Selects rows from `samp_dt` corresponding to valid combinations of 
  # `test_labels`, `param_types`, and `param_names`. Also restricts the output
  # to the iteration range specified by `itr_start` and `itr_stop`.
  #
  # Args:
  #    samp_dt: data.table of MCMC samples. 
  #    itr_start: the minimum iteration number to be included. Named vector 
  #               allows varying the start iteration by test label. See 
  #               `select_mcmc_itr()` for details.
  #    itr_stop: Same as `itr_start` but specifying the upper iteration cutoff.
  #    The remaining arguments are vectors of values used to subset `samp_dt`
  #    by selecting values from the columns "test_label", "param_type", 
  #    "param_name", and "chain_idx". NULL value will include all values found 
  #    in the respective column in `samp_dt`. 
  #
  # Returns:
  #    data.table, containing subset of rows from `samp_dt`. 
  
  samp_dt_subset <- copy(samp_dt)
  
  # If not provided, select all. 
  if(is.null(test_labels)) test_labels <- samp_dt_subset[, unique(test_label)]
  if(is.null(param_types)) param_types <- samp_dt_subset[, unique(param_type)]
  if(is.null(param_names)) param_names <- samp_dt_subset[, unique(param_name)]
  if(is.null(chain_idcs)) chain_idcs <- samp_dt_subset[, unique(chain_idx)]
  
  # Select rows corresponding to label-type-name combinations. 
  samp_dt_subset <- samp_dt[(test_label %in% test_labels) & 
                            (param_type %in% param_types) & 
                            (param_name %in% param_names) &
                            (chain_idx %in% chain_idcs)]
  
  # Restrict to specified iteration bounds.
  samp_dt_subset <- select_mcmc_itr(samp_dt_subset, itr_start, itr_stop)

  return(samp_dt_subset)
}


select_mcmc_samp_mat <- function(samp_dt, test_label, param_type,
                                 param_names=NULL, itr_start=1L, itr_stop=NULL,
                                 chain_idcs=NULL, return_chain_list=FALSE) {
  # A wrapper around `select_mcmc_samp()` that converts the selected samples 
  # to a matrix format, where each row is a sample. Unlike `select_mcmc_samp()`,
  # a single value of `test_label` and `param_type` must be specified. However,
  # `param_names` can be used to return a subset of the parameters for the 
  # specified parameter type. It is recommended to pass this argument 
  # even if all parameters of a specific type are to be returned, as passing 
  # `param_names` will order to the columns according to `param_names`, which 
  # ensures that the ordering is correct. Also allows subsetting by the chain 
  # indices. Can either return a list of matrices (one per chain), or a single 
  # matrix in which all chains have been grouped together.
  #
  # Args:
  #    samp_dt: data.table of MCMC samples. 
  #    test_label: character, the single test label to select. 
  #    param_type: character, the single parameter type to select. Must be a  
  #                type within the specified `test_label`. 
  #    itr_start: the minimum iteration number to be included. Named vector 
  #               allows varying the start iteration by test label. See 
  #               `select_mcmc_itr()` for details.
  #    itr_stop: Same as `itr_start` but specifying the upper iteration cutoff.
  #    param_names: character, if provided then selects the specific parameter 
  #                 names provided, and also orders the columns of the returned 
  #                 matrix in the order they are provided in this argument. The 
  #                 parameters specified here must be of the type specified 
  #                 by `param_type` and contained within the test label specified 
  #                 by `test_label`. If NULL, selects all parameters within the 
  #                 param type, and no explicit ordering is performed. 
  #    chain_idcs: integer, the chain indices to select.
  #    return_chain_list: logical, if TRUE returns a list with each element 
  #                       containing a matrix storing the samples for a specific
  #                       chain. Otherwise, all samples are put in a single 
  #                       matrix. See "Returns" for details.
  #
  # Returns:
  #    If `return_chain_list` is TRUE, a list with each element containing a 
  #    matrix storing the samples for a specific chain. The list names are set 
  #    to the chain indices. The rownames of each matrix is set to the iteration
  #    numbers. If return_chain_list` is FALSE, all samples are put in a single 
  #    matrix. In this case, the rownames of the matrix take the form 
  #    "<chain_idx>_<itr>". In either case, each row of the matrix corresponds 
  #    to a single sample, and each column to a parameter within the selected 
  #    parameter type. The columns are ordered according to `param_names`, if 
  #    this argument is provided. 
  
  assert_that(length(test_label)==1L)
  assert_that(length(param_type)==1L)
  assert_is_samp_dt(samp_dt)
  
  # Select specified subset of `samp_dt`. 
  samp_dt_subset <- select_mcmc_samp(samp_dt, test_labels=test_label, 
                                     param_types=param_type, 
                                     param_names=param_names, itr_start=itr_start,
                                     itr_stop=itr_stop, chain_idcs=chain_idcs)
  if(nrow(samp_dt_subset)==0L) {
    stop("Cant convert zero length data.table to wide matrix format.")
  }
  
  # Create one matrix per chain.
  chain_idcs <- samp_dt_subset[,unique(chain_idx)]
  mat_chain_list <- lapply(chain_idcs, 
                           function(i) convert_samp_to_mat(samp_dt_subset[chain_idx==i]))
  
  # If `param_names` is provided, sort columns based on this order.
  if(!is.null(param_names)) {
    for(i in seq_along(mat_chain_list)) {
      mat_chain_list[[i]] <- mat_chain_list[[i]][,param_names, drop=FALSE]
    }
  }
    
  # If returning a single matrix, stack all matrices into one.
  if(!return_chain_list) {
    for(i in seq_along(chain_idcs)) {
      rownames(mat_chain_list[[i]]) <- paste(chain_idcs[i], 
                                           rownames(mat_chain_list[[i]]),
                                           sep="_")
    }
    return(do.call(rbind, mat_chain_list))
  }
  
  names(mat_chain_list) <- chain_idcs
  return(mat_chain_list)
}


select_mcmc_itr <- function(samp_dt, itr_start=1L, itr_stop=NULL) {
  # Returns a subset of the rows of `samp_dt` obtained by restricting the 
  # iteration column `itr` to rows with values falling in the interval 
  # [itr_start, itr_stop]. Alternatively, such intervals can be specified 
  # by test label, so that different iteration ranges can be selected for 
  # different test labels. See `Args` for details.
  #
  # Args:
  #    samp_dt: data.table of MCMC samples. 
  #    itr_start: integer vector specifying the lower iteration cutoff. If integer 
  #               of length 1 AND unnamed, this is interpreted as the starting 
  #               iteration for all test labels - all earlier iterations
  #               are dropped. If a named vector, must have names set to 
  #               valid test label values. This allows application of a different 
  #               burn-in start iteration for different test labels. If only 
  #               some test labels are specified by the vector, then the others 
  #               will not be affected. 
  #    itr_stop: same as `itr_start` but defining the upper cutoff. A NULL value
  #              implies that no upper cutoff should be enforced.
  #
  # Returns:
  #    data.table, a copy of `samp_dt` which is row subsetted to select the 
  #    specified iteration ranges.
  
  assert_is_samp_dt(samp_dt)
  
  # Argument checking.
  if(is.null(itr_stop)) itr_stop <- samp_dt[, max(itr)]
  assert_that(all(itr_start >= 1L))
  assert_that(is.integer(itr_start) && is.integer(itr_stop))
  
  # If no row subsetting is required just return the data.table.
  if(all(itr_start==1L) && is.null(itr_stop)) return(samp_dt)
  
  # If `itr_stop` is NULL, set to maximum iteration.
  if(is.null(itr_stop)) itr_stop <- samp_dt[,max(itr)]
  
  # Get set of unique test labels.
  test_labels <- samp_dt[,unique(test_label)]
  n_labels <- length(test_labels)
  
  # Helper function for processing the arguments `itr_start` and `itr_stop`.
  get_itr_bound <- function(itr_arg) {
    # If `itr_start` and/or `itr_stop` are unnamed vectors of length 1, then 
    # the iteration cutoff will be applied to all test labels.
    if((length(itr_arg)==1L) && is.null(names(itr_arg))) {
      return(setNames(rep(itr_arg,n_labels), test_labels))
    }
    
    # If `itr_start` and/or `itr_stop` are named vectors, apply values separately 
    # for each test label.
    extra_labels <- setdiff(names(itr_arg), test_labels)
    missing_labels <- setdiff(test_labels, names(itr_arg))
    if(length(extra_labels) > 0) message("Extra labels not used: ", 
                                         paste(extra_labels, collapse=", "))
    if(length(missing_labels) > 0) message("Labels not affected: ", 
                                           paste(missing_labels, collapse=", "))
    return(itr_arg[intersect(test_labels,names(itr_arg))])
  }
  
  itr_start <- get_itr_bound(itr_start)
  itr_stop <- get_itr_bound(itr_stop)

  # Subset `samp_dt` by selecting iteration ranges by test label.
  samp_dt_subset <- copy(samp_dt)
  for(lbl in names(itr_start)) {
    samp_dt_subset <- samp_dt_subset[(test_label != lbl) | (itr >= itr_start[lbl])]
  }
  for(lbl in names(itr_stop)) {
    samp_dt_subset <- samp_dt_subset[(test_label != lbl) | (itr <= itr_stop[lbl])]
  }

  return(samp_dt_subset)
}


convert_samp_to_mat <- function(samp_dt) {
  # A helper function used by `select_mcmc_samp_mat()` to convert sample output
  # in data.table form for a single test label/parameter type combination to 
  # a matrix format.
  #
  # Args:
  #    samp_dt: data.table in the typical form. Only columns "itr", "param_name", 
  #             and "sample" are used here.
  #
  # Returns:
  #    matrix, with samples in the rows and columns corresponding to the 
  #    different parameters. The rownames attribute is set to the iterations.
  samp_wide <- dcast(data=samp_dt[,.(itr, param_name, sample)], 
                     formula=itr~param_name, value.var="sample")
  itrs <- samp_wide$itr
  samp_wide <- as.matrix(samp_wide[, .SD, .SDcols=!"itr"])
  rownames(samp_wide) <- itrs
  return(samp_wide)
}


# ------------------------------------------------------------------------------
# Computing statistics and errors from MCMC samples. 
# ------------------------------------------------------------------------------

compute_mcmc_param_stats <- function(samp_dt, burn_in_start=NULL, test_labels=NULL, param_types=NULL,
                                     param_names=NULL, subset_samp=TRUE, format_long=FALSE) {
  # Currently just computes sample means and variances for the selected parameters/variables in `samp_dt`. 
  #
  # Args:
  #   samp_dt: data.table, must be of the format described in `format_mcmc_output()`.
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

get_trace_plots <- function(samp_dt, test_labels=NULL, param_types=NULL,  
                            param_names=NULL, itr_start=1L, itr_stop=NULL, 
                            chain_idcs=NULL, save_dir=NULL, overlay_chains=TRUE) {
  # Extracts the specified subset of `samp_dt` then generates one trace plot per
  # valid `param_name`-`param_type`-`test_label`-`chain_idx` combination. If 
  # `overlay_chains` is TRUE, then the chains within each unique
  # `param_name`-`param_type`-`test_label` combination are all plotted on the 
  # same plot with different colors.
  #
  # Args:
  #    samp_dt: existing data.table object storing samples.
  #    overlay_chains: If TRUE, chains are overlaid on the same plot, otherwise
  #                    each chain is plotted on a different plot.
  #    save_dir: character(1), if not NULL, a file path to save the plots to. 
  #    Remaining arguments are used to subset `samp_dt`; see `select_mcmc_samp()`.
  #
  # Returns:
  #    list of ggplot objects, the plots as described above. 
  
  assert_is_samp_dt(samp_dt)
  
  # Determine which plots to create by subsetting rows of `samp_dt`. 
  samp_dt_subset <- select_mcmc_samp(samp_dt, test_labels=test_labels,
                                     param_types=param_types, 
                                     param_names=param_names, itr_start=itr_start,
                                     itr_stop=itr_stop, chain_idcs=chain_idcs)
  
  # Determine which combination of columns will uniquely identify a plot.
  id_cols <- c("test_label", "param_type", "param_name")
  if(!overlay_chains) id_cols <- id_cols <- c(id_cols, "chain_idx")
  plt_id_vars <- unique(samp_dt_subset[, ..id_cols])
  
  # Generate plots. 
  plts <- list()
  for(j in 1:nrow(plt_id_vars)) {
    # Title and label for plot.
    id_vals <- plt_id_vars[j]
    plt_label <- paste(as.matrix(id_vals)[1,], collapse="_")
    plt_title <- paste(id_vals$test_label, 
                       paste0(id_vals$param_type, ": ", id_vals$param_name),
                       sep=" | ")
    if(!overlay_chains) plt_title <- paste(plt_title, 
                                           paste0("Chain ", id_vals$chain_idx), 
                                           sep=" | ")
    
    # Select data for plot.
    chain_idx <- NULL
    if(!overlay_chains) chain_idx <- id_vals$chain_idx
    samp_dt_plt <- select_mcmc_samp(samp_dt_subset, test_labels=id_vals$test_label,
                                    param_types=id_vals$param_type, 
                                    param_names=id_vals$param_name, 
                                    chain_idcs=chain_idx)
    
    # Create trace plot.
    if(overlay_chains) {
      plt <- ggplot(samp_dt_plt, aes(x=itr, y=sample, color=as.factor(chain_idx)))
    } else {
      plt <- ggplot(samp_dt_plt, aes(x=itr, y=sample))
    }
    
    plt <- plt + geom_line() + ggtitle(plt_title) + xlab("Iteration")
    if(overlay_chains) plt <- plt + guides(color=guide_legend(title="Chain"))
    plts[[plt_label]] <- plt
  }
  
  # Optionally save to file.
  if(!is.null(save_dir)) save_plots(plts, "trace", save_dir)
  
  return(invisible(plts))
}


get_hist_plots <- function(samp_dt, test_labels=NULL, param_types=NULL,  
                           param_names=NULL, itr_start=1L, itr_stop=NULL, 
                           chain_idcs=NULL, save_dir=NULL, combine_chains=TRUE,
                           plot_type="hist", bins=30) {
  # Plots one histogram plot per unique `param_type`-`param_name` combination.
  # All MCMC runs containing output for this parameter are overlaid on the same 
  # plot (see Args and Returns for specifics on this). These plots can thus get
  # quite cluttered with more than a few histograms on a single plot. Using
  # `plot_type = "freqpoly"` can be helpful to reduce clutter in this case. If 
  # there is a clear baseline against which all other MCMC runs should be 
  # compared, see `get_hist_plot_comparisons()` instead.
  #
  # Args:
  #    samp_dt: existing data.table object storing samples.
  #    combine_chains: logical. If TRUE, samples from all chain indices are 
  #                    pooled together. Otherwise, distinct chains are plotted 
  #                    separately, overlaid on the same plot.
  #    plot_type: character, either "hist" or "freqpoly" to produce a histogram
  #               or frequency polygon.
  #    bins: integer, passed to either geom_histogram() or geom_freqpoly()
  #          `bins` argument.
  #    save_dir: character(1), if not NULL, a file path to save the plots to. 
  #    Remaining arguments are used to subset `samp_dt`; see `select_mcmc_samp()`.
  #
  # Returns:
  # list of ggplot objects, one plot per unique `param_type`-`param_name` 
  # combination.
  
  assert_is_samp_dt(samp_dt)
  assert_that(plot_type %in% c("hist","freqpoly"))
  
  # Determine which plots to create by subsetting rows of `samp_dt`. 
  samp_dt_subset <- select_mcmc_samp(samp_dt, test_labels=test_labels,
                                     param_types=param_types, 
                                     param_names=param_names, itr_start=itr_start,
                                     itr_stop=itr_stop, chain_idcs=chain_idcs)
  
  # Determine which combination of columns will uniquely identify a plot.
  id_cols <- c("param_type", "param_name")
  plt_id_vars <- unique(samp_dt_subset[, ..id_cols])
  legend_label <- ifelse(combine_chains, "test lbl", "test lbl + chain")
  
  # Plot histogram of frequency polygon.
  if(plot_type == "hist") {
    plt_func <- function() geom_histogram(fill="white", alpha=0.2, 
                                          bins=bins, position="identity")
  } else if(plot_type == "freqpoly") {
    plt_func <- function() geom_freqpoly(bins=bins, position="identity")
  }
  
  # Generate plots. 
  plts <- list()
  for(j in 1:nrow(plt_id_vars)) {
    # Title and label for plot.
    id_vals <- plt_id_vars[j]
    plt_label <- paste(as.matrix(id_vals)[1,], collapse="_")
    plt_title <- paste0(id_vals$param_type, ": ", id_vals$param_name)
    param_name <- id_vals$param_name

    # Select data for plot.
    samp_dt_plt <- select_mcmc_samp(samp_dt_subset,
                                    param_types=id_vals$param_type, 
                                    param_names=id_vals$param_name)

    # Create histogram plot.
    if(combine_chains) {
      plt <- ggplot(samp_dt_plt, aes(x=sample, color=test_label))
    } else {
      plt <- ggplot(samp_dt_plt, 
                    aes(x=sample, color=interaction(test_label, 
                                                    as.factor(chain_idx), sep=":")))
    }
    
    plt <- plt + plt_func() + 
                 ggtitle(plt_title) + 
                 guides(color=guide_legend(title=legend_label)) +
                 xlab(param_name)
    plts[[plt_label]] <- plt
  }
  
  # Optionally save to file.
  if(!is.null(save_dir)) save_plots(plts, "hist", save_dir)
  
  return(invisible(plts))
}


get_hist_plot_comparisons <- function(samp_dt, test_label_baseline=NULL, 
                                      test_labels=NULL, param_types=NULL,  
                                      param_names=NULL, itr_start=1L, 
                                      itr_stop=NULL, chain_idcs=NULL, 
                                      save_dir=NULL, combine_chains=TRUE,
                                      plot_type="hist", bins=30) {
  # A convenience wrapper around `get_hist_plots()` that generates figures with 
  # two histograms on each plot: one coming from a "baseline" test label, and 
  # the other from any other test label. This is helpful if comparing algorithms
  # against some ground truth or baseline samples. If the baseline test label 
  # has certain parameters that are not in the other test label, then the 
  # plot will only contain the histogram for the baseline.
  #
  # Args:
  #    samp_dt: existing data.table object storing samples.
  #    test_label_baseline: character(1) or NULL. If non-NULL, then must be a 
  #                         valid test label with associated samples in 
  #                         `samp_dt`. In this case, the histograms 
  #                         corresponding to this baseline test label will be 
  #                         overlaid on the plots for all other test labels.
  #    combine_chains: logical. If TRUE, samples from all chain indices are 
  #                    pooled together. Otherwise, distinct chains are plotted 
  #                    separately, overlaid on the same plot.
  #    plot_type: character, either "hist" or "freqpoly" to produce a histogram
  #               or frequency polygon.
  #    bins: integer, passed to either geom_histogram() or geom_freqpoly()
  #          `bins` argument.
  #    save_dir: character(1), if not NULL, a file path to save the plots to. 
  #    Remaining arguments are used to subset `samp_dt`; see `select_mcmc_samp()`. 
  #
  # Returns: 
  #    list, each element being a ggplot object. 
  
  # Determine which plots to create by subsetting rows of `samp_dt`. 
  if(!is.null(test_label_baseline) && !is.null(test_labels) && 
     !(test_label_baseline %in% test_labels)) {
    test_labels <- c(test_labels, test_label_baseline)
  }
  
  samp_dt_subset <- select_mcmc_samp(samp_dt, test_labels=test_labels, 
                                     param_types=param_types, 
                                     param_names=param_names, 
                                     itr_start=itr_start, itr_stop=itr_stop,
                                     chain_idcs=chain_idcs)
  
  # Determine which combination of columns will uniquely identify a plot.
  id_cols <- c("test_label", "param_type", "param_name")
  plt_id_vars <- unique(samp_dt_subset[, ..id_cols])

  # Separate out baseline label. 
  if(!is.null(test_label_baseline)) {
    plt_id_vars <- plt_id_vars[test_label != test_label_baseline]
  }
  
  # Generate plots. 
  plts <- list()
  for(j in 1:nrow(plt_id_vars)) {
    # Current test label.
    test_lbl_curr <- c(plt_id_vars[j, test_label], test_label_baseline)

    # Produce histograms for current test label.
    plts_curr <- get_hist_plots(samp_dt_subset, test_labels=test_lbl_curr,
                                combine_chains=combine_chains, 
                                plot_type=plot_type, bins=bins)
    plts <- c(plts, plts_curr)
  }
  
  if(!is.null(save_dir)) save_plots(plts, "hist_comparison", save_dir)
  return(invisible(plts))
}


get_1d_kde_plots <- function(samp_dt, test_label_baseline=NULL, 
                             test_labels=NULL, param_types=NULL,  
                             param_names=NULL, itr_start=1L, 
                             itr_stop=NULL, chain_idcs=NULL, 
                             save_dir=NULL, combine_chains=TRUE,
                             N_kde_pts=100, min_q=0.001, max_q =.999,
                             bandwidth_mult=1.0) {
  # Returns plots with kernel density estimates (KDE) for 1-dimensional marginal
  # distributions. Produces one plot per unique param type-param name combination.
  # Each plot contains one line per `test_label` that has samples for that 
  # specific parameter. Each line is produced from a 1d KDE constructed from the
  # samples using the kde1d package. If `test_label_baseline` is provided, then 
  # this test label will be treated as the baseline for comparison and plotted 
  # as a black dashed line.
  #
  # Args:
  #    samp_dt: existing data.table object storing samples.
  #    test_label_baseline: character(1) or NULL. If non-NULL, then must be a 
  #                         valid test label with associated samples in 
  #                         `samp_dt`. Designates a special test label, which 
  #                         typically indicates a baseline/ground truth.
  #    combine_chains: logical. Currently only supports TRUE, which groups all
  #                    chains together.
  #    N_kde_pts: integer, number of points at which the KDE approximation is 
  #               evaluated for each univariate parameter.
  #    min_q, max_q: numeric in (0,1], percentiles used as cutoffs to define the 
  #                  interval over which the KDE is computed.
  #    bandwidth_mult: passed to the `mult` arg of `kde1d()`.
  #    save_dir: character(1), if not NULL, a file path to save the plots to. 
  #    Remaining arguments are used to subset `samp_dt`; see `select_mcmc_samp()`.
  #
  # Returns:
  # list of ggplot objects, one per param type-param name combination.
  
  assert_is_samp_dt(samp_dt)
  
  if(!combine_chains) {
    stop("Currently `get_1d_kde_plots()` only supports `combine_chains = TRUE`.")
  }
  
  # Determine which plots to create by subsetting rows of `samp_dt`.
  test_labels_all <- test_labels
  if(!is.null(test_label_baseline) && !is.null(test_labels) && 
     !(test_label_baseline %in% test_labels)) {
    test_labels_all <- c(test_labels, test_label_baseline)
  }
  
  samp_dt_subset <- select_mcmc_samp(samp_dt, test_labels=test_labels_all, 
                                     param_types=param_types, 
                                     param_names=param_names, 
                                     itr_start=itr_start, itr_stop=itr_stop,
                                     chain_idcs=chain_idcs)
  
  # Separate out baseline label. 
  if(!is.null(test_label_baseline)) {
    samp_dt_baseline <- samp_dt_subset[test_label == test_label_baseline]
    samp_dt_subset <- samp_dt_subset[test_label != test_label_baseline]
  }
  
  # One plot is generated for each unique `param_type`-`param_name` combination.
  id_cols <- c("param_type", "param_name")
  plt_id_vars <- unique(samp_dt_subset[, ..id_cols])

  # Generate plots. 
  plts <- list()
  for(j in 1:nrow(plt_id_vars)) {
    # Title and label for plot.
    id_vals <- plt_id_vars[j]
    plt_label <- paste(as.matrix(id_vals)[1,], collapse="_")
    plt_title <- paste0(id_vals$param_type, ": ", id_vals$param_name)
    param_name <- id_vals$param_name
    
    # Prepare non-baseline data.
    samp_dt_param <- samp_dt_subset[(param_type == id_vals$param_type) & 
                                    (param_name == id_vals$param_name), 
                                    .(sample, test_label)]
    test_labels_curr <- samp_dt_param[,unique(test_label)]
    
    # Prepare baseline data. 
    samp_dt_baseline_param <- NULL
    if(!is.null(test_label_baseline)) {
      samp_dt_baseline_param <- samp_dt_baseline[(param_type == id_vals$param_type) & 
                                                 (param_name == id_vals$param_name), sample] 
    }
    
    # Determine the grid of points at which the KDE will be evaluated. 
    bound_lower <- quantile(c(samp_dt_param$sample, samp_dt_baseline_param), min_q)
    bound_upper <- quantile(c(samp_dt_param$sample, samp_dt_baseline_param), max_q)
    kde_pts <- seq(bound_lower, bound_upper, length.out=N_kde_pts)
    
    # Loop over test labels, constructing KDE for each label. 
    kde_mat <- matrix(nrow=N_kde_pts, ncol=length(test_labels_curr), 
                      dimnames=list(NULL, test_labels_curr))
    for(lbl in test_labels_curr) {
      kde_fit <- kde1d(samp_dt_param[test_label==lbl, sample], mult=bandwidth_mult)
      kde_mat[,lbl] <- dkde1d(kde_pts, kde_fit)
    }
    
    # Add KDE for baseline label. 
    kde_baseline <- NULL
    if(!is.null(test_label_baseline)) {
      kde_fit <- kde1d(samp_dt_baseline_param, mult=bandwidth_mult)
      kde_baseline <- dkde1d(kde_pts, kde_fit)
    }
    
    # Construct KDE comparison plot for the current parameter. 
    plt_curr <- plot_curves_1d_helper(kde_pts, kde_mat, y_new=kde_baseline,
                                      plot_title=plt_title, xlab=param_name, 
                                      ylab="kde")
    plts[[plt_label]] <- plt_curr
  }
  
  if(!is.null(save_dir)) save_plots(plts, "kde1d", save_dir)
  
  return(invisible(plts))
}


get_2d_density_plots <- function(samp_dt, test_labels=NULL, param_types=NULL,
                                 param_names=NULL, itr_start=1L, itr_stop=NULL,
                                 chain_idcs=NULL, save_dir=NULL, 
                                 combine_chains=TRUE) {
  # Uses ggplot2::geom_density_2d() to plot estimated densities between pairs of 
  # parameters from samples. Will produce one plot per pair of parameters within 
  # each value of `test_label`.
  #
  # Args:
  #    samp_dt: existing data.table object storing samples.
  #    combine_chains: logical. Currently only supports TRUE, which groups all
  #                    chains together.
  #    save_dir: character(1), if not NULL, a file path to save the plots to. 
  #    Remaining arguments are used to subset `samp_dt`; see `select_mcmc_samp()`.
  #
  # Returns:
  # list, nested so that the elements of the first level correspond to the 
  # test labels. The sub-lists for each test label then contain the ggplot 
  # objects falling within that test label.

  assert_is_samp_dt(samp_dt)
  
  if(!combine_chains) {
    stop("Currently `get_1d_kde_plots()` only supports `combine_chains = TRUE`.")
  }

  # Determine which plots to create by subsetting rows of `samp_dt`. 
  samp_dt_subset <- select_mcmc_samp(samp_dt, test_labels=test_labels,
                                     param_types=param_types, 
                                     param_names=param_names, itr_start=itr_start,
                                     itr_stop=itr_stop, chain_idcs=chain_idcs)
  
  # Store unique `test_label`-`param_type`-`param_name` combinations.
  id_cols <- c("test_label", "param_type", "param_name")
  plt_id_vars <- unique(samp_dt_subset[, ..id_cols])
  test_lbls <- unique(plt_id_vars$test_label)
  
  # Create plots between pairs of parameters within each test label.
  plts <- list()
  for(lbl in test_lbls) {
    # Get parameters for current test label.
    plt_id_vars_lbl <- plt_id_vars[test_label==lbl]
    plts[[lbl]] <- list()
    n_par <- nrow(plt_id_vars_lbl)
    
    # Subset to current label.
    samp_dt_lbl <- select_mcmc_samp(samp_dt_subset, test_labels=lbl)
    
    for(i in 1:(n_par-1)) {
      par_i <- paste(plt_id_vars_lbl[i,param_type], 
                     plt_id_vars_lbl[i,param_name], sep="-")
      for(j in (i+1):n_par) {
        par_j <- paste(plt_id_vars_lbl[j,param_type], 
                       plt_id_vars_lbl[j,param_name], sep="-")
        plot_tag <- paste(par_i, par_j, sep="_")
        samp_dt_i <- select_mcmc_samp(samp_dt_lbl, 
                                      param_types=plt_id_vars_lbl[i,param_type],
                                      param_names=plt_id_vars_lbl[i,param_name])
        samp_dt_j <- select_mcmc_samp(samp_dt_lbl, 
                                      param_types=plt_id_vars_lbl[j,param_type],
                                      param_names=plt_id_vars_lbl[j,param_name])
        samp_dt_pars <- data.table(par1=samp_dt_i$sample, par2=samp_dt_j$sample)
        plts[[lbl]][[plot_tag]] <- ggplot(samp_dt_pars, aes(x=par1, y=par2)) + 
                                    geom_density_2d() + xlab(par_i) + ylab(par_j)
      }
    }
  }

  if(!is.null(save_dir)) {
    save_plots(unlist(plts,recursive=FALSE), "density2d", save_dir)
  }
  
  return(invisible(plts))
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


get_1d_coverage_plots <- function(samp_dt, test_label_baseline, burn_in_start=NULL, test_labels=NULL,
                                  param_types=NULL, param_names=NULL, xlab="observed", ylab="predicted", 
                                  save_dir=NULL, probs=seq(0.5, 1.0, .1), color_exact="black") {
  # Produces one plot per unique (param type, param name) combination. Each plot will contain one line 
  # per test label that has samples associated with the specific parameter being plotted. These lines 
  # summarize the coverage of the distributions (i.e., the samples for each test label) with respect 
  # to some "baseline" test label, specified by `test_label_baseline` (this typically represents 
  # some notion of the "true" distribution). For a specific plot, let us consider how a single 
  # point is computed: 
  #    (1) The "nominal coverage" interval of probability `p` is estimated from samples for 
  #        the distribution associated with each test label. This interval is centered at the 
  #        empirical median, with the endpoints computed by excluding 100*(1-p)/2 % of the 
  #        samples from each tail, so that the interval contains 100*p % of the samples. 
  #    (2) Let [a,b] be the bounds of this estimated nominal interval for a specific test label. Next we 
  #        estimate the probability of the baseline distribution falling within the interval [a,b]. 
  #        This is accomplished by counting the number of samples of the baseline distribution 
  #        falling within [a,b] and dividing by the number of baseline samples. Let's call the 
  #        computed fraction `q`. 
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
  
  # Probabilities for lower/upper quantiles for nominal coverage interval.
  probs_lower <- (1 - probs) / 2
  probs_upper <- (1 + probs) / 2
  
  # Create one plot per unique (param_type, param_name) combination. 
  plt_list <- list()
  dt_plt <- data.table(param_name=character(), param_type=character(),  
                       test_label=character(), prob=numeric(), 
                       nominal_lower=numeric(), median=numeric(), 
                       nominal_upper=numeric(), N_sample=integer(), actual_coverage=numeric())
  for(i in 1:nrow(plt_id_vars)) {
    # (param_type, param_name) combination for current plot. 
    param_type_curr <- plt_id_vars[i, param_type]
    param_name_curr <- plt_id_vars[i, param_name]
    samp_dt_param <- samp_dt_subset[(param_type==param_type_curr) & (param_name==param_name_curr)]
    lbl_curr <- paste(param_type_curr, param_name_curr, sep="-")
    test_labels_curr <- unique(samp_dt_param$test_label)
    samp_baseline_curr <- samp_dt_baseline[(param_type==param_type_curr) & (param_name==param_name_curr), sample]
    
    # Estimate highest posterior density interval for each test label and probability. 
    dt_plt_param <- data.table(test_label=character(), prob=numeric(), nominal_lower=numeric(),
                               nominal_upper=numeric(), median=numeric(), N_sample=integer())
    for(lbl in test_labels_curr) {
      
      # Nominal coverage for each test label. 
      samp_dt_param_lbl <- samp_dt_param[test_label==lbl, sample]
      median_lbl <- median(samp_dt_param_lbl)
      param_lbl_nominal_interval <- rbind(quantile(samp_dt_param_lbl, probs=probs_lower), 
                                          quantile(samp_dt_param_lbl, probs=probs_upper))
      dt_plt_param_lbl <- data.table(test_label=lbl, prob=probs, 
                                     nominal_lower=quantile(samp_dt_param_lbl, probs=probs_lower), 
                                     nominal_upper=quantile(samp_dt_param_lbl, probs=probs_upper), 
                                     median=median_lbl, N_sample=length(samp_dt_param_lbl))
      
      # Append to current data.table. 
      dt_plt_param <- rbindlist(list(dt_plt_param, dt_plt_param_lbl), use.names=TRUE)
    }
    
    # Compute actual coverage (probability mass of baseline distribution that is captured
    # by the nominal interval).
    # Actual coverage for each test label.
    N_baseline_curr <- length(samp_baseline_curr)
    compute_actual_coverage <- function(l, u) sum((samp_baseline_curr >= l) & (samp_baseline_curr <= u)) / N_baseline_curr
    dt_plt_param[, actual_coverage := compute_actual_coverage(nominal_lower, nominal_upper), by=seq_len(nrow(dt_plt_param))]

    # Generate plot. 
    plt <- ggplot(dt_plt_param) + 
            geom_point(aes(x=prob, y=actual_coverage, color=test_label)) + 
            geom_line(aes(x=prob, y=actual_coverage, color=test_label)) + 
            geom_abline(slope=1, intercept=0, color=color_exact, linetype="dashed") + 
            ggtitle(param_name_curr) + xlab("Nominal Coverage") + ylab("Actual Coverage")
        
    plt_list[[lbl_curr]] <- plt
    
    # Append plot data.
    dt_plt_param[, param_name := param_name_curr]
    dt_plt_param[, param_type := param_type_curr]
    dt_plt <- rbindlist(list(dt_plt, dt_plt_param), use.names=TRUE)
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














