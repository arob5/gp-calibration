# sim_study_functions.r
# 
# Functions related to conducting a large-scale simulation study aimed at 
# testing methods for posterior approximation and sequential design using 
# surrogate models.
#

# ------------------------------------------------------------------------------
# MCMC Processing
# ------------------------------------------------------------------------------

process_mcmc_run <- function(samp_list, rhat_threshold=1.05, 
                             min_itr_threshold=500L, ...) {
  # Post-processing for a single MCMC run (with potentially multiple chains).
  # A chain is defined as "valid" if the within-chain split R-hat values for 
  # all parameters do not exceed `rhat_threshold`. `n_itr` is the total number  
  # of iterations with non-NA values across all valid chains. An MCMC run 
  # (test_label) is marked as valid if it has at least one valid chain and 
  # at least one valid iteration.
  #
  # `samp_list` is a list returned by the MCMC code, with elements "samp", 
  # "info", and "output_list". `min_itr_threshold` is the minimum number of 
  # samples we allow for a single chain. 

  # Check to see if error occurred during run.
  err_occurred <- !is.null(samp_list$output_list[[1]]$condition)
  if(err_occurred) {
    mcmc_summary <- data.table(n_chains=0L, max_rhat=NA, status="mcmc_err")
    return(list(summary=mcmc_summary, chain_info=NULL, rhat=NULL))
  }
  
  if(length(unique(samp_list$samp$test_label)) > 1L) {
    stop("`process_mcmc_run` is defined to operate on a single test label.")
  }
  
  # MCMC samples.
  samp_dt <- samp_list$samp
  
  # Process chain by chain: increase burn-in until each chain individually 
  # satisfies rhat threshold. If a chain fails to meet threshold for any 
  # tried burnin, then `process_mcmc_chain` returns empty data.table.
  chains <- unique(samp_dt$chain_idx)
  chain_list <- lapply(chains, function(i) process_mcmc_chain(samp_dt, i, 
                                                              rhat_threshold=rhat_threshold, 
                                                              min_itr_threshold=min_itr_threshold, ...))
  samp_dt <- rbindlist(chain_list, use.names=TRUE)
  if(nrow(samp_dt) == 0L) {
    mcmc_summary <- data.table(n_chains=0L, max_rhat=NA, status="processing_failed")
    return(list(summary=mcmc_summary, chain_info=NULL, rhat=NULL))
  }
  
  # Parameter level summary: within-chain Rhat statistics.
  id_cols <- c("test_label", "chain_idx", "param_type", "param_name")
  rhat_dt <- calc_R_hat(samp_dt, within_chain=TRUE)$R_hat_vals
  rhat_dt <- rhat_dt[, .SD, .SDcols=c(id_cols, "R_hat")]
  
  # Chain level summary: Maximum within within-chain Rhat over all parameters 
  # within the chain.
  rhat_max_dt <- rhat_dt[, .(rhat=max(R_hat, na.rm=TRUE)), by=.(test_label, chain_idx)]
  other_chain_info <- samp_dt[, .(itr_min=min(itr), itr_max=max(itr)),
                              by=.(test_label, chain_idx)]
  chain_info <- data.table::merge.data.table(rhat_max_dt, other_chain_info,
                                             by=c("test_label", "chain_idx"))

  # Calculate chain weights. Note that any invalid chains have been dropped at 
  # this point, so the weights will sum to one only over the remaining valid
  # chains. First need to subset `info_dt` to align with the adjustments
  # to `samp_dt`.
  info_dt <- samp_list$info
  info_dt <- select_mcmc_samp(info_dt, chain_idcs=chain_info$chain_idx)
  info_list <- list()
  for(i in chain_info$chain_idx) {
    info_list[[i]] <- select_mcmc_samp(info_dt, chain_idcs=i,
                                       itr_start=chain_info[chain_idx==i,itr_min],
                                       itr_stop=chain_info[chain_idx==i,itr_max])
  }
  
  info_dt <- rbindlist(info_list, use.names=TRUE)
  chain_weights <- calc_chain_weights(info_dt)
  chain_info <- data.table::merge.data.table(chain_info, chain_weights,
                                             by=c("test_label", "chain_idx"))
  
  
  # Compute MCMC run summary.
  mcmc_summary <- data.table(n_chains = nrow(chain_info),
                             max_rhat = max(chain_info$rhat),
                             status = "valid")
  
  return(list(summary=mcmc_summary, chain_info=chain_info, rhat=rhat_dt))
}


process_mcmc_chain <- function(samp_dt, chain_idx, rhat_threshold=1.05, 
                               min_itr_threshold=500L, itr_start=1L,
                               n_tries=10L) {
  # `itr_start` can be used to specify an initial burn-in, before computing 
  # any Rhat statistics. `n_tries` is the maximum number of adjustments 
  # that will be made to the burn-in before giving up.

  chain_dt <- select_mcmc_samp(samp_dt, itr_start=itr_start, chain_idcs=chain_idx)
  
  if(nrow(unique(chain_dt[, .(test_label, chain_idx)])) > 1L) {
    stop("`process_mcmc_chain` is defined to operate on a single chain and test label.")
  }
  
  # If the chain passed as argument already doesn't meet the min itr 
  # threshold, return empty data.table.
  if(!chain_meets_min_itr_threshold(chain_dt, min_itr_threshold)) {
    return(get_empty_samp_dt())
  }
  
  # Maximum rhat over all parameters.
  rhat <- calc_R_hat(chain_dt, within_chain=TRUE)$R_hat_vals[,max(R_hat)]
  
  # Compute candidate burn-in values to try sequentially until rhat threshold
  # is met
  itr_cutoffs <- round(seq(chain_dt[,min(itr)], chain_dt[,max(itr)], 
                       length.out=n_tries))
  i <- 2L
  
  # Note that the condition on minimum iteration threshold prevents infinite 
  # loop here, since the final cutoff will be the final iteration and thus 
  # violate the minimum iteration threshold.
  while(isTRUE(rhat > rhat_threshold)) {
    
    # Increase burn-in size.
    chain_dt <- select_mcmc_itr(chain_dt, itr_start=itr_cutoffs[i])
    
    # If the chain doesn't meet the min itr threshold, return empty data.table.
    if(!chain_meets_min_itr_threshold(chain_dt, min_itr_threshold)) {
      return(get_empty_samp_dt())
    }
    
    # Maximum rhat over all parameters.
    rhat <- calc_R_hat(chain_dt, within_chain=TRUE)$R_hat_vals[,max(R_hat)]
    
    i <- i + 1L
  } 
  
  return(chain_dt)
}


chain_meets_min_itr_threshold <- function(chain_dt, min_itr_threshold) {
  if(nrow(unique(chain_dt[, .(test_label, chain_idx)])) > 1L) {
    stop("`process_mcmc_chain` is defined to operate on a single chain and test label.")
  }
  
  n_itr_by_par <- chain_dt[, .N, by=.(param_type, param_name)]
  n_itr_by_par <- unique(n_itr_by_par$N)
  
  if(length(n_itr_by_par) > 1) {
    stop("Different parameters within MCMC chain have different number of iterations.")
  }
  
  n_itr_by_par >= min_itr_threshold
}


# ------------------------------------------------------------------------------
# Loading/Assembling Processed MCMC samples
# ------------------------------------------------------------------------------

get_samp_dt_reps <- function(experiment_dir, round, mcmc_tags=NULL) {
  
}




# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------


get_coverage_plot_reps <- function() {
  
  
}












