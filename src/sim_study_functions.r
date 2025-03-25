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

get_mcmc_ids <- function(experiment_dir, round, mcmc_tag, em_tag, design_tag, 
                         only_valid=TRUE) {
  # Assembles an ID map table with columns "mcmc_id", "em_id", "em_tag",
  # "design_id", and "design_tag" containing the set of MCMC IDs associated
  # with the specified experiment, round, MCMC tag, emulator tag, and design
  # tag specified in the function arguments. If `only_valid` is TRUE, then 
  # the MCMC iterations and chains are subset using the MCMC post-processing
  # results in the "summary_files" directory. This also implies that MCMC runs
  # with no valid chains will be dropped entirely. Setting `only_valid=TRUE`
  # also attaches additional information as columns in the returned data.table,
  # including Rhat values and chain weights.
  #
  # NOTE: 
  # Currently this function returns a data.table that is unique by `mcmc_id`
  # if `only_valid=FALSE` and is unique by (`mcmc_id`, `chain_idx`) if 
  # `only_valid=TRUE`. This is a bit confusing and should probably be changed.
  
  # Assemble set of MCMC tags that are associated with the specified experiment,
  # round, MCMC tag, and emulator tag.
  round_tag <- paste0("round", round)
  mcmc_dir <- file.path(experiment_dir, round_tag, "mcmc", mcmc_tag)
  mcmc_id_map <- fread(file.path(mcmc_dir, "id_map.csv"))
  
  # Subset MCMC IDs to those associated with specified emulator tag.
  em_tag_curr <- em_tag
  mcmc_id_map <- mcmc_id_map[em_tag==em_tag_curr]
  
  # Subset MCMC IDs to those associated with specified design tag.
  design_tag_curr <- design_tag
  em_id_map <- fread(file.path(experiment_dir, round_tag, "em", em_tag, "id_map.csv"))
  mcmc_id_map <- data.table::merge.data.table(mcmc_id_map, em_id_map, by="em_id",
                                              all.x=TRUE, all.y=FALSE)
  mcmc_id_map <- mcmc_id_map[design_tag==design_tag_curr]
  
  # Use MCMC post-processing results to select only valid runs/chains/itrs.
  if(only_valid) {
    
    # Retain only valid runs.
    mcmc_tag_curr <- mcmc_tag
    summary_dir <- file.path(experiment_dir, round_tag, "mcmc", "summary_files")
    run_summary <- fread(file.path(summary_dir, "mcmc_summary.csv"))
    run_summary <- run_summary[(status=="valid") & (mcmc_tag == mcmc_tag_curr)]
    mcmc_id_map <- data.table::merge.data.table(mcmc_id_map, run_summary, 
                                                all=FALSE, by="mcmc_id")
    
    # Retain only valid chains, and set burn-ins determined by preprocessing.
    chain_summary <- fread(file.path(summary_dir, "chain_summary.csv"))
    chain_summary <- chain_summary[mcmc_tag==mcmc_tag_curr]
    mcmc_id_map <- data.table::merge.data.table(mcmc_id_map, chain_summary,
                                                all=FALSE, by=c("mcmc_id", "mcmc_tag"))
  }
  
  return(mcmc_id_map)
}


get_samp_dt_reps <- function(experiment_dir, round, mcmc_tag, em_tag, design_tag, 
                             only_valid=TRUE) {
  # Loads all MCMC samples found within the given experiment, round, MCMC tag,
  # emulator tag, and design tag. If `only_valid = TRUE`, drops invalid 
  # runs/chains/iterations - see `get_mcmc_ids()` for details. All MCMC runs
  # within the given set of tags are viewed as random replications of a single
  # experimental setup (e.g., replications stemming from different initial 
  # design samples). All of these are compiled into a single data.table, with 
  # a `rep_id` column added, which is set to the corresponding design_id for
  # each run.

  mcmc_id_dt <- get_mcmc_ids(experiment_dir, round, mcmc_tag, em_tag, design_tag, 
                             only_valid=TRUE)
  
  mcmc_ids <- unique(mcmc_id_dt$mcmc_id)
  mcmc_tag_dir <- file.path(experiment_dir, paste0("round", round), "mcmc",
                            mcmc_tag)
  
  # Read MCMC samples, and subset chains/itrs based on `mcmc_ids`.
  dt_list <- vector(mode="list", length=length(mcmc_ids))

  for(i in seq_along(dt_list)) {
    id <- mcmc_ids[i]
    design_id <- mcmc_id_dt[mcmc_id==id, design_id][1]
    samp_list <- readRDS(file.path(mcmc_tag_dir, id, "samp.rds"))
    samp_dt <- samp_list$samp
    
    # Subset to include only valid iterations.
    if(only_valid) {
      chain_itr_dt <- mcmc_id_dt[mcmc_id==id, .(chain_idx, itr_start=itr_min, 
                                                itr_stop=itr_max)]
      samp_dt <- select_itr_by_chain(samp_dt, chain_itr_dt)
    }
    
    samp_dt[, rep_id := design_id]
    dt_list[[i]] <- samp_dt
  }
  
  samp_dt_reps <- rbindlist(dt_list, use.names=TRUE)
  return(list(samp=samp_dt_reps, ids=mcmc_id_dt))
}


get_samp_dt_reps_agg <- function(experiment_dir, round, mcmc_tag, em_tag, 
                                 design_tag, only_valid=TRUE, format_long=FALSE,
                                 interval_probs=NULL) {
  # This is similar to `get_samp_dt_reps()`, but here each MCMC run is 
  # aggregated after it is loaded, producing one-dimensional summaries of each
  # variable. Currently, this includes mean, variance, and coverage. If such 
  # summaries are all that is needed, this function is a much more 
  # space-efficient option compared to `get_samp_dt_reps()`.
  #
  # TODO:
  # 1.) add support for using chain weights.
  
  # Define grouping columns for aggregation.
  group_cols <- c("test_label", "param_type", "param_name")
  
  # Determine set of MCMC runs to read.  
  mcmc_id_dt <- get_mcmc_ids(experiment_dir, round, mcmc_tag, em_tag, design_tag, 
                             only_valid=TRUE)
  
  mcmc_ids <- unique(mcmc_id_dt$mcmc_id)
  mcmc_tag_dir <- file.path(experiment_dir, paste0("round", round), "mcmc",
                            mcmc_tag)
  
  # Separate data.tables will store means/vars and marginal credible intervals.
  stats_list <- vector(mode="list", length=length(mcmc_ids))
  cred_interval_list <- vector(mode="list", length=length(mcmc_ids))
  
  # Read MCMC samples, and subset chains/itrs based on `mcmc_ids`.
  for(i in seq_along(stats_list)) {
    
    # Read samples from run.
    id <- mcmc_ids[i]
    design_id <- mcmc_id_dt[mcmc_id==id, design_id][1]
    samp_list <- readRDS(file.path(mcmc_tag_dir, id, "samp.rds"))
    samp_dt <- samp_list$samp
    
    # Subset to include only valid iterations.
    if(only_valid) {
      chain_itr_dt <- mcmc_id_dt[mcmc_id==id, .(chain_idx, itr_start=itr_min, 
                                                itr_stop=itr_max)]
      samp_dt <- select_itr_by_chain(samp_dt, chain_itr_dt)
    }
    
    # Compute univariate aggregate statistics (mean, variance).
    stat_info <- compute_mcmc_param_stats(samp_dt, subset_samp=FALSE, 
                                          format_long=format_long,
                                          group_cols=group_cols,
                                          interval_probs=interval_probs)
    
    # Means/variances.
    par_stats <- stat_info$par_stats
    par_stats[, rep_id := design_id]
    stats_list[[i]] <- par_stats
    
    # Marginal credible intervals.
    cred_intervals <- stat_info$cred_intervals
    cred_intervals[, rep_id := design_id]
    cred_interval_list[[i]] <- cred_intervals
  }
  
  par_stats_dt <- rbindlist(stats_list, use.names=TRUE)
  cred_intervals_dt <- rbindlist(cred_interval_list, use.names=TRUE)
  
  return(list(par_stats=par_stats_dt, cred_intervals=cred_intervals_dt,
              ids=mcmc_id_dt))
}


# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------


get_coverage_plot_reps <- function() {
  .NotYetImplemented()
}












