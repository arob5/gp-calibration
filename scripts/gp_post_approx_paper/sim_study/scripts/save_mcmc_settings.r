# 
# save_mcmc_settings.r
# Defines the approximate MCMC algorithms that will be run and saves the 
# settings to file. The saved CSV file is read in by `run_approx_mcmc.r`. 
#
# Andrew Roberts
#

library(data.table)

# Output directory and path for saving settings file.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
out_dir <- file.path(base_dir, "output", "gp_inv_prob", "vsem")
out_path <- file.path(out_dir, "mcmc_approx_settings.rds")

# Common settings that will be applied to all MCMC algorithms.
common_settings <- list(n_itr=100000L, try_parallel=FALSE, n_chain=4L)

# Common settings that will be applied to all algorithms using the BayesianTools
# wrapper function.
common_bt_settings <- list(mcmc_func_name="mcmc_bt_wrapper", sampler="DEzs", 
                           defer_ic=TRUE, settings_list=list(consoleUpdates=25000))

# List of MCMC settings.
mcmc_settings_list <- list(
  c(list(mcmc_tag="quantile7", approx_type="quantile", alpha=0.7),
    common_settings, common_bt_settings),
  c(list(mcmc_tag="quantile8", approx_type="quantile", alpha=0.8),
    common_settings, common_bt_settings),
  c(list(mcmc_tag="quantile9", approx_type="quantile", alpha=0.9),
    common_settings, common_bt_settings),
  c(list(mcmc_tag="mean", approx_type="mean"),
    common_settings, common_bt_settings),
  c(list(mcmc_tag="marginal", approx_type="marginal"),
    common_settings, common_bt_settings),
  c(list(mcmc_tag="mcwmh-joint", mcmc_func_name="mcmc_noisy_llik", use_joint=TRUE),
    common_settings),
  c(list(mcmc_tag="mcwmh-ind", mcmc_func_name="mcmc_noisy_llik", use_joint=FALSE),
    common_settings),
  c(list(mcmc_tag="marg-acc-prob-joint", mcmc_func_name="mcmc_gp_acc_prob_approx", 
       approx_type="joint-marginal"), common_settings),
  c(list(mcmc_tag="marg-acc-prob-ind", mcmc_func_name="mcmc_gp_acc_prob_approx", 
       approx_type="marginal"), common_settings)
)

# Set names of list to the MCMC tags.
tags <- sapply(mcmc_settings_list, function(x) x$mcmc_tag)
names(mcmc_settings_list) <- tags

# Save settings to file.
saveRDS(mcmc_settings_list, out_path)





