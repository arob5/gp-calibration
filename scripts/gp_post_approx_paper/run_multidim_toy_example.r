#
# run_multidim_toy_example.r
# Bayesian inversion for synthetic examples with multidimensional 
# parameter spaces. This file is written to be able to run on a 
# computing cluster. 
#
# Andrew Roberts
#

# -----------------------------------------------------------------------------
# Setup 
# -----------------------------------------------------------------------------

print("---------------------------------- Setup ----------------------------------")

library(lhs)
library(ggplot2)
library(viridis)
library(parallel)
library(gridExtra)
library(data.table)

base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
src_dir <- file.path(base_dir, "src")
output_dir <- file.path(base_dir, "output", "gp_post_approx_paper", "multidim_toy_examples")

source(file.path(src_dir, "gpWrapper.r"))
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "statistical_helper_functions.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(src_dir, "mcmc_helper_functions.r"))
source(file.path(src_dir, "llikEmulator.r"))
source(file.path(src_dir, "gp_emulator_functions.r"))
source(file.path(src_dir, "sim_study_functions.r"))
source(file.path(src_dir, "mcmc_calibration_functions.r"))
source(file.path(src_dir, "seq_design_gp.r"))
source(file.path(src_dir, "seq_design.r"))
source(file.path(src_dir, "gp_mcmc_functions.r"))

# Seeds for random number generator for different portions of analysis. 
seed_data <- 82
seed_init_design <- 500

# Global likelihood normalization settings. 
default_conditional <- FALSE
default_normalize <- TRUE

# GP-accelerated MCMC algorithms to run. 
mcmc_tags <- c("gp-mean", "gp-marg", "mcwmh-joint", "mcwmh-ind", "acc-prob-marg")

# Settings for MCMC runs.  
N_mcmc <- 50000

# If `output_dir` already exists, append timestep and create new directory.
# Otherwise create the directory. 
if(dir.exists(output_dir)) {
  timestamp <- as.character(Sys.time())
  output_dir <- paste(output_dir, timestamp, sep="_")
  message("`output_dir` already exists. Creating new `output_dir:`")
  message(output_dir)
  dir.create(output_dir)
} else {
  dir.create(output_dir, recursive=TRUE)
}


# -----------------------------------------------------------------------------
# Bayesian inverse problem setup. 
# -----------------------------------------------------------------------------

print("--------------------- Bayesian Inverse Problem Setup -------------------------")

set.seed(seed_data)

#
# Forward model: linear combination of basis functions. 
#

# Basis functions. 
basis_func_list <- list(phi0=function(t) rep(1, length(t)), phi1=function(t) t, 
                        phi2=function(t) cos(2*pi*t), phi3=function(t) cos(2*pi*2*t), 
                        phi4=function(t) sin(2*pi*t), phi5=function(t) sin(2*pi*t/3))
dim_par <- length(basis_func_list)

# Discretize the unit time interval to define the dimension of the output space. 
dim_output <- 10
times <- seq(0, 1, length.out=dim_output)

# Forward model. 
eval_with_args <- function(FUN, ...) FUN(...)
Phi <- mapply(eval_with_args, basis_func_list, list(times))

fwd <- function(par) {
  Phi %*% matrix(par, ncol=1)
}

fwd_vect <- function(par_mat) {
  # `par_mat` assumed to have dimension (num inputs, dim_par). Returns 
  # output of dimension (num inouts, dim_output). 
  par_mat %*% t(Phi)
}

#
# Likelihood. 
#

# Observation covariance. 
cov_obs <- diag(0.5, nrow=dim_output)

#
# Prior. 
#

# Prior distribution. 
par_prior <- data.frame(dist="Gaussian", param1=rep(1,dim_par), param2=rep(1,dim_par), 
                        row.names=paste0("u", 1:dim_par))

#
# Synthetic data: ground truth and noisy observations. 
# 

# Ground truth. 
par_true <- sample_prior_theta(par_prior)
output_true <- fwd(par_true)

# Noisy observations. 
y <- output_true + t(chol(cov_obs)) %*% matrix(rnorm(dim_output), ncol=1)

#
# Collect inverse problem components in list, save plots summarizing data.  
#

# List containing inverse problem components. 
inv_prob <- list(y=y, fwd=fwd, fwd_vect=fwd_vect, par_prior=par_prior, cov_obs=cov_obs, 
                 par_true=par_true, output_true=output_true)

# Plot ground truth and observations. 
df_inv_prob_data <- data.frame(t=times, y_true=output_true, y=y)
plt_outputs <- ggplot(data=df_inv_prob_data) +
               geom_point(aes(x=t, y=y), color="black") + 
               geom_point(aes(x=t, y=y_true), color="red") + 
               ylab("output") + ggtitle("True and Observed y")

plot(plt_outputs)
ggsave(file.path(output_dir, "output_data.png"), plt_outputs)


# -----------------------------------------------------------------------------
# Inversion using exact forward model. 
# -----------------------------------------------------------------------------

print("--------------------- Bayesian Inversion with True Forward Model -------------------------")

# Store exact log-likelihood object. 
llik_exact <- llikEmulatorExactGaussDiag(llik_lbl="exact", fwd_model=inv_prob$fwd, 
                                         fwd_model_vectorized=inv_prob$fwd_vect,
                                         y_obs=t(inv_prob$y), dim_par=dim_par,
                                         use_fixed_lik_par=TRUE, sig2=diag(cov_obs),
                                         par_names=rownames(inv_prob$par_prior), 
                                         default_conditional=default_conditional, 
                                         default_normalize=default_normalize)
inv_prob$llik_obj <- llik_exact
save(inv_prob, file=file.path(output_dir, "inv_prob_list.RData"))

# MCMC sampling using exact likelihood. 
mcmc_par_init <- sample_prior_theta(inv_prob$par_prior)
mcmc_exact_list <- mcmc_gp_noisy(inv_prob$llik_obj, inv_prob$par_prior, N_itr=N_mcmc, 
                                 mode="MCMH", par_init=mcmc_par_init)
samp_dt <- format_mcmc_output(mcmc_exact_list$samp, test_label="exact")
mcmc_list <- list(exact=mcmc_exact_list[setdiff(names(mcmc_exact_list), "samp")])


# -----------------------------------------------------------------------------
# Initial design and fit emulators. 
# -----------------------------------------------------------------------------

print("--------------------- Fitting GP Emulators -------------------------")

seed_init_design <- 500
set.seed(seed_init_design)

# Generate design. 
N_design <- 10 * dim_par
design_method <- "LHS"
design_info <- list(design_method=design_method, N_design=N_design, seed=seed_init_design)
design_info$input <- get_batch_design("LHS", N_design, prior_params=inv_prob$par_prior)
design_info$fwd <- llik_exact$run_fwd_model(design_info$input)
design_info$llik <- llik_exact$assemble_llik(design_info$input)
design_info$lprior <- calc_lprior_theta(design_info$input, inv_prob$par_prior)
save(design_info, file=file.path(output_dir, "design_info.RData"))

# Test points. 
N_design_test <- 5000
test_design_method <- "LHS"
u_grid <- get_batch_design(test_design_method, N_design_test, prior_params=inv_prob$par_prior)
test_info <- list(input=u_grid, design_method=test_design_method, N_design=N_design_test)
test_info$fwd <- llik_exact$run_fwd_model(test_info$input) 
test_info$llik <- llik_exact$assemble_llik(test_info$input)
test_info$lprior <- calc_lprior_theta(test_info$input, inv_prob$par_prior)

# Log-likelihood emulator. 
em_llik_gp <- gpWrapperHet(design_info$input, matrix(design_info$llik, ncol=1), normalize_output=TRUE, scale_input=TRUE)
em_llik_gp$fit("Gaussian", "constant", estimate_nugget=FALSE)
em_llik <- llikEmulatorGP("em_llik", em_llik_gp, default_conditional=default_conditional, 
                          default_normalize=default_normalize, lik_par=diag(cov_obs), 
                          use_fixed_lik_par=TRUE)
emulator_pred_llik <- em_llik$predict_emulator(test_info$input, return_cov=TRUE)
test_info$pred_llik <- emulator_pred_llik

# Forward model emulator.
em_fwd_gp <- gpWrapperHet(design_info$input, design_info$fwd, normalize_output=TRUE, scale_input=TRUE)
em_fwd_gp$fit("Gaussian", "constant", estimate_nugget=FALSE)
em_fwd <- llikEmulatorFwdGaussDiag("em_fwd", em_fwd_gp, t(inv_prob$y), sig2=diag(cov_obs), 
                                   default_conditional=default_conditional, 
                                   default_normalize=default_normalize, use_fixed_lik_par=TRUE, 
                                   par_names=rownames(inv_prob$par_prior))
emulator_pred_fwd <- em_fwd$predict_emulator(test_info$input, return_cov=TRUE)
test_info$pred_fwd <- emulator_pred_fwd

# Save test info with emuilator predictions. 
save(test_info, file=file.path(output_dir, "test_info.RData"))


#
# Log-likelihood predictions. 
#

em_llik_pred_plot <- em_llik$plot_pred_validation(test_info$input, true_llik=test_info$llik, 
                                                  emulator_pred_list=emulator_pred_llik, plot_type="llik",
                                                  conditional=default_conditional, normalize=default_normalize,
                                                  include_interval=TRUE, interval_method="pm_std_dev", N_std_dev=1, 
                                                  plot_title="Log-Likelihood Predictions [llik]")
em_fwd_pred_plot <- em_fwd$plot_pred_validation(test_info$input, true_llik=test_info$llik, 
                                                emulator_pred_list=emulator_pred_fwd, plot_type="llik",
                                                conditional=default_conditional, normalize=default_normalize,
                                                include_interval=TRUE, interval_method="pm_std_dev", N_std_dev=1, 
                                                plot_title="Log-Likelihood Predictions [fwd]")

ggsave(file.path(output_dir, "llik_em_pred_plot_llik.png"), em_llik_pred_plot)
ggsave(file.path(output_dir, "llik_em_pred_plot_fwd.png"), em_fwd_pred_plot)


# -----------------------------------------------------------------------------
# Run MCMC Samplers.  
# -----------------------------------------------------------------------------

print("--------------------- Running GP-Accelerated MCMC: llik emulation -------------------------")

mcmc_tags <- "gp-mean"

mcmc_info_list_llik <- run_approx_mcmc_comparison(inv_prob_list=inv_prob, llik_em_obj=em_llik, 
                                                  mcmc_tags=mcmc_tags, save_dir=output_dir, 
                                                  samp_dt=samp_dt, mcmc_list=mcmc_list, 
                                                  test_label_suffix="llik", N_itr=N_mcmc, 
                                                  par_init=mcmc_par_init, overwrite=TRUE)


print("--------------------- Running GP-Accelerated MCMC: forward model emulation -------------------------")

mcmc_info_list_fwd <- run_approx_mcmc_comparison(inv_prob_list=inv_prob, llik_em_obj=em_fwd, 
                                                 mcmc_tags=mcmc_tags, save_dir=output_dir, 
                                                 samp_dt=samp_dt, mcmc_list=mcmc_list, 
                                                 test_label_suffix="fwd", N_itr=N_mcmc, 
                                                 par_init=mcmc_par_init, overwrite=TRUE)


print("-------------------------- End of Script ------------------------------")




