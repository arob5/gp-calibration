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

library(lhs)
library(ggplot2)
library(viridis)
library(parallel)
library(gridExtra)
library(data.table)

base_dir <- getwd()
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

# Settings for MCMC runs.  
N_mcmc <- 50000

# Forward model parameters. 
N_time_step <- 365*3

# Global likelihood normalization settings. 
default_conditional <- FALSE
default_normalize <- TRUE


# -----------------------------------------------------------------------------
# Bayesian inverse problem setup. 
# -----------------------------------------------------------------------------

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
inv_prob <- list(y=y, fwd=fwd, fwd_vect=fwd_vect, par_prior=par_prior, cov_obs=cov_obs)

# Plot ground truth and observations. 
df_inv_prob_data <- data.frame(t=times, y_true=output_true, y=y)
plt_outputs <- ggplot(data=df_inv_prob_data) +
               geom_point(aes(x=t, y=y), color="black") + 
               geom_point(aes(x=t, y=y_true), color="red") + 
               ylab("output") + ggtitle("True and Observed y")

plot(plt_outputs)
ggsave(file.path(output_dir, "output_data.png"))


# -----------------------------------------------------------------------------
# Inversion using exact forward model. 
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# Initial design and fit emulators. 
# -----------------------------------------------------------------------------

seed_init_design <- 500
set.seed(seed_init_design)

# Generate design. 
N_design <- 10 * dim_par
design_info <- list()
design_info$input <- get_batch_design("LHS", N_design, prior_params=inv_prob$par_prior)
design_info$fwd <- llik_exact$run_fwd_model(design_info$input)
design_info$llik <- llik_exact$assemble_llik(design_info$input)
design_info$lprior <- calc_lprior_theta(design_info$input, inv_prob$par_prior)

# Test points. 
N_grid <- 5000
u_grid <- get_batch_design("LHS", N_grid, prior_params=inv_prob$par_prior)
test_info <- list(input=u_grid)
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

# Forward model emulator.
em_fwd_gp <- gpWrapperHet(design_info$input, design_info$fwd, normalize_output=TRUE, scale_input=TRUE)
em_fwd_gp$fit("Gaussian", "constant", estimate_nugget=FALSE)
em_fwd <- llikEmulatorFwdGaussDiag("em_fwd", em_fwd_gp, t(inv_prob$y), sig2=diag(cov_obs), 
                                   default_conditional=default_conditional, 
                                   default_normalize=default_normalize, use_fixed_lik_par=TRUE, 
                                   par_names=rownames(inv_prob$par_prior))
emulator_pred_fwd <- em_fwd$predict_emulator(test_info$input, return_cov=TRUE)


#
# Log-likelihood predictions. 
#

em_llik$plot_pred_validation(test_info$input, true_llik=test_info$llik, 
                             emulator_pred_list=emulator_pred_llik, plot_type="llik",
                             conditional=default_conditional, normalize=default_normalize,
                             include_interval=TRUE, interval_method="pm_std_dev", N_std_dev=1, 
                             plot_title="Log-Likelihood Predictions [llik]")


em_fwd$plot_pred_validation(test_info$input, true_llik=test_info$llik, 
                            emulator_pred_list=emulator_pred_fwd, plot_type="llik",
                            conditional=default_conditional, normalize=default_normalize,
                            include_interval=TRUE, interval_method="pm_std_dev", N_std_dev=1, 
                            plot_title="Log-Likelihood Predictions [fwd]")









