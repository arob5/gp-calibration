#
# run_vsem_example.r
# Bayesian inversion for Very Simple Ecosystem Model (VSEM). This file 
# is written to be able to run on a computing cluster. 
#
# Andrew Roberts
#

# TODO:
#    - Need to check llikEmulatorGP and forward emulator and ensure that the 
#      normalization/conditional parameters are properly defined and resulting 
#      in the correct posterior approximations.
#    - Write a function that acts on `par_prior` and returns bounds for each parameter. 
#      Will be -Inf/Inf for unbounded. 
#    - Talk with Mike about reasonable priors to use in this experiment. 
#      e.g., VSEM default parameter range for LAR is enormous, and would
#      probably benefit from something like a log-normal/half-normal/Gamma 
#      prior rather than a uniform one. 
#    - Set up the settings here so that the VSEM parametrer names can be passed
#      in from the bash file; right now `dim_par` and `dim_output` do nothing. 
#    - Add hyperparameter constraints on other kergp models (not just Gaussian+Quadratic). 
#    - Add test points along coordinate axes.
#    - Add basis function forward model emulator.
#    - Figure out weird sampling behavior for approximate posteriors. 

# -----------------------------------------------------------------------------
# docopt string for parsing command line arguments.  
# -----------------------------------------------------------------------------

"Usage:
  test_docopt.r <run_id> <output_dir> [options]
  test_docopt.r (-h | --help)

Options:
  -h --help                                 Show this screen.
  --dim_par=<dim_par>                       Dimension of parameter space. 
  --dim_output=<dim_output>                 Dimension of output space. 
  --N_design=<N_design>                     Number of design points.
  --design_method=<design_method>           Algorithm to generate design points. 
  --N_design_test=<N_design_test>           Number of validation points for emulators. 
  --design_method_test=<design_method_test> Algorithm to generate validation points. 
  --mcmc_tags=<mcmc_tag,mcmc_tag>           GP-approx MCMC algorithms. 
  --N_mcmc=<N_mcmc>                         Number of MCMC iterations. 
" -> doc


# -----------------------------------------------------------------------------
# Setup 
# -----------------------------------------------------------------------------

print("---------------------------------- Setup ----------------------------------")

run_approx_mcmc <- FALSE

# Base directory (i.e., project directory). Required for loading the R project (renv). 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")

# Load R project. 
renv::load(base_dir)
print(".libPaths()")
print(.libPaths())
renv::status()

library(lhs)
library(ggplot2)
library(viridis)
library(parallel)
library(gridExtra)
library(data.table)
library(docopt)

# Read command line arguments. 
arguments <- docopt(doc)
run_id <- arguments$run_id
output_dir <- arguments$output_dir
required_settings <- c("dim_par", "dim_output", "N_design", "design_method", 
                       "N_design_test", "design_method_test", "N_mcmc")
settings <- arguments[required_settings]

# Source and output paths.  
src_dir <- file.path(base_dir, "src")

# Source scripts. 
source(file.path(src_dir, "gp_helper_functions.r"))
source(file.path(src_dir, "gpWrapper.r"))
source(file.path(src_dir, "inv_prob_test_functions.r"))
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "statistical_helper_functions.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(src_dir, "llikEmulator.r"))
source(file.path(src_dir, "gp_emulator_functions.r"))
source(file.path(src_dir, "mcmc_helper_functions.r"))
source(file.path(src_dir, "sim_study_functions.r"))
source(file.path(src_dir, "mcmc_calibration_functions.r"))
source(file.path(src_dir, "seq_design.r"))
source(file.path(src_dir, "seq_design_gp.r"))
source(file.path(src_dir, "seq_design_for_post_approx.r"))
source(file.path(src_dir, "gp_mcmc_functions.r"))

# Seeds for random number generator for different portions of analysis. 
seed_data <- 623434
seed_design <- 26423

# Number of parameters to calibrate. Will determine the number of basis 
# functions in the forward model. 
if(!is.null(settings$dim_par)) settings$dim_par <- as.integer(settings$dim_par)

# Number of forward model outputs. For forward model emulation, this means 
# one GP must be fit per output. 
if(is.null(settings$dim_output)) {settings$dim_output <- 10L
} else {settings$dim_output <- as.integer(settings$dim_output)}


# Design settings for emulator construction. 
if(is.null(settings$N_design)) {
  settings$N_design <- 10L*settings$dim_par
} else {settings$N_design <- as.integer(settings$N_design)}

if(is.null(settings$design_method)) settings$design_method <- "LHS"

# Design settings for emulator validation. 
if(is.null(settings$N_design_test)) {
  settings$N_design_test <- min(100L*settings$dim_par, 500L)
} else {settings$N_design_test <- as.integer(settings$N_design_test)}

if(is.null(settings$design_method_test)) settings$design_method_test <- "LHS"

# Global likelihood normalization settings. 
default_conditional <- FALSE
default_normalize <- TRUE

# GP-accelerated MCMC algorithms to run. 
if(is.null(settings$mcmc_tags)) {
  settings$mcmc_tags <- c("gp-mean", "gp-marg", "gp-quantile", "mcwmh-joint", "mcwmh-ind")
} else {
  settings$mcmc_tags <- unlist(strsplit(settings$mcmc_tags, ","))
}

# Settings for MCMC runs.  
if(is.null(settings$N_mcmc)) {settings$N_mcmc <- 50000L
} else {settings$N_mcmc <- as.integer(settings$N_mcmc)}

# Create global variables with the setting values. 
missing_settings <- names(settings)[sapply(settings, is.null)]
if(length(missing_settings) > 0) stop("Missing settings: ", missing_settings)
dim_par <- settings$dim_par
dim_output <- settings$dim_output
N_design <- settings$N_design
design_method <- settings$design_method
N_design_test <- settings$N_design_test
design_method_test <- settings$design_method_test
mcmc_tags <- settings$mcmc_tags
N_mcmc <- settings$N_mcmc


# If `output_dir` already exists, append timestep and create new directory.
# Otherwise create the directory. 
if(dir.exists(output_dir)) {
  stop("`output_dir` already exists.")
} else {
  dir.create(output_dir, recursive=TRUE)
}
save(settings, file=file.path(output_dir, "settings.RData"))

# Print settings. 
print("--------------------- User specified settings -------------------------")

print("--> General settings")
print(paste0("run ID: ", run_id))
print(paste0("output directory: ", output_dir))
print(paste0("Data generation seed (seed_data): ", seed_data))
print(paste0("Design data seed (seed_design): ", seed_design))
print(paste0("Likelihood densities conditional: ", default_conditional))
print(paste0("Likelihood densities normalized: ", default_normalize))

print("--> Inverse problem")
print(paste0("Parameter space dimension: ", dim_par))
print(paste0("Forward model output dimension: ", dim_output))

print("--> Design")
print(paste0("Number design points: ", N_design))
print(paste0("Design method: ", design_method))
print(paste0("Number test points: ", N_design_test))
print(paste0("Test design method: ", design_method_test))

print("--> MCMC settings")
print(paste0("Number iterations: ", N_mcmc))
print(paste0("GP-accelerated MCMC algs: ", paste(mcmc_tags, collapse=", ")))


# -----------------------------------------------------------------------------
# Bayesian inverse problem setup. 
# -----------------------------------------------------------------------------

print("--------------------- Bayesian Inverse Problem Setup -------------------------")

set.seed(seed_data)

#
# Forward model. 
#

# VSEM default parameter values.
par_names <- get_vsem_par_names()
par_default <- get_vsem_default_pars()
print("Parameter defaults:")
print(data.frame(par_name=par_names, default_value=par_default))

# Parameters to calibrate.
par_cal_names <- c("KEXT", "LUE", "tauV")
dim_par <- length(par_cal_names)

# Define prior on calibration parameter.
par_cal_idx <- which(par_names %in% par_cal_names)
par_prior_params <- get_vsem_default_priors()[par_cal_idx,,drop=FALSE]
rownames(par_prior_params) <- par_prior_params$par_name
print("Prior on calibration parameter:")
print(par_prior_params)

# Number of time points (days) and model driver. 
n_year <- 3
n_time <- 365*n_year
time_points <- 1:n_time
driver <- BayesianTools::VSEMcreatePAR(days=time_points)
driver_plt <- ggplot(data.frame(t=time_points,PAR=driver)) + 
                geom_point(aes(x=time_points,y=PAR)) + 
                ggtitle("Model Driver (PAR)") + 
                xlab("Day") + ylab("PAR")
ggsave(file.path(output_dir, "vsem_driver_data.png"), driver_plt)

# Map from calibration parameter to VSEM outputs.
output_names <- get_vsem_output_names()
lai_idx <- which(output_names=="LAI") 
param_to_output_map <- get_vsem_fwd_model(driver, dim_par, par_cal_idx, par_default, simplify=FALSE)

# Observation operator (map from model outputs to observable): daily LAI observations. 
obs_op <- function(model_outputs) {
  single_run <- (dim(model_outputs)[1]==1L)
  lai_trajectory <- model_outputs[,,lai_idx]
  if(single_run) lai_trajectory <- matrix(lai_trajectory, nrow=1L)
  lai_trajectory
}

# Forward model: defined here as the calibration parameter-to-observable map. Returns 
# matrix of dimension (N_run,N_time).
fwd_model <- function(par_cal) {
  obs_op(param_to_output_map(par_cal))
}

#
# Likelihood, ground truth, and observed data. 
#

signal_to_noise_ratio <- 20
par_true <- par_default
par_cal_true <- par_true[par_cal_idx]
model_outputs_true <- param_to_output_map(par_cal_true)
y_true <- obs_op(model_outputs_true)
sig_eps <- mean(y_true) / signal_to_noise_ratio
sig2_eps <- sig_eps^2
sig2_true <- sig2_eps # Currently obs variance is fixed, 
                      # so no distinction between exact and estimated variance. 
y <- y_true + sig_eps*rnorm(n=length(drop(y_true)))

# Plot ground truth and observed data.
df_model_outputs_true <- data.frame(model_outputs_true[1,,])
colnames(df_model_outputs_true) <- output_names
df_model_outputs_true$time <- time_points
for(output_var in output_names) {
  output_var <- sym(output_var)
  plt <- ggplot(df_model_outputs_true) + geom_line(aes(x=time,y=!!output_var)) + 
         ggtitle("Ground Truth Trajectory") + xlab("Days") + ylab(output_var)
  ggsave(file.path(output_dir, paste0("true_",output_var,".png")), plt)
}

df_model_outputs_true$y <- drop(y)
df_model_outputs_true$y_true <- drop(y_true)
plt_obs <- ggplot(df_model_outputs_true) + 
            geom_line(aes(x=time, y=y), color="red") + 
            geom_line(aes(x=time, y=y_true)) +
            ggtitle("Observable: true signal vs. obs") + xlab("Days") + ylab("LAI")
ggsave(file.path(output_dir, "true_vs_obs.png"), plt_obs)


#
# Collect inverse problem components in list, save plots summarizing data.

# List containing inverse problem components.
inv_prob <- list(y=y, fwd=fwd_model, fwd_vect=fwd_model, par_prior=par_prior_params,
                 sig2_eps=sig2_eps, par_true=par_cal_true, output_true=y_true,
                 times=time_points, par_names=rownames(par_prior_params), 
                 dim_par=dim_par, param_to_output_map=param_to_output_map, 
                 obs_op=obs_op)

# -----------------------------------------------------------------------------
# Exact log-likelihood object. 
# -----------------------------------------------------------------------------

print("--------------------- Creating Exact Log-Likelihood Object -------------------------")

# Store exact log-likelihood object. 
llik_exact <- llikEmulatorExactGaussDiag(llik_lbl="exact", fwd_model=inv_prob$fwd, 
                                         fwd_model_vectorized=inv_prob$fwd_vect,
                                         y_obs=inv_prob$y, dim_par=as.integer(dim_par),
                                         use_fixed_lik_par=TRUE, sig2=inv_prob$sig2_eps,
                                         par_names=rownames(inv_prob$par_prior), 
                                         default_conditional=default_conditional, 
                                         default_normalize=default_normalize)
inv_prob$llik_obj <- llik_exact
save(inv_prob, file=file.path(output_dir, "inv_prob_list.RData"))


# -----------------------------------------------------------------------------
# Initial design and test points.  
# -----------------------------------------------------------------------------

print("--------------------- Constructing design -------------------------")

set.seed(seed_design)

# Generate design. 
design_info <- list(design_method=design_method, N_design=N_design, seed=seed_design)
design_info$input <- get_batch_design(design_info$design_method, N_design, 
                                      prior_params=inv_prob$par_prior)
design_info$fwd <- llik_exact$run_fwd_model(design_info$input)
design_info$llik <- llik_exact$assemble_llik(design_info$input)
design_info$lprior <- calc_lprior_theta(design_info$input, inv_prob$par_prior)
design_info$bounds <- get_bounds(design_info$input)
save(design_info, file=file.path(output_dir, "design_info.RData"))

# Test points. 
u_grid <- get_batch_design(design_method_test, N_design_test, prior_params=inv_prob$par_prior)
test_info <- list(input=u_grid, design_method=design_method_test, N_design=N_design_test)
test_info$fwd <- llik_exact$run_fwd_model(test_info$input) 
test_info$llik <- llik_exact$assemble_llik(test_info$input)
test_info$lprior <- calc_lprior_theta(test_info$input, inv_prob$par_prior)
save(test_info, file=file.path(output_dir, "test_info.RData"))

# -----------------------------------------------------------------------------
# Exact MCMC. 
# -----------------------------------------------------------------------------

# Initial condition and covariance proposal for MCMC. 
mcmc_par_init <- design_info$input[which.max(design_info$llik + design_info$lprior),]
cov_prop_init <- cov(design_info$input)
print(paste0("MCMC initial value:", mcmc_par_init))
print("Initial proposal covariance:")
print(cov_prop_init)

# MCMC sampling using exact likelihood. 
mcmc_exact_list <- mcmc_gp_noisy(inv_prob$llik_obj, inv_prob$par_prior, N_itr=N_mcmc, 
                                 mode="MCMH", par_init=mcmc_par_init, cov_prop=cov_prop_init)
samp_dt <- format_mcmc_output(mcmc_exact_list$samp, test_label="exact")
mcmc_list <- list(exact=mcmc_exact_list[setdiff(names(mcmc_exact_list), "samp")])


# -----------------------------------------------------------------------------
# Fitting Emulators.   
# -----------------------------------------------------------------------------

print("--------------------- Fitting Emulators -------------------------")

llik_em_list <- list()

# Log-likelihood emulator, constant mean. 
print("-> Fitting em_llik_const")
em_llik_gp <- gpWrapperHet(design_info$input, matrix(design_info$llik, ncol=1), 
                           normalize_output=TRUE, scale_input=TRUE)
em_llik_gp$fit("Gaussian", "constant", estimate_nugget=FALSE)
llik_em_list[["em_llik_const"]] <- llikEmulatorGP("em_llik_const", em_llik_gp, default_conditional=default_conditional, 
                                                  default_normalize=default_normalize, lik_par=inv_prob$sig2_eps, 
                                                  use_fixed_lik_par=TRUE)


# Log-likelihood emulator, quadratic mean.
print("-> Fitting em_llik_quad")
em_llik_gp2 <- gpWrapperKerGP(design_info$input, matrix(design_info$llik, ncol=1), 
                              normalize_output=TRUE, scale_input=TRUE)
em_llik_gp2$fit("Gaussian", "quadratic", estimate_nugget=FALSE, 
                optimFun="nloptr::nloptr", trace=TRUE, multistart=10)
llik_em_list[["em_llik_quad"]] <- llikEmulatorGP("em_llik_quad", em_llik_gp2, default_conditional=default_conditional, 
                                                 default_normalize=default_normalize, lik_par=inv_prob$sig2_eps, 
                                                 use_fixed_lik_par=TRUE)

# Log-likelihood emulator, constant mean, Gaussian plus quadratic kernel. 
print("-> Fitting em_llik_const_GaussQuad")
em_llik_gp3 <- gpWrapperKerGP(design_info$input, matrix(design_info$llik, ncol=1), 
                              normalize_output=TRUE, scale_input=TRUE)
em_llik_gp3$fit("Gaussian_plus_Quadratic", "constant", estimate_nugget=FALSE, 
                optimFun="nloptr::nloptr", trace=TRUE, multistart=10)
llik_em_list[["em_llik_const_GaussQuad"]] <- llikEmulatorGP("em_llik_const_GaussQuad", em_llik_gp3, 
                                                            default_conditional=default_conditional, 
                                                            default_normalize=default_normalize, 
                                                            lik_par=inv_prob$sig2_eps, use_fixed_lik_par=TRUE)

# Save log-likelihood emulators. 
save(llik_em_list, file=file.path(output_dir, "llik_em_list.RData"))

# -----------------------------------------------------------------------------
# Log-Likelihood predictions at test points.   
# -----------------------------------------------------------------------------

# Emulator predictions. 
print("--------------------- Emulator Predictions at Test Points -------------------------")
emulator_pred_list <- list()
for(em_name in names(llik_em_list)) {
  print(paste0("Predicting with emulator: ", em_name))
  emulator_pred_list[[em_name]] <- llik_em_list[[em_name]]$predict_emulator(test_info$input)
}
save(emulator_pred_list, file=file.path(output_dir, "emulator_pred_list.RData"))


# Log-likelihood predictions. 
print("--------------------- Log-Likelihood Predictions at Test Points -------------------------")
for(em_name in names(llik_em_list)) {
  plt <- llik_em_list[[em_name]]$plot_pred_validation(test_info$input, true_llik=test_info$llik, 
                                                      emulator_pred_list=emulator_pred_list[[em_name]], 
                                                      plot_type="llik", conditional=default_conditional,
                                                      normalize=default_normalize, include_interval=TRUE,
                                                      interval_method="pm_std_dev", N_std_dev=1.64, 
                                                      plot_title=paste0("llik predictions: ", em_name))
  ggsave(file.path(output_dir, paste0("llik_pred_", em_name, ".png")), plt)
}

# -----------------------------------------------------------------------------
# Run MCMC Samplers.  
# -----------------------------------------------------------------------------

print("--------------------- GP-Accelerated MCMC -------------------------")

if(run_approx_mcmc) {
  for(em_name in names(llik_em_list)) {
    print(paste0("----- ", em_name, " -----"))
    mcmc_info_list <- run_approx_mcmc_comparison(inv_prob_list=inv_prob, 
                                                 llik_em_obj=llik_em_list[[em_name]], 
                                                 mcmc_tags=mcmc_tags, save_dir=output_dir, 
                                                 samp_dt=samp_dt, mcmc_list=mcmc_list, 
                                                 test_label_suffix=em_name, N_itr=N_mcmc, 
                                                 par_init=mcmc_par_init, cov_prop=cov_prop_init,
                                                 overwrite=TRUE)
    
    samp_dt <- mcmc_info_list$samp
    mcmc_list <- mcmc_info_list$mcmc_list
  }
}

print("-------------------------- End of Script ------------------------------")




