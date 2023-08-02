
# -------------------------------------------------------------------
# Create lpost emulator objects. 
# -------------------------------------------------------------------

# Emulator trained on all data. 
gp_fits <- fit_independent_GPs(X_train = design_info$inputs_scaled, 
                               Y_train = design_info$outputs_normalized, 
                               gp_lib = emulator_settings$gp_lib, 
                               gp_kernel = emulator_settings$kernel)$fits
emulator_info_list <- list(gp_fits = gp_fits, 
                           input_bounds = design_info$input_bounds, 
                           output_stats = design_info$output_stats, 
                           settings = emulator_settings, 
                           emulator_target = "SSR")
lpost_emulator <- get_lpost_emulator_obj(emulator_info_list, design_info, computer_model_data, sig_eps_estimates, theta_prior_params)

# Emulator trained on all but the last design point. 
N_design <- nrow(design_info$inputs_scaled)
design_info_nm1 <- design_info
design_info_nm1$inputs <- design_info_nm1$inputs[1:(N_design-1),,drop=FALSE]
design_info_nm1$inputs_scaled <- design_info_nm1$inputs_scaled[1:(N_design-1),,drop=FALSE]
design_info_nm1$outputs <- design_info_nm1$outputs[1:(N_design-1),,drop=FALSE]
design_info_nm1$outputs_normalized <- design_info_nm1$outputs_normalized[1:(N_design-1),,drop=FALSE]

# Use emulators fit on full data so that the hyperparameter estimates are the same, but use the reduced design set here. 
lpost_emulator_nm1 <- get_lpost_emulator_obj(emulator_info_list, design_info_nm1, computer_model_data, sig_eps_estimates, theta_prior_params)


# -------------------------------------------------------------------
# Test lpost emulator predictions.   
# -------------------------------------------------------------------

# Test data. 
inputs_test <- get_input_design(N_points = 2, theta_prior_params = theta_prior_params, design_method = "LHS", scale_inputs = TRUE, 
                                param_ranges = design_info$input_bounds)

# Baseline predictions. 
pred1 <- predict_lpost_GP_approx(theta_vals_scaled = inputs_test$inputs_scaled, 
                                 theta_vals_unscaled = inputs_test$inputs, 
                                 emulator_info_list = emulator_info_list, 
                                 sig2_eps = sig_eps_estimates, theta_prior_params = theta_prior_params, 
                                 N_obs = computer_model_data$n_obs, include_nugget = TRUE)

# Predictions using lpost emulator function. 
pred2 <- predict_lpost_emulator(inputs_new_scaled = inputs_test$inputs_scaled, lpost_emulator = lpost_emulator, inputs_new_unscaled = inputs_test$inputs)

print("Variances:")
print(paste0("Baseline: ", paste0(pred1$var, collapse = ", ")))
print(paste0("lpost: ", paste0(pred2$var, collapse = ", ")))


print("Means:")
print(paste0("Baseline: ", paste0(pred1$mean, collapse = ", ")))
print(paste0("lpost: ", paste0(pred2$mean, collapse = ", ")))

# -------------------------------------------------------------------
# Test lpost inverse kernel matrix update.    
# -------------------------------------------------------------------

# Update inverse kernel matrix using Nth data point. 
K_inv_update <- update_lpost_inverse_kernel_matrix(lpost_emulator_nm1, input_new_scaled = design_info$inputs_scaled[N_design,,drop=FALSE], 
                                                   include_nugget = TRUE)

print(K_inv_update)
print(lpost_emulator$K_inv)
print(paste0("Max absolute error: ", max(abs(K_inv_update - lpost_emulator$K_inv))))


# -------------------------------------------------------------------
# Test lpost emulator update.    
# -------------------------------------------------------------------

lpost_emulator_update <- update_lpost_emulator(lpost_emulator_nm1, inputs_new_scaled=design_info$inputs_scaled[N_design,,drop=FALSE], 
                                               outputs_lpost_new = lpost_emulator$outputs_lpost[N_design])

print(all.equal(lpost_emulator$inputs_lpost, lpost_emulator_update$inputs_lpost))
print(all.equal(lpost_emulator$outputs_lpost, lpost_emulator_update$outputs_lpost))
print(paste0("Max absolute error: ", max(abs(lpost_emulator_update$K_inv - lpost_emulator$K_inv))))

# -------------------------------------------------------------------
# Code speed tests.    
# -------------------------------------------------------------------

t0 <- proc.time()
test <- acquisition_EIVAR_lpost(theta_vals = round1_candidate_samp_scaled[1:10,,drop=FALSE], lpost_emulator = lpost_emulator, 
                                theta_grid_integrate = round1_integrate_samp_scaled)
time_elapsed <- (proc.time() - t0)[["elapsed"]]
print(paste0("Time: ", time_elapsed))

print(paste0("Extrapolated time for 10000 evaluations: ", time_elapsed * (10000/10)))


## My update. 
t0_lpost <- proc.time()
for(i in 1:100) {
  # Update variance by conditioning on theta evaluation value. Should not affect `lpost_emulator` outside of local function scope.
  lpost_emulator_temp <- update_lpost_emulator(lpost_emulator, inputs_new_scaled = round1_candidate_samp_scaled[i,,drop=FALSE], outputs_lpost_new = 0)
}
t1_lpost <- proc.time()
time_elapsed_lpost <- (t1_lpost - t0_lpost)[["elapsed"]]
print(time_elapsed_lpost)


## hetGP update. 
hetgp_obj <- emulator_info_list$gp_fits[[1]]

t0_hetGP <- proc.time()
for(i in 1:100) {
  # Update variance by conditioning on theta evaluation value. Should not affect `lpost_emulator` outside of local function scope.
  hetgp_obj_temp <- update(object = hetgp_obj, Xnew = round1_candidate_samp_scaled[i,,drop=FALSE], Znew = NULL, maxit = 0)
}
t1_hetGP <- proc.time()
time_elapsed_hetGP <- (t1_hetGP - t0_hetGP)[["elapsed"]]
print(time_elapsed_hetGP)


## My predict 
mu_new <- calc_lpost_mean(lpost_emulator, inputs_scaled = round1_candidate_samp_scaled[1:100,,drop=FALSE])

t0_lpost <- proc.time()
for(i in 1:100) {
  test <- predict_lpost_emulator(inputs_new_scaled = round1_candidate_samp_scaled[1:100,,drop=FALSE],  
                                 lpost_emulator = lpost_emulator, include_nugget = TRUE,
                                 prior_mean_vals_new = mu_new)
}
t1_lpost <- proc.time()
time_elapsed_lpost <- (t1_lpost - t0_lpost)[["elapsed"]]
print(time_elapsed_lpost)

## hetGP predict.  
t0_hetGP <- proc.time()
for(i in 1:100) {
  test <- predict(hetgp_obj, round1_candidate_samp_scaled[1:100,,drop=FALSE])
}
t1_hetGP <- proc.time()
time_elapsed_hetGP <- (t1_hetGP - t0_hetGP)[["elapsed"]]
print(time_elapsed_hetGP)













