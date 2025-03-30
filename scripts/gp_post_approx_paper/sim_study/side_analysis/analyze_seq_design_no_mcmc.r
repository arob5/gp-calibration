#
# analyze_seq_design_no_mcmc.r
#
# Formats and produces analysis of the outputs saved by 
# `run_seq_design_no_mcmc_reps.r`
#
# Andrew Roberts
# 

# ------------------------------------------------------------------------------
# Plot Results 
# ------------------------------------------------------------------------------

comp_quant <- acq_results$tracking_list$computed_quantities

dt_results <- data.table(itr = integer(),
                         metric = numeric(),
                         dataset = character(),
                         value = numeric())

for(i in seq_along(comp_quant)) {
  
  itr <- as.integer(strsplit(names(comp_quant)[i], "_", fixed=TRUE)[[1]][2])
  for(j in seq_along(comp_quant[[i]])) {
    comp_quant[[i]][[j]][, dataset := names(comp_quant[[i]])[j]]
  }
  
  dt <- rbindlist(comp_quant[[i]], use.names=TRUE)
  dt[, itr := itr]
  setnames(dt, c("func", "mean"), c("metric", "value"))
  
  dt_results <- rbindlist(list(dt_results, dt), use.names=TRUE)
}


plt <- ggplot(dt_results[metric=="crps" & dataset=="gp_eval_post"]) + geom_line(aes(x=itr, y=value))


plot(plt)


acq_val <- acq_results$tracking_list$acq_val
plot(seq_along(acq_val), acq_val, type="l")
