#
# run_seq_design_no_mcmc_reps.r
#
# The main starting point for running the sequential design (no MCMC) workflow.
# The file is set up so that the user specifies the experiment tag, emulator
# tag, and design tag. A separate job is then submitted for each emulator ID
# found within that emulator tag.
#
# Andrew Roberts
#

library(data.table)

# ------------------------------------------------------------------------------
# Settings.
# ------------------------------------------------------------------------------

# Determines which emulator models will be used for sequential design.
experiment_tag <- "vsem"
em_tag <- "llik_quad_mean"
design_tag <- "LHS_200"

# `n_candidates` is the number of points over which the discrete optimization
# will take place; `n_batch` is how many new design points will be acquired.
n_candidates <- 1000L
n_batch <- 200L

print(paste0("n_batch: ", n_batch))
print(paste0("n_candidates: ", n_candidates))

# Directories.
base_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration")
experiment_dir <- file.path(base_dir, "output", "gp_inv_prob", experiment_tag)
em_dir <- file.path(experiment_dir, "round1", "em", em_tag)
base_design_dir <- file.path(experiment_dir, "round1", "design")
inv_prob_dir <- file.path(experiment_dir, "inv_prob_setup")

# Path to bash file, which runs the script.
bash_path <- file.path(base_dir, "scripts", "gp_post_approx_paper",
                       "sim_study", "side_analysis", "run_seq_design_no_mcmc.sh")


# ------------------------------------------------------------------------------
# Start jobs.
# ------------------------------------------------------------------------------

# The ID map file is used to select only emulator IDs associated with the 
# specified design tag.
id_map <- fread(file.path(em_dir, "id_map.csv"))
design_tag_curr <- design_tag
id_map <- id_map[design_tag==design_tag_curr]

# Bash command for running bash script, will be appended to.
base_cmd <- paste("qsub", bash_path, experiment_tag, design_tag, em_tag)

# Dispatch one remote job per emulator ID.
print("Starting jobs with em_ids:")
for(em_id in id_map$em_id) {
  print(em_id)
  cmd <- paste(base_cmd, em_id, n_candidates, n_batch)
  system(cmd)
}
  
  









