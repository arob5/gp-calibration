#!/bin/bash -l
#$ -N run_seq_design_no_mcmc
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.

# This script is intended to run the R script `seq_design_no_mcmc.r` remotely on  
# the cluster. It is typically called from `run_seq_design_no_mcmc_reps.r` in 
# order to dispatch one remote job per emulator ID.
# The following commandline arguments to this script should 
# be provided (in order), which are in turn passed as commandline arguments to 
# `seq_design_no_mcmc.r`.
#
# Arguments:
#   1: experiment tag.
#   2: design tag.
#   3: emulator tag.
#   4: emulator ID.
#   5: number of candidate points for discrete optimization.
#   6: number of points to acquire in the optimization.

# NOTE: 
# The R library path is not set explicitly in this file. Therefore, if using an 
# R project, one should ensure that renv::load() is called within the R script 
# being run so that the correct library path is set. 

# Read commandline arguments.
EXPERIMENT_TAG=$1
DESIGN_TAG=$2
EM_TAG=$3
EM_ID=$4
N_CANDIDATES=$5
N_BATCH=$6

# Print arguments.
echo "Commandline arguments:"
echo "Experiment tag: ${EXPERIMENT_TAG}"
echo "Design tag: ${DESIGN_TAG}"
echo "Emulator tag: ${EM_TAG}"
echo "Emulator ID: ${EM_ID}"
echo "Number of candidates for optimization: ${N_CANDIDATES}"
echo "Number of points to acquire in optimization: ${N_BATCH}"

# For logging: log file will be moved from its default location.
LOG_FILENAME="${JOB_NAME}.o${JOB_ID}"
echo "Log file: ${LOG_FILENAME}"

# Main output directory for this run.
OUT_DIR="/projectnb/dietzelab/arober/gp-calibration/output/gp_inv_prob/${EXPERIMENT_TAG}/side_analysis/seq_design_no_mcmc/${EM_TAG}/${EM_ID}/${DESIGN_TAG}"
echo "Creating output directory: ${OUT_DIR}"
mkdir -p ${OUT_DIR};

# Path to R script to execute.
R_SCRIPT_PATH="/projectnb/dietzelab/arober/gp-calibration/scripts/gp_post_approx_paper/sim_study/side_analysis/seq_design_no_mcmc.r"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Execute R script.
Rscript ${R_SCRIPT_PATH} --experiment_tag=${EXPERIMENT_TAG} \
    --design_tag=${DESIGN_TAG} \
    --em_tag=${EM_TAG} \
    --em_id=${EM_ID} \
    --n_candidates=${N_CANDIDATES} \
    --n_batch=${N_BATCH}

# Move output file to the output directory (if it exists). 
if [ -d ${OUT_DIR} ]; then
mv "${LOG_FILENAME}" "${OUT_DIR}/${LOG_FILENAME}"
fi
