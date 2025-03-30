#!/bin/bash -l
#$ -N run_emulator_seq_design
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.

# This script is intended to run the R script `emulator_seq_design.r` remotely   
# on the cluster. It is set up to be called from `run_emulator_seq_design_reps.r`.
# See the notes in these two R files for explanation.
#
# Arguments:
#   1: experiment tag.
#   2: round number
#   3: MCMC tag (from previous round)
#   4: MCMC ID (from previous round)
#   5: Number of points to acquire in the batch
#   6: Acquisition ID

# NOTE: 
# The R library path is not set explicitly in this file. Therefore, if using an 
# R project, one should ensure that renv::load() is called within the R script 
# being run so that the correct library path is set. 

# Read commandline arguments.
EXPERIMENT_TAG=$1
ROUND=$2
MCMC_TAG=$3
MCMC_ID=$4
N_BATCH=$5
ACQ_ID=$6

# Print arguments.
echo "Commandline arguments:"
echo "Experiment tag: ${EXPERIMENT_TAG}"
echo "Round: ${ROUND}"
echo "MCMC Tag: ${MCMC_TAG}"
echo "MCMC ID: ${MCMC_ID}"
echo "Batch Size: ${N_BATCH}"
echo "Acquisition ID: ${ACQ_ID}"

# For logging: log file will be moved from its default location.
LOG_FILENAME="${JOB_NAME}.o${JOB_ID}"
echo "Log file: ${LOG_FILENAME}"

# Main output directory for this run.
OUT_DIR="/projectnb/dietzelab/arober/gp-calibration/output/gp_inv_prob/${EXPERIMENT_TAG}/round${ROUND}/design/acq_${ACQ_ID}/${MCMC_TAG}/${MCMC_ID}"
echo "Creating output directory: ${OUT_DIR}"
mkdir -p ${OUT_DIR};

# Path to R script to execute.
R_SCRIPT_PATH="/projectnb/dietzelab/arober/gp-calibration/scripts/gp_post_approx_paper/sim_study/scripts/emulator_seq_design.r"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Execute R script.
Rscript ${R_SCRIPT_PATH} --experiment_tag=${EXPERIMENT_TAG} \
    --round=${ROUND} \
    --mcmc_tag=${MCMC_TAG} \
    --mcmc_id=${MCMC_ID} \
    --n_batch=${N_BATCH} \
    --acq_id=${ACQ_ID}

# Move output file to the output directory (if it exists). 
if [ -d ${OUT_DIR} ]; then
mv "${LOG_FILENAME}" "${OUT_DIR}/${LOG_FILENAME}"
fi






