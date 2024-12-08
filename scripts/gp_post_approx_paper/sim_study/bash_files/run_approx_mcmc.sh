#!/bin/bash -l
#$ -N run_rmarkdown
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.
#$ -m beas                      # Email when job begins/ends/aborted/suspended

# This script is intended to run the R script `run_approx_mcmc.r` remotely on  
# the cluster. It can be executed directly to target a single 
# design/log-likelihood emulator, or can be called from `run_reps_approx_mcmc.sh`,
# which runs this script for different log-likelihood emulators across many 
# nodes in parallel. The following commandline arguments to this script should 
# be provided (in order), which are in turn passed as commandline arguments to 
# run_approx_mcmc.r`.
#
# Example executions of this script:
#    bash run_approx_mcmc.sh "vsem" "test_run" "init_emulator/LHS_250" "1018157756"
#    qsub run_approx_mcmc.sh "vsem" "test_run" "init_emulator/LHS_250" "1018157756"
#
# Arguments:
#   1: experiment tag.
#   2: run id, used to define output directory. 
#   3: llikEmulator base directory.
#   4: llikEmulator sub-directory.

# NOTE: 
# The R library path is not set explicitly in this file. Therefore, if using an 
# R project, one should ensure that renv::load() is called within the R script 
# being run so that the correct library path is set. 

# Read commandline arguments.
EXPERIMENT_TAG=$1
RUN_ID=$2
EM_DIR=$3
EM_ID=$4

# Print arguments.
echo "Commandline arguments:"
echo "Experiment tag: ${EXPERIMENT_TAG}"
echo "Run ID: ${RUN_ID}"
echo "llikEmulator base directory: ${EM_DIR}"
echo "llikEmulator sub directory: ${EM_ID}"

# For logging: log file will be moved from its default location.
LOG_FILENAME="${JOB_NAME}.o${JOB_ID}"
echo "Log file: ${LOG_FILENAME}"

# Main output directory for this run.
OUT_DIR="/projectnb/dietzelab/arober/gp-calibration/output/gp_inv_prob/${EXPERIMENT_TAG}/${RUN_ID}/${EM_DIR}/${EM_ID}"
echo "Creating output directory: ${OUT_DIR}"
mkdir -p ${OUT_DIR};

# Path to R script to execute.
R_SCRIPT_PATH="/projectnb/dietzelab/arober/gp-calibration/scripts/gp_post_approx_paper/sim_study/scripts/run_approx_mcmc.r"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Execute R script.
Rscript ${R_SCRIPT_PATH} --experiment_tag=${EXPERIMENT_TAG} \
    --run_id=${RUN_ID} \
    --em_dir=${EM_DIR} \
    --em_id=${EM_ID}

# Move output file to the output directory (if it exists). 
if [ -d ${OUT_DIR} ]; then
mv "${LOG_FILENAME}" "${OUT_DIR}/${LOG_FILENAME}"
fi






