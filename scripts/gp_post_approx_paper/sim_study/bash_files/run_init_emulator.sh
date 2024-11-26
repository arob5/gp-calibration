#!/bin/bash -l
#$ -N run_rmarkdown
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.
#$ -m beas                      # Email when job begins/ends/aborted/suspended

# This script is intended to run the R script `init_emulator.r` remotely on the 
# cluster. This script is intended to be executed from `run_reps_init_emulator.sh`
# to run replications of `init_emulator.r` across many nodes in parallel. The 
# following commandline arguments to this script should be provided (in order), 
# which are in turn passed as commandline arguments to init_emulator.r`.
#
# Example call to this script:
#    bash run_init_emulator.sh "test_run" 100 "LHS" 10
#
# Arguments:
#   1: experiment tag.
#   2: number of design points.
#   3: design method. 
#   4: number of replicate designs created by init_emulator.r`.

# NOTE: 
# The R library path is not set explicitly in this file. Therefore, if using an 
# R project, one should ensure that renv::load() is called within the R script 
# being run so that the correct library path is set. 

# Read commandline arguments.
EXPERIMENT_TAG=$1
N_DESIGN=$2
DESIGN_METHOD=$3
N_REP=$4

# Print arguments.
echo "Commandline arguments:"
echo "Experiment tag: ${EXPERIMENT_TAG}"
echo "Number of design points: ${N_DESIGN}"
echo "Design method: ${DESIGN_METHOD}"
echo "Number of replicate designs: ${N_REP}"

# Define a design ID used to create the directory in which all replicate 
# designs will be saved (each in their own subdirectory).
DESIGN_ID="${DESIGN_METHOD}_${N_DESIGN}"
echo "Created design ID: ${DESIGN_ID}"

# For logging: log file will be moved from its default location.
LOG_FILENAME="${JOB_NAME}.o${JOB_ID}"
echo "Log file: ${LOG_FILENAME}"

# Main output directory for this run.
OUT_DIR="/projectnb/dietzelab/arober/gp-calibration/output/gp_inv_prob/${EXPERIMENT_TAG}/init_emulator/${DESIGN_ID}"
echo "Creating output directory: ${OUT_DIR}"
mkdir -p ${OUT_DIR};

# Path to R script to execute.
R_SCRIPT_PATH="/projectnb/dietzelab/arober/gp-calibration/scripts/gp_post_approx_paper/sim_study/scripts/init_emulator.r"

# Load modules, including specification of the R version. 
module load gcc/8.3.0
module load R/4.3.1

# Required for the `Rscript` command to be found. 
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

# Execute R script.
Rscript ${R_SCRIPT_PATH} --run_id=${DESIGN_ID} \
    --n_design=${N_DESIGN} \
    --design_method=${DESIGN_METHOD} \
    --n_rep=${N_REP}

# Move output file to the output directory (if it exists). 
if [ -d ${OUT_DIR} ]; then
mv "${LOG_FILENAME}" "${OUT_DIR}/${LOG_FILENAME}"
fi






