#!/bin/bash -l
#$ -N multidim_toy_example
#$ -P dietzelab                 # Specify project. 
#$ -l buyin                     # Request buyin node. 
#$ -j y                         # Merge the error and output streams into a single file.
#$ -m beas                      # Email when job begins/ends/aborted/suspended

module load gcc/8.3.0
module load R/4.3.1
export PATH="$PATH:/share/pkg.8/r/4.3.1/install/lib64/R"

Rscript ../run_multidim_toy_example.r


