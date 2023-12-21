#!/bin/bash -l

# Hard time limit (default 12 hours)
#$ -l h_rt=20:00:00

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

# Job name
#$ -N mcmc_pecan

# Request 4 cores
#$ -pe omp 4

# When a parallel environment is requested the environment variable NSLOTS is 
# set to the number of cores requested. This variable can be used within a 
# program to setup an appropriate number of threads or processors to use.
# For example, some programs rely on the environment variable OMP_NUM_THREADS for parallelization:
OMP_NUM_THREADS=$NSLOTS

Rscript /projectnb2/dietzelab/arober/test_code/run.stan.test.R

