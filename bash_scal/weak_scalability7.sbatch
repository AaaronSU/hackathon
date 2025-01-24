#!/bin/bash
#SBATCH --job-name=weak_scalability
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --partition=c7g


THREADS=$SLURM_CPUS_PER_TASK
NUM_RUNS=1000000  

export OMP_NUM_THREADS=$THREADS

TOTAL_SIM=$(($1 * $THREADS))

./kk $TOTAL_SIM $NUM_RUNS
