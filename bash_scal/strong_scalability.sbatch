#!/bin/bash
#SBATCH --job-name=strong_scalability
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --partition=c8g


NUM_RUNS=1000000        
TOTAL_SIM=$1          
THREADS=$SLURM_CPUS_PER_TASK

export OMP_NUM_THREADS=$THREADS

./kk $TOTAL_SIM $NUM_RUNS
