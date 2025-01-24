#!/bin/bash

g++ cg_opti2.cxx -o kk -Ofast -mcpu=neoverse-v1 -funroll-loops -ftree-vectorize -flto -ffp-contract=fast -fopenmp -DNDEBUG -g0 -lm -fomit-frame-pointer


sbatch --exclusive --cpus-per-task=64 --output=weak_scalability7_64.log --error=weak_scalability7_64.err weak_scalability7.sbatch 10000     
sbatch --exclusive --cpus-per-task=32 --output=weak_scalability7_32.log --error=weak_scalability7_32.err weak_scalability7.sbatch 10000     

sbatch --exclusive --cpus-per-task=64 --output=strong_scalability7_64.log --error=strong_scalability7_64.err strong_scalability7.sbatch 100000      
sbatch --exclusive --cpus-per-task=32 --output=strong_scalability7_32.log --error=strong_scalability7_32.err strong_scalability7.sbatch 100000      

sbatch --exclusive --cpus-per-task=16 --output=weak_scalability7_16.log --error=weak_scalability7_16.err weak_scalability7.sbatch 10000     
sbatch --exclusive --cpus-per-task=8 --output=weak_scalability7_8.log --error=weak_scalability7_8.err weak_scalability7.sbatch 10000    
sbatch --exclusive --cpus-per-task=4 --output=weak_scalability7_4.log --error=weak_scalability7_4.err weak_scalability7.sbatch 10000    
sbatch --exclusive --cpus-per-task=2 --output=weak_scalability7_2.log --error=weak_scalability7_2.err weak_scalability7.sbatch 10000    
sbatch --exclusive --cpus-per-task=1 --output=weak_scalability7_1.log --error=weak_scalability7_1.err weak_scalability7.sbatch 10000    

sbatch --exclusive --cpus-per-task=16 --output=strong_scalability7_16.log --error=strong_scalability7_16.err strong_scalability7.sbatch 100000      
sbatch --exclusive --cpus-per-task=8 --output=strong_scalability7_8.log --error=strong_scalability7_8.err strong_scalability7.sbatch 100000     
sbatch --exclusive --cpus-per-task=4 --output=strong_scalability7_4.log --error=strong_scalability7_4.err strong_scalability7.sbatch 100000     
sbatch --exclusive --cpus-per-task=2 --output=strong_scalability7_2.log --error=strong_scalability7_2.err strong_scalability7.sbatch 100000     
sbatch --exclusive --cpus-per-task=1 --output=strong_scalability7_1.log --error=strong_scalability7_1.err strong_scalability7.sbatch 100000     

