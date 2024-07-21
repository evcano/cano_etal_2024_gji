#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --time=20:00
#SBATCH --mem=10000
#SBATCH --nodes=1

mpirun -np 8 ./bin/xsmooth_sem 14000.0 14000.0 gradient_mask ./DATA ./ false
