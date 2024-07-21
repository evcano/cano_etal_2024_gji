#!/bin/bash
#SBATCH --ntasks=8
#SBATCH --time=20:00
#SBATCH --mem=10000
#SBATCH --nodes=1

mpirun -np 8 ./bin/xsmooth_sem 20000.0 20000.0 vs ./DATA ../target_velocity_model/complete_mesh/smoothed/ false
