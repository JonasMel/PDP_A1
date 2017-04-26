#!/bin/bash -l
#SBATCH -A g2017012
#SBATCH -t 3:00
#SBATCH -n 1

module load gcc openmpi
mpirun -np 1 ./assign1 1 100 
