#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 30:00:00
#SBATCH -J tempExp

srun --output slurm/job%j-%2t.out --ntasks=2 --multi-prog makeCoaddTempExp_slurm.conf
