#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 60
#SBATCH -t 20:00:00
#SBATCH -J tempExp

srun --output slurm/job%j-%2t.out --ntasks=60 --multi-prog makeCoaddTempExp_slurm.conf
