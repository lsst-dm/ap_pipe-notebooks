#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 2
#SBATCH -n 29
#SBATCH -t 20:00:00
#SBATCH -J assemble

srun --output slurm/job%j-%2t.out --ntasks=29 --multi-prog assemble_comparewarp_slurm.conf
