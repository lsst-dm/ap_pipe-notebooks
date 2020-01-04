#!/bin/bash -l

#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 60
#SBATCH -t 12:00:00
#SBATCH -J ap_pipe

# To submit a slurm job:
#       $ sbatch run_ap_pipe.sl

VISITS=`sqlite3 /project/mrawls/hits2015/registry.sqlite3 "select distinct visit from raw where filter = 'g';"`

for VISIT in ${VISITS};
do
    echo "Processing ${VISIT}"
    VISITSTR=$( printf '%06d' $VISIT )
    srun --output slurm/ap_pipe%j-${VISITSTR}-%2t.out --multi-prog run_ap_pipe.conf visit=${VISIT}
done
wait
