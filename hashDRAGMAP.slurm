#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

echo $HOSTNAME
export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

ref=$2
refdir=$(echo ${ref%/*})

module load conda
conda activate $1

if [[ ! -f ${refdir}/hash_table.cmp ]]
then
    echo "Creating hash tables for DRAGMAP..."
    dragen-os \
        --build-hash-table true \
        --ht-reference $2 \
        --output-directory $refdir \
        --ht-num-threads 4
fi

(exit) && echo success
