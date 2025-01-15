#!/bin/bash

#SBATCH --job-name=drgmp_hash
#SBATCH --ntasks-per-node=32

start=`date +%s`
echo $HOSTNAME

dc=$(cat $1 | awk '{print $1}')
ec=$(cat $1 | awk '{print $2}')

source ${dc}/conda.sh
conda activate ${ec}
ref=$2
refdir=$(echo ${ref%/*})

echo "Creating hash tables for DRAGMAP..."
dragen-os \
    --build-hash-table true \
    --ht-reference $2 \
    --output-directory $refdir \
    --ht-num-threads 32

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
