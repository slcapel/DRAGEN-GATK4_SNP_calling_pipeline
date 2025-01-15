#!/bin/bash

#SBATCH --job-name=genome_win
#SBATCH --nodes=1

start=`date +%s`
echo $HOSTNAME

ref=$1
refpre=$(echo ${ref%.*})

module load bedtools2

echo "breaking genome into 50kb windows..."
cat $1.fai | awk '{print $1"\t"0"\t"$2-1}' > $refpre.bed
bedtools makewindows -b $refpre.bed -w 50000 -s 45000 > ${refpre}_50kb_win.bed

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime seconds"
