#!/bin/bash

#SBATCH --nodes=1

ref=$1
refpre=$(echo ${ref%.*})

module load bedtools2

if [[ ! -f ${refpre}_50kb_win.bed ]]
   then
echo "breaking genome into 50kb windows..."
cat $1.fai | awk '{print $1"\t"0"\t"$2-1}' > $refpre.bed
bedtools makewindows -b $refpre.bed -w 50000 -s 45000 > ${refpre}_50kb_win.bed
fi

(exit) && echo success
