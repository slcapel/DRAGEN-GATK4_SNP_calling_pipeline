#!/bin/bash

#SBATCH --job-name=merge
#SBATCH --nodes=1

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`

start=`date +%s`
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${sample}"

module load samtools
echo "Merging sample reads across lanes..."
samtools merge -f -o $2/${sample}/${sample}.sorted.rg.md.bam $2/${sample}/${sample}*$3*.sorted.rg.md.bam
rm $2/${sample}/*$3*.bam

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime seconds"

