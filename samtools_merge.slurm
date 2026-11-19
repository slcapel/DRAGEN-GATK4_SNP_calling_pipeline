#!/bin/bash

#SBATCH --nodes=1

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`
echo "${sample}"

module load samtools

if [[ ! -f $2/${sample}/${sample}.sorted.rg.md.bam ]]
then
    echo "Merging sample reads across lanes..."
    samtools merge -f -o $2/${sample}/${sample}.sorted.rg.md.bam $2/${sample}/${sample}*$3*.sorted.rg.md.bam
    rm $2/${sample}/${sample}*$3*.bam
fi

(exit) && echo success
