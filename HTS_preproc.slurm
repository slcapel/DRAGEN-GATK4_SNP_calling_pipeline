#!/bin/bash

#SBATCH --job-name=htstream
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --cpus-per-task=1

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`
[[ -d $2/${sample} ]] || mkdir -p $2/${sample}

start=`date +%s`
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${sample}"

module load htstream

echo "Cleaning raw sequences..."
hts_Stats -L $2/${sample}/${sample}.json -N 'initial stats' \
	  -1 $3/${sample}_R1_001.fastq.gz \
	  -2 $3/${sample}_R2_001.fastq.gz | \
    hts_SeqScreener -A $2/${sample}/${sample}.json -N 'screen phix' | \
        hts_AdapterTrimmer -A $2/${sample}/${sample}.json -N 'trim adapters' | \
    hts_QWindowTrim -A $2/${sample}/${sample}.json -N 'quality trim the ends of reads' | \
    hts_NTrimmer -A $2/${sample}/${sample}.json -N 'remove any remaining N characters' | \
    hts_LengthFilter -A $2/${sample}/${sample}.json -N "remove reads < 100bp" -n -m 100 | \
    hts_Stats -A $2/${sample}/${sample}.json -N 'final stats' -F -f $2/${sample}/${sample}_preproc

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
