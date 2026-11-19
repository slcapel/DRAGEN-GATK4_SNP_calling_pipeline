#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`
ec=$7
ref=$3
refdir=$(echo ${ref%/*})
sm=$(cat $6 | grep -w ${sample} | awk '{print $2}')
lib=$(cat $6 | grep -w ${sample} | awk '{print $3}')
pu=$(cat $6 | grep -w ${sample} | awk '{print $4}')

echo $HOSTNAME
echo "${sample}"

[[ -d $2/${sm} ]] || mkdir -p $2/${sm}

module load conda
conda activate ${ec}
module load samtools
export TMPDIR=$4

if [[ ! -f $2/${sm}/${sample}.sorted.bam && ! -f $2/${sm}/${sample}.sorted.rg.bam && ! -f $2/${sm}/${sample}.sorted.rg.md.bam ]]
   then
echo "Aligning reads..."
dragen-os -v \
    -r $refdir \
    -1 $5/${sample}/${sample}_preproc_R1.fastq.gz \
    -2 $5/${sample}/${sample}_preproc_R2.fastq.gz \
    --num-threads 8 | samtools view -bh -@ 2 -m 3G - | samtools sort -@ 2 -O bam -o $2/${sm}/${sample}.sorted.bam
fi

if [[ ! -f $2/${sm}/${sample}.sorted.rg.bam && ! -f $2/${sm}/${sample}.sorted.rg.md.bam ]]
then
    echo "Replacing read groups..."
    picard AddOrReplaceReadGroups -I $2/${sm}/${sample}.sorted.bam -O $2/${sm}/${sample}.sorted.rg.bam -PL Illumina -LB ${lib} -SM ${sm} -PU ${pu} &&
	rm $2/${sm}/${sample}.sorted.bam
fi

if [[ ! -f $2/${sm}/${sample}.sorted.rg.md.bam ]]
   then
       module load gatk
       echo "Marking duplicates..."
       gatk MarkDuplicates \
	    -I $2/${sm}/${sample}.sorted.rg.bam \
	    -O $2/${sm}/${sample}.sorted.rg.md.bam \
	    -M $2/${sm}/${sample}.sorted.rg.md.metric.txt &&
	   rm $2/${sm}/${sample}.sorted.rg.bam
fi

if [[ $8 == "T" ]]
then
    rm $5/${sample}/${sample}_preproc_R*.fastq.gz
fi

(exit) && echo success
