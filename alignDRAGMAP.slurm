#!/bin/bash

#SBATCH --job-name=DRAGMAP
#SBATCH --ntasks-per-node=32

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`
ind=$(echo ${sample} | sed -E "s/$8.+//")
ref=$3
refdir=$(echo ${ref%/*})
sm=$(cat $6 | grep -w ${sample} | awk '{print $2}')
lib=$(cat $6 | grep -w ${sample} | awk '{print $3}')
pu=$(cat $6 | grep -w ${sample} | awk '{print $4}')
dc=$(cat $4/conda.txt | awk '{print $1}')
ec=$(cat $4/conda.txt | awk '{print $2}')

start=`date +%s`
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${sample}"

[[ -d $2/${ind} ]] || mkdir -p $2/${ind}

module load samtools
aklog
source ${dc}/conda.sh
conda activate ${ec}
export TMPDIR=$4

if [[ ! -f $2/${ind}/${sample}.sorted.bam ]]
   then
echo "Aligning reads..."
dragen-os -v \
    -r $refdir \
    -1 $5/${sample}/${sample}_preproc_R1.fastq.gz \
    -2 $5/${sample}/${sample}_preproc_R2.fastq.gz \
    --num-threads 32 | \
    samtools view -bh - | \
    samtools sort -@ 2 -O bam -o $2/${ind}/${sample}.sorted.bam
fi

echo "Replacing read groups..."
picard AddOrReplaceReadGroups -I $2/${ind}/${sample}.sorted.bam -O $2/${ind}/${sample}.sorted.rg.bam -PL Illumina -LB ${lib} -SM ${sm} -PU ${pu} &&
rm $2/${ind}/${sample}.sorted.bam &&


echo "Marking duplicates..."        
$7/gatk MarkDuplicates \
     -I $2/${ind}/${sample}.sorted.rg.bam \
     -O $2/${ind}/${sample}.sorted.rg.md.bam \
     -M $2/${ind}/${sample}.sorted.rg.md.metric.txt &&
rm $2/${ind}/${sample}.sorted.rg.bam

if [ $9 == "T" ]
then
    rm $5/${sample}/${sample}_preproc_R*.fastq.gz
fi

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
