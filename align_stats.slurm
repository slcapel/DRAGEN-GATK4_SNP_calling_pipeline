#!/bin/bash

#SBATCH --nodes=1

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1`
echo "${sample}"
[[ -d $2/${sample} ]] || mkdir -p $2/${sample}
ref=$3
refpre=$(echo ${ref%.*})

module load samtools
module load bedtools2

if [[ ! -f $4/${sample}/${sample}.sorted.rg.md.bam.bai ]]
then
    echo "indexing bam..."
    samtools index $4/${sample}/${sample}*.bam
fi

if [[ ! -f $2/${sample}/${sample}.dragmap_coverage.tsv ]]
then
    echo "Calculating samtools stats:"
    echo "  flagstat..."
    samtools flagstat $4/${sample}/${sample}*.bam > $2/${sample}/${sample}.dragmap_flagstat.txt
    echo "  stats..."
    samtools stats -c 1,10000,100 $4/${sample}/${sample}*.bam > $2/${sample}/${sample}.dragmap_stats.txt
    echo "  idxstats..."
    samtools idxstats $4/${sample}/${sample}*.bam > $2/${sample}/${sample}.dragmap_IDstats.tsv
    echo "  depth..."
    samtools depth $4/${sample}/${sample}*.bam | awk '{sum+=$3; cnt++}END{print "no_positions\tmean_depth"; print sum"\t"sum/cnt}' > $2/${sample}/${sample}.dragmap_depth.txt
    echo "  coverage..."
    samtools coverage $4/${sample}/${sample}*.bam > $2/${sample}/${sample}.dragmap_coverage.tsv
fi

if [[ ! -f $2/${sample}/${sample}.dragmap.gatk.stat.txt ]]
then
    module load gatk
    echo "Calculating GATK metrics..."
    echo "checking mapping quality"
    gatk CollectAlignmentSummaryMetrics \
	 -R $3 \
	 -I $4/${sample}/${sample}*.bam \
	 -O $2/${sample}/${sample}.dragmap.gatk.stat.txt
fi

if [[ ! -f $2/${sample}/${sample}.dragmap_50kbwin_meancov.tsv ]]
then
    echo "Calculating coverage for 50Kb windows..."
    echo -e "scaff\twin_start\twin_end\tmean_depth" > $2/${sample}/${sample}.dragmap_50kbwin_meancov.tsv
    bedtools coverage -mean -sorted -a ${refpre}_50kb_win.bed -b $4/${sample}/${sample}*.bam >> $2/${sample}/${sample}.dragmap_50kbwin_meancov.tsv
fi

(exit) && echo success
