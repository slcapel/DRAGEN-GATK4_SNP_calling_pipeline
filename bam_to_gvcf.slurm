#!/bin/bash

#SBATCH --job-name=GATK4hapcall
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1/indiv_scaff.tsv`
ref=$2
refdir=$(echo ${ref%/*})
samp=$(echo ${sample} | awk '{print $1}')
scaff=$(echo ${sample} |  awk '{print $2}')
lib=$(echo ${sample} |  awk '{print $3}')
[[ -d $1/${samp}_${scaff} ]] || mkdir -p $1/${samp}_${scaff}

start=`date +%s`
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "${sample}"
echo "${lib}"

echo "Calibrating DRAGEN STR model for haplotype calling..."
$3/gatk CalibrateDragstrModel \
     -R $2 \
     -I $4/${lib}/${samp}/*rg.md.bam \
     -L $1/${scaff}.bed \
     -str $refdir/str_table.tsv \
     -O $5/dragstr_model/${samp}_${scaff}.dragstr_model.txt

if [ $7 ]
then
    echo "Calling individual haplotypes over specified intervals..."
else
    echo "Calling individual haplotypes..."
fi
$3/gatk HaplotypeCaller \
     -R $2 \
     -I $4/${lib}/${samp}/*rg.md.bam \
     -O $5/indiv_scaff_gvcfs/${samp}_${scaff}.g.vcf.gz \
     -L $1/${scaff}.bed \
     -ERC GVCF \
     --dragen-mode true \
     --dragstr-params-path $5/dragstr_model/${samp}_${scaff}.dragstr_model.txt \
     --use-jdk-deflater \
     --use-jdk-inflater \
     --tmp-dir $1/${samp}_${scaff}

echo "Compressing gvcf..."
$3/gatk ReblockGVCF \
     -R $2 \
     -V $5/indiv_scaff_gvcfs/${samp}_${scaff}.g.vcf.gz \
     -O $5/indiv_scaff_gvcfs/${samp}_${scaff}.rb.g.vcf.gz \
     --use-jdk-deflater \
     --use-jdk-inflater &&

if [ $6 == "T" ]
then
    rm $4/${lib}/${samp}/*rg.md.bam
fi
rm $5/indiv_scaff_gvcfs/${samp}_${scaff}.g.vcf.gz*

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
