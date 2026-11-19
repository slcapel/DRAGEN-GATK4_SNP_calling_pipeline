#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1/indiv_scaff.tsv`
ref=$2
refdir=$(echo ${ref%/*})
samp=$(echo ${sample} | awk '{print $1}')
scaff=$(echo ${sample} |  awk '{print $2}')
lib=$(echo ${sample} |  awk '{print $3}')
[[ -d $1/${samp}_${scaff} ]] || mkdir -p $1/${samp}_${scaff}

module load gatk
echo "${sample}"
echo "${lib}"

if [[ ! -f $4/dragstr_model/${samp}_${scaff}.dragstr_model.txt ]]
then
    echo "Calibrating DRAGEN STR model for haplotype calling..."
    gatk CalibrateDragstrModel \
	 --parallel true \
	 --threads 4 \
	 -R $2 \
	 -I $3/${lib}/${samp}/*rg.md.bam \
	 -L $1/${scaff}.bed \
	 -str $refdir/str_table.tsv \
	 -O $4/dragstr_model/${samp}_${scaff}.dragstr_model.txt
fi

if [[ ! -f $4/indiv_scaff_gvcfs/${samp}_${scaff}.rb.g.vcf.gz ]]
then
    echo "Calling individual haplotypes..."
    gatk HaplotypeCaller \
	 --native-pair-hmm-threads 4 \
	 -R $2 \
	 -I $3/${lib}/${samp}/*rg.md.bam \
	 -O $4/indiv_scaff_gvcfs/${samp}_${scaff}.g.vcf.gz \
	 -L $1/${scaff}.bed \
	 -ERC GVCF \
	 --dragen-mode true \
	 --dragstr-params-path $4/dragstr_model/${samp}_${scaff}.dragstr_model.txt \
	 --use-jdk-deflater \
	 --use-jdk-inflater \
	 --tmp-dir $1/${samp}_${scaff}

    echo "Compressing gvcf..."
    gatk ReblockGVCF \
	 -R $2 \
	 -V $4/indiv_scaff_gvcfs/${samp}_${scaff}.g.vcf.gz \
	 -O $4/indiv_scaff_gvcfs/${samp}_${scaff}.rb.g.vcf.gz \
	 --use-jdk-deflater \
	 --use-jdk-inflater
fi && rm $4/indiv_scaff_gvcfs/${samp}_${scaff}.g.vcf.gz*

if [ $5 == "T" ]
then
    rm $3/${lib}/${samp}/*rg.md.bam
fi

(exit) && echo success
