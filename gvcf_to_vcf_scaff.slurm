#!/bin/bash

#SBATCH --nodes=1

scaff=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1/scaffolds.txt`
echo ${scaff}

module load bcftools
module load gatk

ls -d $2/indiv_scaff_gvcfs/*${scaff}.rb.g.vcf.gz | sed -E "s/(.+\/)(.+)\_${scaff}\.rb.g.vcf.gz/\2\t\1\2\_${scaff}.rb.g.vcf.gz/" > $1/samples_${scaff}.gvcf.txt

# GenomicsDBImport
if [[ ! -f $2/genomicDB/${scaff}.gvcf.db ]]
then
   echo "Importing single-sample gvcfs into database..."
   gatk GenomicsDBImport \
	--genomicsdb-workspace-path $2/genomicDB/${scaff}.gvcf.db \
	--intervals $1/${scaff}.bed \
	--sample-name-map $1/samples_${scaff}.gvcf.txt \
	--use-jdk-deflater \
	--use-jdk-inflater \
	--tmp-dir $1 \
	--batch-size 26
fi

# GenotypeGVCFs
if [[ ! -f $2/scaff_vcfs/${scaff}.gatk.vcf.gz ]]
then
    echo "Joint genotyping all individuals..."
    gatk --java-options "-Xmx12g" GenotypeGVCFs \
	 --reference $3 \
	 --variant gendb://$2/genomicDB/${scaff}.gvcf.db \
	 --output $2/scaff_vcfs/${scaff}.gatk.vcf.gz \
	 --tmp-dir $1
fi

(exit) && echo success
