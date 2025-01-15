#!/bin/bash

#SBATCH --job-name=GATK4genotype
#SBATCH --nodes=1

scaff=`sed "${SLURM_ARRAY_TASK_ID}q;d" $1/scaffolds.txt`

start=`date +%s`
echo $HOSTNAME
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo ${scaff}

module load bcftools

ls -d $2/indiv_scaff_gvcfs/*${scaff}.rb.g.vcf.gz | sed -E "s/(.+\/)(.+)\_${scaff}\.rb.g.vcf.gz/\2\t\1\2\_${scaff}.rb.g.vcf.gz/" > $1/samples_${scaff}.gvcf.txt

# GenomicsDBImport
echo "Importing single-sample gvcfs into database..."
$3/gatk GenomicsDBImport \
     --genomicsdb-workspace-path $2/genomicDB/${scaff}.gvcf.db \
     --intervals $1/${scaff}.bed \
     --sample-name-map $1/samples_${scaff}.gvcf.txt \
     --use-jdk-deflater \
     --use-jdk-inflater \
     --tmp-dir $1 \
     --batch-size 26

# GenotypeGVCFs
echo "Joint genotyping all individuals..."
$3/gatk GenotypeGVCFs \
     --reference $4 \
     --variant gendb://$2/genomicDB/${scaff}.gvcf.db \
     --output $2/scaff_vcfs/${scaff}.gatk.vcf.gz \
     --tmp-dir $1
zcat $2/scaff_vcfs/${scaff}.gatk.vcf.gz | grep -m 1 "#CHR" | cut -f 10- | tr "\t" "\n" | sed -E "s/\_$scaff//" > $1/${scaff}.vcf.rehead

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
