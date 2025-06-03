#!/bin/bash

#SBATCH --job-name=gvcf2snp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

start=`date +%s`
echo $HOSTNAME


module load bcftools
nind=$0
minMQ=166.67
minDP=10
maxDP=$(cat $5/*/*depth.txt | grep -E [0-9] | cut -f 2 | sort -k1,1gr | head -n 1 | awk '{print $1*3}')
minmaf=$(awk -v nind=$nind 'BEGIN {printf "%.6f\n", 5 / (2 * nind) }')
minmafi=$(awk -v minmaf=$minmaf 'BEGIN {printf "%.6f\n", 1 - minmaf }')
sc1=$(cat $6 | awk '{print $1}')
sc2=$(cat $6 | awk '{print $2}')

echo "Concatenating vcfs..."
bcftools concat -o $1/gatk.all_indivs.vcf.gz -O z $1/scaff_vcfs/*.gatk.vcf.gz &&

if [ $7 == "T"]
then
    rm $1/indiv_scaff_gvcfs/*.rb.g.vcf.gz $1/scaff_vcfs/*.gatk.vcf.gz
fi

echo "Indexing concatenated vcf..."
bcftools index $1/gatk.all_indivs.vcf.gz
$2/gatk IndexFeatureFile \
     -I $1/gatk.all_indivs.vcf.gz

echo "Selecting SNPs..."
$2/gatk SelectVariants \
    -V $1/gatk.all_indivs.vcf.gz \
    -select-type SNP \
    -O $1/gatk.snp.vcf.gz \
    --use-jdk-inflater \
    --use-jdk-deflater \
    --tmp-dir $1

echo "Creating finalized VCF with iterative soft filtering of SNPs using following filters:"
echo "  multiallelic sites"
echo "  per-individual genotype read depth < ${minDP} & > ${maxDP}"
echo "  QUAL < 10.4139 (DRAGENHardQUAL filter) - (1)DRAGENHardQUAL"
echo "  % individual missing > 20% - (2)MISSING0.2"
echo "  invariable sites - (3)INVARIABLE"
echo "  minimum allele frequency < ${minmaf} - (4)MINMAF${minmaf}"
echo "  minimum mapping quality < ${minMQ} - (5)MINMQ${minMQ}"
echo "  minimum quality by depth < 2 - (6)MINQD2"
echo "  maximum strand odds ratio > 3 - (7)MAXSOR3"
if [ $8 ]
then
    echo "  masked regions - (8)MASK"
    echo "filtering..."
    bcftools view --threads 16 -m2 -M2 -v snps $1/gatk.snp.vcf.gz | \
        bcftools filter --threads 16 -S . -e "FORMAT/DP<${minDP} | FORMAT/DP>${maxDP}" -O z | \
        bcftools filter -s "(1)DRAGENHardQUAL" -e 'QUAL < 10.4139' -O z | \
        bcftools filter --threads 16 -m + -s "(2)MISSING0.2" -e 'F_MISSING > 0.2' -O z | \
        bcftools filter --threads 16 -m + -s "(3)INVARIABLE" -e 'AC==0 || AC==AN' -O z | \
        bcftools filter --threads 16 -m + -s "(4)MINMAF${minmaf}" -e "AF<${minmaf} || AF>${minmafi}" -O z | \
        bcftools filter --threads 16 -m + -s "(5)MINMQ166.67" -e "MQ < 40" -O z | \
        bcftools filter --threads 16 -m + -s "(6)MINQD2" -e "QD < 2" -O z | \
        bcftools filter --threads 16 -m + -s "(7)MAXSOR3" -e "SOR > 3" -O z | \
        bcftools filter --threads 16 -m + -s "(8)MASK" --mask-file $8 -O z -o $1/gatk.snp.filtered_soft.vcf.gz
else
    echo "filtering..."
    bcftools view --threads 16 -m2 -M2 -v snps $1/gatk.snp.vcf.gz | \
        bcftools filter --threads 16 -S . -e "FORMAT/DP<${minDP} | FORMAT/DP>${maxDP}" -O z | \
        bcftools filter -s "(1)DRAGENHardQUAL" -e 'QUAL < 10.4139' -O z | \
        bcftools filter --threads 16 -m + -s "(2)MISSING0.2" -e 'F_MISSING > 0.2' -O z | \
        bcftools filter --threads 16 -m + -s "(3)INVARIABLE" -e 'AC==0 || AC==AN' -O z | \
        bcftools filter --threads 16 -m + -s "(4)MINMAF${minmaf}" -e "AF<${minmaf} || AF>${minmafi}" -O z | \
        bcftools filter --threads 16 -m + -s "(5)MINMQ166.67" -e "MQ < 40" -O z | \
        bcftools filter --threads 16 -m + -s "(6)MINQD2" -e "QD < 2" -O z | \
        bcftools filter --threads 16 -m + -s "(7)MAXSOR3" -e "SOR > 3" -O z -o $1/gatk.snp.filtered_soft.vcf.gz
fi

echo "Creating hard filtered vcf..."
bcftools filter --threads 16 -i 'FILTER="PASS"' -O z -o $1/gatk.snp.filtered_hard.vcf.gz $1/gatk.snp.filtered_soft.vcf.gz &&
# rm $1/gatk.snp.vcf.gz*

echo "Indexing final vcfs..."
bcftools index $1/gatk.snp.filtered_soft.vcf.gz
bcftools index $1/gatk.snp.filtered_hard.vcf.gz

echo "Calculating distribution of filters..."
bcftools query --threads 16 -f '%CHROM\t%FILTER\n' $1/gatk.snp.filtered_soft.vcf.gz | \
    sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' > $1/gatk.snp.filtered_soft.FILT_DIST
echo -e "scaff\tfilter\tno.snps" > $1/gatk.snp.filtered_soft.FILT_DIST.tsv
cat $1/gatk.snp.filtered_soft.FILT_DIST | cut -f 1 | uniq | while read scaff
do
    cat $1/gatk.snp.filtered_soft.FILT_DIST | \
        grep $scaff | \
        sed -E 's/\;\(.+(\t[0-9]+)/\1/' | sed -E 's/\([0-9]+\)//' | \
        cut -f 2 | sort | uniq | while read filt
    do
        cat $1/gatk.snp.filtered_soft.FILT_DIST | \
            grep $scaff | \
            sed -E 's/MQ_filter\;//; s/RPRS_filter\;//; s/FS_SOR_filter\;//' | \
            sed -E 's/\;\(.+(\t[0-9]+)/\1/' | sed -E 's/\([0-9]+\)//' | \
            grep $filt | \
            awk -v scaff=$scaff -v filt=$filt '{sum += $3} END {print scaff"\t"filt"\t"sum}'
    done
done >> $1/gatk.snp.filtered_soft.FILT_DIST.tsv

echo "Collecting summary and per-sample metrics for final VCF..."
bcftools stats -v -d 5,100,5 -S $3 $1/gatk.snp.filtered_hard.vcf.gz > $1/gatk.snp.filtered.stats
cat $1/gatk.snp.filtered.stats | grep PSC | tail -n +2 > $1/gatk.snp.filtered.per-indiv.stats

echo "Selecting INDELs..."
$2/gatk SelectVariants \
    -V $1/gatk.all_indivs_rehead.vcf.gz \
    -select-type INDEL \
    -O $1/gatk.indel.vcf.gz
bcftools index $1/gatk.indel.vcf.gz
bcftools index -s $1/gatk.indel.vcf.gz | sort -k2,2gr | cut -f 1,3 > gatk.indel.no_snps

echo "Creating vcf with multialleleic SNPs..."
bcftools view -m3 -v snps -O z -o $1/gatk.snp.multiallelic.vcf.gz $1/gatk.snp.vcf.gz
bcftools index $1/gatk.multiallelic.vcf.gz
bcftools index -s $1/gatk.multiallelic.vcf.gz | sort -k2,2gr | cut -f 1,3 > gatk.multialleleic.no_snps


(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
