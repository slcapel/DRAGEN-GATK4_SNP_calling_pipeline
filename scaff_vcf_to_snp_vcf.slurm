#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

dpmin=$3
dpmax=$4
minmaf=$(echo $5 | awk '{print 5/(2*$0)}')
minmafi=$(echo $5 | awk '{print 1-5/(2*$0)}')

module load bcftools
if [[ ! -f $1/gatk.all_indivs.vcf.gz.csi ]]
then
    echo "Concatenating vcfs..."
    bcftools concat -o $1/gatk.all_indivs.vcf.gz -O z $1/scaff_vcfs/*.gatk.vcf.gz

    echo "Indexing concatenated vcf..."
    bcftools index $1/gatk.all_indivs.vcf.gz
    gatk IndexFeatureFile -I $1/gatk.all_indivs.vcf.gz
fi

if [ $6 == "T" ]
then
    rm -r $1/indiv_scaff_gvcfs $1/scaff_vcfs
fi

if [[ ! -f $1/gatk.indel.vcf.gz ]]
then
    module load gatk
    echo "Selecting SNPs..."
    gatk SelectVariants \
         -V $1/gatk.all_indivs.vcf.gz \
         -select-type SNP \
         -O $1/gatk.snp.vcf.gz

    echo "Selecting INDELs..."
    gatk SelectVariants \
         -V $1/gatk.all_indivs.vcf.gz \
         -select-type INDEL \
         -O $1/gatk.indel.vcf.gz

    module load bcftools
    echo "Creating vcf with multialleleic SNPs..."
    bcftools view -m3 -v snps -O z -o $1/gatk.snp.multiallelic.vcf.gz $1/gatk.snp.vcf.gz

    echo "Indexing SNP vcf..."
    bcftools index $1/gatk.snp.vcf.gz
fi

module load vcftools
module load bcftools
if [ ! -f "$1/gatk.snps.sub1M.GT_DP.tsv" ] || [ $(wc -l < "$1/gatk.snps.sub1M.GT_DP.tsv") -lt 2 ]
then
    echo "Subsampling pre-filtering quality metrics for QC assessment..."
    bcftools view --threads 8 -H $1/gatk.snp.vcf.gz | cut -f 1,2 | shuf | head -n 1000000 | sort -k1,1 -k2,2g > $1/gatk.snps.sub1M.pos
    echo -e "QD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum" > $1/gatk.snps.sub1M.qual.tsv
    bcftools view --threads 8 -R $1/gatk.snps.sub1M.pos -O u $1/gatk.snp.vcf.gz | bcftools query -f '%QD\t%FS\t%SOR\t%MQ\t%MQRankSum\t%ReadPosRankSum\n' >> $1/gatk.snps.sub1M.qual.tsv
    echo -e "SAMPLE\tGT\tDP" > $1/gatk.snps.sub1M.GT_DP.tsv
    bcftools view --threads 8 -R $1/gatk.snps.sub1M.pos -O u $1/gatk.snp.vcf.gz | bcftools query -f '[%SAMPLE\t%GT\t%DP\n]' | sed -E 's/\|/\//' >> $1/gatk.snps.sub1M.GT_DP.tsv
fi

if [[ ! -f $8/raw_SNP_GT_read_depth_distribution.svg ]]
then
    module load R
    Rscript raw_SNP_QC_metric_dists.R $8 $1/gatk.snps.sub1M.qual.tsv
    Rscript raw_SNP_GT_read_depth_dist.R $8 $1/gatk.snps.sub1M.GT_DP.tsv
fi

if [[ ! -f $1/gatk.snp.soft_filtered.vcf.gz.csi ]]
then
    echo "Creating VCF with iterative soft filtering of SNPs using following filters:"
    echo "  QUAL < 10.4139 (DRAGENHardQUAL filter) - (1)DRAGENHardQUAL"
    echo "  % individual missing > 20% - (2)MISSING>0.2"
    echo "  invariable sites - (3)INVARIABLE"
    echo "  minor allele frequency < ${minmaf} - (4)MINMAF<${minmaf}"
    echo "  RMS mapping quality < 166.67 - (5)MQ<166.67"
    echo "  quality by depth < 2 - (6)QD<2"
    echo "  symmetric odds ratio > 3 - (7)SOR>3"
    echo "  Fisher's strandedness > 60 - (8)FS>60"
    echo "  MQ rank sum test Z-score < -12.5 - (9)MQRankSum<-12.5"
    echo "  Read position bias test Z-score < -8 - (10)ReadPosRankSum<-8"
    echo "filtering..."
    bcftools view --threads 8 -m2 -M2 -v snps -Ou $1/gatk.snp.vcf.gz | \
	bcftools filter --threads 8 -S . -e "FORMAT/DP<${dpmin} | FORMAT/DP>${dpmax}" -Ou | \
	bcftools +fill-tags --threads 8 -Ou -- -t AC,AN,AF | \
	bcftools filter --threads 8 -s "(1)DRAGENHardQUAL" -e 'QUAL < 10.4139' -Ou | \
	bcftools filter --threads 8 -m + -s "(2)MISSING0.2" -e 'F_MISSING > 0.2' -Ou | \
	bcftools filter --threads 8 -m + -s "(3)INVARIABLE" -e 'AC==0 || AC==AN' -Ou | \
	bcftools filter --threads 8 -m + -s "(4)MINMAF${minmaf}" -e "AF<${minmaf} || AF>${minmafi}" -Ou | \
	bcftools filter --threads 8 -m + -s "(5)MINMQ166.67" -e 'MQ < 166.67' -Ou | \
	bcftools filter --threads 8 -m + -s "(6)MINQD2" -e 'QD < 2' -Ou | \
	bcftools filter --threads 8 -m + -s "(7)MAXSOR3" -e 'SOR > 3' -Ou | \
	bcftools filter --threads 8 -m + -s "(8)MAXFS60" -e 'FS > 60' -Ou | \
	bcftools filter --threads 8 -m + -s "(9)MINMQRankSum-12.5" -e 'MQRankSum < -12.5' -Ou | \
	bcftools filter --threads 8 -m + -s "(10)MINReadPosRankSum-8" -e 'ReadPosRankSum < -8' -O z -o $1/gatk.snp.soft_filtered.vcf.gz
    bcftools index $1/gatk.snp.soft_filtered.vcf.gz
fi

if [[ ! -f $1/gatk.snp.hard_filtered.vcf.gz.csi ]]
then
    echo "Creating hard filtered vcf..."
    bcftools filter --threads 8 -i 'FILTER="PASS"' -O z -o $1/gatk.snp.hard_filtered.vcf.gz $1/gatk.snp.soft_filtered.vcf.gz
    bcftools index $1/gatk.snp.hard_filtered.vcf.gz
fi

if [ ! -f "$1/gatk.snp.soft_filtered.FILT_DIST.tsv" ] || [ $(wc -l < "$1/gatk.snp.soft_filtered.FILT_DIST.tsv") -lt 2 ]
then
    echo "Calculating distribution of filters..."
    bcftools query -f '%CHROM\t%FILTER\n' $1/gatk.snp.soft_filtered.vcf.gz | \
	sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' > $1/gatk.snp.soft_filtered.FILT_DIST
    echo -e "scaff\tfilter\tno.snps" > $1/gatk.snp.soft_filtered.FILT_DIST.tsv
    cat $1/gatk.snp.soft_filtered.FILT_DIST | cut -f 1 | uniq | while read scaff
    do
	cat $1/gatk.snp.soft_filtered.FILT_DIST | \
            grep $scaff | \
            sed -E 's/\;\(.+(\t[0-9]+)/\1/' | sed -E 's/\([0-9]+\)//' | \
            cut -f 2 | sort | uniq | while read filt
	do
            cat $1/gatk.snp.soft_filtered.FILT_DIST | \
		grep $scaff | \
		sed -E 's/MQ_filter\;//; s/RPRS_filter\;//; s/FS_SOR_filter\;//' | \
		sed -E 's/\;\(.+(\t[0-9]+)/\1/' | sed -E 's/\([0-9]+\)//' | \
		grep $filt | \
		awk -v scaff=$scaff -v filt=$filt '{sum += $3} END {print scaff"\t"filt"\t"sum}'
	done
    done >> $1/gatk.snp.soft_filtered.FILT_DIST.tsv
fi

if [ ! -f "$1/gatk.snp.hard_filtered.thin20Kbp_hetVAF.tsv.gz" ] || [ $(wc -l < "$1/gatk.snp.hard_filtered.thin20Kbp_hetVAF.tsv.gz") -lt 2 ]
then
    echo "Collecting summary and per-sample metrics for final VCF..."
    bcftools stats -v -d 5,100,5 -S $2 $1/gatk.snp.hard_filtered.vcf.gz > $1/gatk.snp.filtered.stats
    cat $1/gatk.snp.filtered.stats | grep -w PSC | tail -n +2 | cut -f 3- | sed -E 's/\[[0-9]+\]//g' > $1/gatk.snp.filtered.per-ind_stats.tsv
    vcftools --gzvcf $1/gatk.snp.hard_filtered.vcf.gz --thin 20000 --recode --recode-INFO-all --stdout | \
	bcftools +fill-tags -- -t VAF | \
	bcftools query -f '[%SAMPLE\t%GT\t%VAF\t%DP\n]' | \
	sed -E 's/\|/\//' | \
	grep "0/1" | \
	cut -f 1,3,4 | gzip > $1/gatk.snp.hard_filtered.thin20Kbp_hetVAF.tsv.gz
fi

module load R
Rscript filtered_SNP_QC.R $8 $1/gatk.snp.filtered.per-ind_stats.tsv $1/gatk.snp.hard_filtered.thin20Kbp_hetVAF.tsv.gz $1/gatk.snp.soft_filtered.FILT_DIST.tsv ${minmaf}

(exit) && echo success
