#!/bin/bash

#SBATCH --job-name=gvcf2snp
#SBATCH --nodes=1

start=`date +%s`
echo $HOSTNAME

module load bcftools
minmaf=$(echo $5 | awk '{print 5/(2*$0)}')
minmafi=$(echo $5 | awk '{print 1-5/(2*$0)}')

echo "Concatenating vcfs..."
bcftools concat -o $1/gatk.all_indivs.vcf.gz -O z $1/scaff_vcfs/*.gatk.vcf.gz &&

if [ $6 == "T"]
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
    -O $1/gatk.snp.vcf.gz

echo "Selecting INDELs..."
$2/gatk SelectVariants \
    -V $1/gatk.all_indivs.vcf.gz \
    -select-type INDEL \
    -O $1/gatk.indel.vcf.gz

echo "Creating mtDNA-specific VCF..."
bcftools filter -r $4 -O z -o $1/gatk.mtDNA.vcf.gz $1/gatk.all_indivs.vcf.gz

echo "Creating vcf with multialleleic SNPs..."
bcftools view -m3 -v snps -O z -o $1/gatk.snp.multiallelic.vcf.gz $1/gatk.snp.vcf.gz

echo "Creating finalized VCF with iterative soft filtering of SNPs using following filters:"
echo "  multiallelic sites"
echo "  mtDNA sites"
echo "  per-individual genotype read depth < 10 & > 200"
echo "  QUAL < 10.4139 (DRAGENHardQUAL filter) - (1)DRAGENHardQUAL"
echo "  % individual missing > 20% - (2)MISSING0.2"
echo "  invariable sites - (3)INVARIABLE"
echo "  minimum allele frequency < ${minmaf} - (4)MINMAF${minmaf}"
if [ $7 ]
then
    echo "  masked regions - (5)MASK"
    echo "filtering..."
    bcftools view -m2 -M2 -v snps $1/gatk.snp.vcf.gz | bcftools filter -m $4 -O z | bcftools filter -S . -e 'FORMAT/DP<10 | FORMAT/DP>200' -O z |  bcftools filter -s "(1)DRAGENHardQUAL" -e '%QUAL<10.4139' -O z | bcftools filter -m + -s "(2)MISSING0.2" -e 'F_MISSING > 0.2' -O z | bcftools filter -m + -s "(3)INVARIABLE" -e 'AC==0 || AC==AN' -O z | bcftools filter -m + -s "(4)MINMAF${minmaf}" -e "AF<${minmaf} || AF>${minmafi}" -O z | bcftools filter -m + -s "(5)MASK" --mask-file $6 -O z -o $1/gatk.snp.soft_filtered.vcf.gz
else
    echo "filtering..."
    bcftools view -m2 -M2 -v snps $1/gatk.snp.vcf.gz | bcftools filter -m $4 -O z | bcftools filter -S . -e 'FORMAT/DP<10 | FORMAT/DP>200' -O z | bcftools filter -s "(1)DRAGENHardQUAL" -e '%QUAL<10.4139' -O z | bcftools filter -m + -s "(2)MISSING0.2" -e 'F_MISSING > 0.2' -O z | bcftools filter -m + -s "(3)INVARIABLE" -e 'AC==0 || AC==AN' -O z | bcftools filter -m + -s "(4)MINMAF${minmaf}" -e "AF<${minmaf} || AF>${minmafi}" -O z -o $1/gatk.snp.soft_filtered.vcf.gz
fi

echo "Creating hard filtered vcf..."
bcftools filter -i 'FILTER="PASS"' -O z -o $1/gatk.snp.hard_filtered.vcf.gz $1/gatk.snp.soft_filtered.vcf.gz &&
rm $1/gatk.snp.vcf.gz*

echo "Indexing final vcfs..."
bcftools index $1/gatk.snp.soft_filtered.vcf.gz
bcftools index $1/gatk.snp.hard_filtered.vcf.gz

echo "Calculating distribution of filters..."
bcftools query -f '%CHROM\t%FILTER\n' $1/gatk.snp.soft_filtered.vcf.gz | sort -g | uniq -c > $1/gatk.snp.filtered.FILTER_DIST

echo "Collecting summary and per-sample metrics for final VCF..."
if [ $7 ]
then
    bcftools filter -i 'FILTER="PASS" || FILTER="(5)MASK"' $1/gatk.snp.soft_filtered.vcf.gz | bcftools stats -v -d 5,100,5 -S $3 > $1/gatk.snp.filtered+masked.stats
fi
bcftools stats -v -d 5,100,5 -S $3 $1/gatk.snp.hard_filtered.vcf.gz > $1/gatk.snp.filtered.stats

if [ $7 ]
then
   echo -e "Sample\tnRefHom\tnNonRefHom\tnHets\tmean_depth\tnMissing\tnSites" > $1/gatk.snp.filtered+masked.per-indiv.stats
cat $1/gatk.snp.filtered+masked.stats | grep ^PSC | cut -f 3-6,10,14 | awk '{print $0"\t"$2+$3+$4+$6}' >> $1/gatk.snp.filtered+masked.per-indiv.stats
fi
echo -e "Sample\tnRefHom\tnNonRefHom\tnHets\tmean_depth\tnMissing\tnSites" > $1/gatk.snp.filtered.per-indiv.stats
cat $1/gatk.snp.filtered.stats | grep ^PSC | cut -f 3-6,10,14 | awk '{print $0"\t"$2+$3+$4+$6}' >> $1/gatk.snp.filtered.per-indiv.stats

if [ $7 ]
then
    cat $1/gatk.snp.filtered+masked.per-indiv.stats | sed -E 's/Sites/Sites\tSet/' | awk '{if ($1 != "Sample") print $0"\tUnmasked"; else print $0}' > temp.tsv && mv temp.tsv $1/gatk.snp.filtered+masked.per-indiv.stats
    tail -n +2 $1/gatk.snp.filtered.per-indiv.stats | awk '{print $0"\tMasked"}' >> $1/gatk.snp.filtered+masked.per-indiv.stats
fi

(exit) && echo success
end=`date +%s`
runtime=$((end-start))
echo "RUNTIME: $runtime"
