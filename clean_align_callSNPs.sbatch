#!/bin/bash

#SBATCH --job-name=callSNPs_pipeline
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=1MB
#SBATCH --output=./00_CLARCS_pipeline_%A.out
#SBATCH --error=./00_CLARCS_pipeline_%A.err

# by Samantha L.R. Capel <Samantha.Capel@Wildlife.ca.gov> <slr.capel2@gmail.com>

###############
# DESCRIPTION #
###############
# This pipeline performs whole-genome SNP calling of whole genome resequencing data following the DRAGEN-GATK best practices protocol (https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939) 
#  with the exception of extra preprocessing utilizing HTStream. Preprocessing laregely follows the recommendations of the UC Davis Bioinformatics Core 
#  (https://ucdavis-bioinformatics-training.github.io/2020-mRNA_Seq_Workshop/data_reduction/01-preproc_htstream_mm) with the addition of filtering reads with short (< 100bp) insert lengths, a necessary 
#  step for avoidance of segmenation fault errors for DRAGMAP v1.2.1. This pipeline works with samples sequenced on individual lanes or samples sequenced across multiple lanes. Restriction of SNP-calling 
#  to specific genomic intervals as well as masking of unwanted regions is also possible by specifying whitelist or blacklist files (see below). This pipeline was developed using the job scheduling 
#  manager SLURM on an HPC cluster.


################
# DEPENDENCIES #
################
# HTStream       (https://s4hts.github.io/HTStream/#hts_QWindowTrim; conda install -c bioconda htstream)
# DRAGMAP v1.2.1 (https://github.com/Illumina/DRAGMAP; conda install -c bioconda dragmap=1.2.1)
# picard         (https://github.com/broadinstitute/picard; conda install -c bioconda picard)
# GATK4          (https://github.com/broadinstitute/gatk; conda install -c bioconda gatk4)
# samtools       (https://github.com/samtools/samtools; conda install -c bioconda samtools)
# bedtools       (https://github.com/arq5x/bedtools2; conda install -c bioconda bedtools)
# bcftools       (https://github.com/samtools/bcftools; conda install -c bioconda bcftools)
## NOTE: conda is used to load DRAGMAP & picard


##################
# REQUIRED FILES #
##################
# raw fastq files
# reference genome fasta
# sample list(s): file containing sample prefixes (one per line); if reads need to be merged across lanes create two separate prefix lists - one with and one without lane designation
# read group information file: tab delimited file containing the following read group information for each (unmerged) sample (one per line) in the following order - ID SM LB PL PU
## ID = sample read group ID (e.g. sample prefix including lane designation; identical to SM if sequenced on a single lane)
## SM = sample name (excluding lane designation)
## LB = DNA prep library identifyier (only imporatant if multiple libraries were sequenced)
## PU = platform unit - code(s) for the sequencing unit/machine and lane used; if the same sequencing unit was used for all libraries just put lane designation (e.g. L001)
# OPTIONAL: bed specifying genomic interval(s) to include OR bed file specifying genomic interval(s) to exclude/mask


#################
# SET VARIABLES #
#################
dir=""  # path to root directory for all input and output files

dirconda=".../etc/profile.d"  # path to directory containing conda.sh
envconda=""  # conda environment containing DRAGMAP & picard (& fastq-filter?)
gatk=""  # path and name of gatk wrapper (e.g. /software/gatk/4.3.0.0/static/gatk-4.3.0.0)

ref=""  # path to reference genome
mtDNA=""  # name of mtDNA scaffold if known
include=""  # OPTIONAL bed file with intervals to detect SNPs within; for whole genome leave blank
exclude=""  # OPTIONAL path to masking bed file for SNP calling

raw=""  # directory containing raw reads
samples=""  # path to file containing fastq sample prefixes (including individual lane designations if applicable)
samples_merged=""  # path to file containing merged sample prefixes (i.e. excluding lane designations); if samples sequenced on single lane set to same file path as ${samples} 
rg=""  # read group information file; see README.md for details

## SET THE FOLLOWING TO 'T' OR 'F' to delete temporary and intermediate files##
rmtemp="" # set to T to delete the directory containing temporary files (i.e. ${dir}/temp)
rmhts=""  # set to T to delete cleaned FASTQs (i.e. ${dir}/01_HTS_PreProc/*fastq)
rmbam=""  # set to T to delete bams (i.e. ${dir}/02_AlignDRAGMAP/*.bam)
rmivcf="" # set to T to delete gvcfs and individual scaffold vcfs (i.e. ${dir}/04_GATKvcfs/indiv_scaff_gvcfs AND ${dir}/04_GATKvcfs/scaff_vcfs)

## IF READS NEED TO BE MERGED ACROSS LANES ##
merge=""  # set to any string if fastqs need to be merged across multiple lanes, if not leave empty
lane=""  # shared lane designation prefix (e.g. "_L00" for sample_prefix_L001, sample_prefix_L002, etc.)

## NOTE: time limits and/or memory allocation may need to be adjusted depending on the size of your data


############
# PIPELINE #
############
temp="${dir}/temp"
std="${dir}/slurmout"
preproc="${dir}/01_HTS_PreProc"
align="${dir}/02_AlignDRAGMAP"
stats="${dir}/03_AlignStats"
vcfs="${dir}/04_GATKvcfs"
[[ -d ${std}/01_htstream ]] || mkdir -p ${std}/01_htstream
[[ -d ${std}/02_DRAGMAP ]] || mkdir -p ${std}/02_DRAGMAP
[[ -d ${std}/03_merge ]] || mkdir -p ${std}/03_merge
[[ -d ${std}/04_stats ]] || mkdir -p ${std}/04_stats
[[ -d ${std}/05_hapcall ]] || mkdir -p ${std}/05_hapcall
[[ -d ${std}/06_genotype ]] || mkdir -p ${std}/06_genotype
[[ -d ${dir}/temp ]] || mkdir -p ${dir}/temp
[[ -d ${preproc} ]] || mkdir -p ${preproc}
[[ -d ${align} ]] || mkdir -p ${align}
[[ -d ${stats} ]] || mkdir -p ${stats}
[[ -d ${vcfs} ]] || mkdir -p ${vcfs}
[[ -d ${vcfs}/dragstr_model ]] || mkdir -p ${vcfs}/dragstr_model
[[ -d ${vcfs}/indiv_scaff_gvcfs ]] || mkdir -p ${vcfs}/indiv_scaff_gvcfs
[[ -d ${vcfs}/scaff_vcfs ]] || mkdir -p ${vcfs}/scaff_vcfs
[[ -d ${vcfs}/genomicDB ]] || mkdir -p ${vcfs}/genomicDB

total=$(cat $samples | wc -l)
nindiv=$(cat $samples_merged | wc -l)
echo -e "${dirconda}\t${envconda}" > ${temp}/conda.txt
module load samtools

echo "Indexing reference genome..."
samtools faidx ${ref}

echo "Creating lists for job arrays..."
if [ ${include} ]
then
    cat ${include} | cut -f 1 | sort | uniq > ${temp}/scaffolds.txt
else
    cat ${ref}.fai | cut -f 1 > ${temp}/scaffolds.txt
fi

cat $samples_merged | while read indiv
do
    cat ${temp}/scaffolds.txt | while read scaff
    do
	echo -e "${indiv}\t${scaff}"
    done
done > ${temp}/indiv_scaff.tsv

cat ${rg} | cut -f 2-3 | sort | uniq | while read line
do
    i=$(echo $line | awk '{print $1}')
    l=$(echo $line | awk '{print $2}')
    cat ${temp}/indiv_scaff.tsv | grep -w ${i} | awk -v lib=$l '{print $0"\t"lib}'
done > ${temp}/temp.tsv && mv ${temp}/temp.tsv ${temp}/indiv_scaff.tsv

echo "Creating bed files for individual scaffolds..."
cat ${ref}.fai | while read line
do
    s=$(echo $line | awk '{print $1}')
    l=$(echo $line | awk '{print $2-1}')
    cat ${ref}.fai | grep -w $s | awk -v scaff=$s -v len=$l '{print scaff"\t0\t"len}' > ${temp}/$s.bed
done

nscaffs=$(cat ${temp}/indiv_scaff.tsv | cut -f 2 | sort | uniq | wc -l)
ind_scaf=$((${nindiv}*${nscaffs}))

#                 #
## CLEAN & ALIGN ##
#                 #
# Trim & clean raw reads
htsid=$(sbatch -t 05-00:00:00 \
	       --array=1-${total} \
               --mem=200MB \
               --output=${std}/01_htstream/01_htstream_%A_%a.out \
               --error=${std}/01_htstream/01_htstream_%A_%a.err \
               HTS_preproc.slurm ${samples} ${preproc} ${raw} | sed 's/Submitted batch job //')

# Create a DRAGMAP hash table for the reference genome
htid=$(sbatch -t 10-00:00:00 \
	      --mem=85GB \
	      --output=${std}/dragmap_hash_%A.out \
	      --error=${std}/dragmap_hash_%A.err \
	      hashDRAGMAP.slurm ${temp}/conda.txt ${ref} | sed 's/Submitted batch job //')

# Align reads to genome, assign read groups, and mark duplicates
dmid=$(sbatch -t 05-00:00:00 \
	      --array=1-${total} \
              --mem=45GB \
              --output=${std}/02_DRAGMAP/02_DRAGMAP_%A_%a.out \
              --error=${std}/02_DRAGMAP/02_DRAGMAP_%A_%a.err \
	      --dependency=afterok:${htid}:${htsid} \
              alignDRAGMAP.slurm ${samples} ${align} ${ref} ${temp} ${preproc} ${rg} ${gatk} ${lane} ${rmhts} | sed 's/Submitted batch job //')

# Break genome into 50 Kb windows for aligment stats
gwid=$(sbatch -t 01-00:00:00 \
              --mem=1G \
              --output=${std}/genome_win_%A.out \
              --error=${std}/genome_win_%A.err \
              genome_wins.slurm ${ref} | sed 's/Submitted batch job //')

if [ ${merge} ]
then
   # Merge files across lanes
   mgid=$(sbatch -t 03-00:00:00 \
   		 --array=1-${nindiv} \
   	         --mem=100MB \
   	         --output=${std}/03_merge/03_merge_%A_%a.out \
   	         --error=${std}/03_merge/03_merge_%A_%a.err \
		 --dependency=afterok:${dmid} \
   	         samtools_merge.slurm ${samples_merged} ${align} ${lane} | sed 's/Submitted batch job //')

   # Calculate alignment statistics
   asid=$(sbatch -t 3-00:00:00 \
	      --array=1-${nindiv} \
              --mem=5G \
              --output=${std}/04_stats/04_stats_%A_%a.out \
              --error=${std}/04_stats/04_stats_%A_%a.err \
	      --dependency=afterok:${mgid}:${gwid} \
	      align_stats.slurm ${samples_merged} ${stats} ${ref} ${align} ${gatk} | sed 's/Submitted batch job //')
else
    asid=$(sbatch -t 3-00:00:00 \
	      --array=1-${nindiv} \
              --mem=5G \
              --output=${std}/04_stats/04_stats_%A_%a.out \
              --error=${std}/04_stats/04_stats_%A_%a.err \
	      --dependency=afterok:${dmid}:${gwid} \
	      align_stats.slurm ${samples_merged} ${stats} ${ref} ${align} ${gatk} | sed 's/Submitted batch job //')
fi

#             #
## CALL SNPS ##
#             #
# Calibrate DRAGEN STR model
stid=$(sbatch -t 02-00:00:00 \
              --mem=1GB \
              --output=${std}/str_table_%A.out \
              --error=${std}/str_table_%A.err \
	      STRtable.slurm ${ref} ${temp}/conda.txt ${gatk} ${include} | sed 's/Submitted batch job //')

# Call individual haplotypes
htcid=$(sbatch -t 10-00:00:00 \
              --mem=2GB \
              --array=1-${ind_scaf}%1000 \
              --error=${std}/05_hapcall/05_hapcall_%A_%a.err \
              --output=${std}/05_hapcall/05_hapcall_%A_%a.out \
	      --dependency=afterok:${stid}:${asid} \
	      bam_to_gvcf.slurm ${temp} ${ref} ${gatk} ${align} ${vcfs} ${rmbam} ${include} | sed 's/Submitted batch job //')

# Joint genotype individuals by scaffold
jgtid=$(sbatch -t 10-00:00:00 \
               --mem=5GB \
               --array=1-${nscaffs} \
               --error=${std}/06_genotype/06_genotype_%A_%a.err \
               --output=${std}/06_genotype/06_genotype_%A_%a.out \
	       --dependency=afterok:${htcid} \
	       gvcf_to_vcf_scaff.slurm ${temp} ${vcfs} ${gatk} ${ref} ${samples_merged} | sed 's/Submitted batch job //')

# Concatenate vcfs & filter SNPs
csid=$(sbatch -t 05-00:00:00 \
	      --mem=5GB \
	      --error=${std}/07_concat_vcfs_filter_SNPs_%A.err \
	      --output=${std}/07_concat_vcfs_filter_SNPs_%A.out \
	      --dependency=afterok:${jgtid} \
	      vcf_scaff_to_snp.vcf.slurm ${vcfs} ${gatk} ${samples_merged} ${mtDNA} ${nindiv} ${rmivcf} ${exclude} | sed 's/Submitted batch job //')

# Report % successful jobs
echo -e "${htsid}\n${htid}\n${dmid}\n${gwid}\n${asid}\n${stid}\n${htcid}\n${jgtid}\n${csid}\n${mgid}" > ${temp}/job_IDs.txt
sucid=$(sbatch -t 01-00:00:00 \
	  --mem=1MB \
	  --error=${std}/JOB_REPORT_%A.err \
	  --output=${std}/JOB_REPORT_%A.out \
	  --dependency=after:${csid} \
	  job_report.slurm ${std} ${temp} ${rmtemp} | sed 's/Submitted batch job //')

(exit) && echo success
