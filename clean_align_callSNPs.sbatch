#!/bin/bash

#SBATCH --job-name=callSNPs_pipeline
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=1MB
#SBATCH --output=/path/to/standard/output/00_CLARCS_pipeline_%A.out
#SBATCH --error=/path/to/standard/output/00_CLARCS_pipeline_%A.err

# by Samantha L.R. Capel <Samantha.Capel@Wildlife.ca.gov> <slr.capel2@gmail.com>

###############
# DESCRIPTION #
###############
# This pipeline performs whole-genome SNP calling of whole genome resequencing data following the DRAGEN-GATK best practices protocol (https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939) 
# with the exception of extra preprocessing utilizing HTStream. Preprocessing laregely follows the recommendations of the UC Davis Bioinformatics Core 
# (https://ucdavis-bioinformatics-training.github.io/2020-mRNA_Seq_Workshop/data_reduction/01-preproc_htstream_mm) with the addition of filtering reads with short (< 100bp) insert lengths, a necessary 
# step for avoidance of segmenation fault errors for DRAGMAP v1.2.1. This pipeline works with samples sequenced on individual lanes or samples sequenced across multiple lanes. Restriction of SNP-calling 
# to specific genomic intervals as well as masking of unwanted regions is also possible by specifying whitelist or blacklist bed files (see below). This pipeline was developed using the job scheduling 
# manager SLURM on an HPC cluster.


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
# R              (https://cran.r-project.org/bin/windows/base)
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
dir="/path/to/dir"  # path to root directory for all input and output files
envconda="/quobyte/yotegrp/cdfwwildlife/Capel_Dedicated/miniconda/envs/dragmap"  # conda environment containing DRAGMAP

ref="${dir}/path/to/ref.fna"  # path to reference genome fasta
include="${dir}/path/to/whitelist.bed"  # OPTIONAL bed file with intervals for final regions to include in VCF
exclude=""  # OPTIONAL bed file for masking regions for SNP calling

raw="${dir}/00_Raw"  # directory containing raw reads
samples="${dir}/Sample_Lists/all_samples_unmerged.txt"  # path to file containing fastq sample prefixes (including individual lane designations if applicable)
samples_merged="${dir}/Sample_Lists/all_samples_merged.txt"  # REQUIRED path to file containing merged sample prefixes (i.e. excluding lane designations); !!IMPORTANT!! if samples sequenced on single lane set to same file path as ${samples} 
rg="${dir}/Sample_Lists/RG_info.tsv"  # read group information file; does not have to be in same order as sample name files

dpmin=0  # minimum read depth per genotype; set to 0 to omit genotype filtering in final VCF
dpmax=1000000  # maximum read depth per genotype; set to 1000000 to omit genotype filtering

## SET THE FOLLOWING TO 'T' or 'F' to dictate intermediate file deletion ##
rmtemp="F" # delete the directory containing temporary files (i.e. ${dir}/temp)
rmhts="F"  # delete cleaned FASTQs (i.e. ${dir}/01_HTS_PreProc/*fastq)
rmbam="F"  # delete bams (i.e. ${dir}/02_AlignDRAGMAP/*.bam)
rmivcf="F" # delete gvcfs and individual scaffold vcfs (i.e. ${dir}/04_GATKvcfs/indiv_scaff_gvcfs AND ${dir}/04_GATKvcfs/scaff_vcfs)

## IF READS NEED TO BE MERGED ACROSS LANES ##
merge="T"  # set to any string if fastqs need to be merged across multiple lanes, if not leave empty
lane="_L0"  # shared lane designation prefix (e.g. "_L00" for sample_prefix_L001, sample_prefix_L002, etc.); only needed if merging bams

## NOTE: time limits and/or memory allocation may need to be adjusted depending on the size of your data


############
# PIPELINE #
############
temp="${dir}/temp"
std="${dir}/SCRIPTS/slurmout/DRAGMAP-GATK_pipeline"
preproc="${dir}/01_HTS_PreProc"
align="${dir}/02_AlignDRAGMAP"
stats="${dir}/03_AlignStats"
vcfs="${dir}/04_GATKvcfs"
figs="${dir}/Figures"
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
[[ -d ${figs} ]] || mkdir -p ${figs}

total=$(cat $samples | wc -l)
nindiv=$(cat $samples_merged | wc -l)
module load samtools

echo "Reference genome: ${ref}"

if [[ ! -f ${ref}.fai ]]
then
    echo "Indexing reference genome..."
    samtools faidx ${ref}
fi

if [[ ! -f ${temp}/indiv_scaff.tsv ]]
then
    echo "Creating lists for job arrays..."
    if [[ -n ${include} ]]
    then
	cat ${include} | cut -f 1 > ${temp}/scaffolds.txt
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

    echo "Creating bed files for individual scaffolds..."
    cat ${ref}.fai | while read line
    do
	s=$(echo $line | awk '{print $1}')
	l=$(echo $line | awk '{print $2-1}')
	cat ${ref}.fai | grep -w $s | awk -v scaff=$s -v len=$l '{print scaff"\t0\t"len}' > ${temp}/$s.bed
    done
fi

nscaffs=$(cat ${temp}/indiv_scaff.tsv | cut -f 2 | sort | uniq | wc -l)
ind_scaf=$((${nindiv}*${nscaffs}))

echo -e "total\t$total"
echo -e "nindiv\t$nindiv"
echo -e "nscaffs\t$nscaffs"
echo -e "ind_scaf\t$ind_scaf"

echo -e "job_name\tscript_name\tjob_ID" > ${temp}/job_IDs_${SLURM_JOB_ID}.txt


#                 #
## CLEAN & ALIGN ##
#                 #
# Trim & clean raw reads
htsid=$(sbatch --account=yotegrp \
	       --job-name=htstream \
	       -t 02-00:00:00 \
	       --array=1-${total}%500 \
               --mem=200MB \
               --output=${std}/01_htstream/01_htstream_%A_%a.out \
               --error=${std}/01_htstream/01_htstream_%A_%a.err \
               HTS_preproc.slurm ${samples} ${preproc} ${raw} | sed 's/Submitted batch job //')
echo -e "htstream\tHTS_preproc.slurm\t${htsid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

# Create a DRAGMAP hash table for the reference genome
htid=$(sbatch --account=yotegrp \
	      --job-name=drgmp_hash \
              -t 03-00:00:00 \
	      --mem=85GB \
	      --output=${std}/dragmap_hash_%A.out \
	      --error=${std}/dragmap_hash_%A.err \
	      hashDRAGMAP.slurm ${envconda} ${ref} | sed 's/Submitted batch job //')
echo -e "drgmp_hash\thashDRAGMAP.slurm\t${htid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

# Align reads to genome, assign read groups, and mark duplicates
dmid=$(sbatch --account=yotegrp \
	      -t 03-00:00:00 \
	      --job-name=DRAGMAP \
	      --array=1-${total}%100 \
              --mem=60GB \
	      --dependency=afterok:${htid}:${htsid} \
              --output=${std}/02_DRAGMAP/02_DRAGMAP_%A_%a.out \
              --error=${std}/02_DRAGMAP/02_DRAGMAP_%A_%a.err \
              alignDRAGMAP.slurm ${samples} ${align} ${ref} ${temp} ${preproc} ${rg} ${envconda} ${rmhts} ${exclude} | sed 's/Submitted batch job //')
echo -e "DRAGMAP\talignDRAGMAP.slurm\t${dmid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

# Break genome into 50 Kb windows for aligment stats
gwid=$(sbatch --account=yotegrp \
	      --job-name=genome_win \
              -t 01-00:00:00 \
              --mem=1G \
              --output=${std}/genome_win_%A.out \
              --error=${std}/genome_win_%A.err \
              genome_wins.slurm ${ref} | sed 's/Submitted batch job //')
echo -e "genome_win\tgenome_wins.slurm\t${gwid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

if [ ${merge} ]
then
   # Merge files across lanes
    mgid=$(sbatch --account=yotegrp \
		  --job-name=bam_merge \
	          -t 03-00:00:00 \
		  --array=1-${nindiv} \
	          --mem=100MB \
	          --output=${std}/03_merge/03_merge_%A_%a.out \
	          --error=${std}/03_merge/03_merge_%A_%a.err \
		  --dependency=afterok:${dmid} \
	          samtools_merge.slurm ${samples_merged} ${align} ${lane} | sed 's/Submitted batch job //')
    echo -e "bam_merge\tsamtools_merge.slurm\t${mgid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

   # Calculate alignment statistics
    asid=$(sbatch --account=yotegrp \
		  --job-name=align_stats \
	          -t 03-00:00:00 \
	          --array=1-${nindiv} \
                  --mem=5G \
		  --dependency=afterok:${mgid}:${gwid} \
                  --output=${std}/04_stats/04_stats_%A_%a.out \
                  --error=${std}/04_stats/04_stats_%A_%a.err \
	          align_stats.slurm ${samples_merged} ${stats} ${ref} ${align} | sed 's/Submitted batch job //')
    echo -e "align_stats\talign_stats.slurm\t${asid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt
	  
else
    asid=$(sbatch --account=yotegrp \
		  --job-name=align_stats \
	          -t 03-00:00:00 \
	          --array=1-${nindiv} \
                  --mem=5G \
                  --output=${std}/04_stats/04_stats_%A_%a.out \
                  --error=${std}/04_stats/04_stats_%A_%a.err \
	          --dependency=afterok:${dmid}:${gwid} \
	          align_stats.slurm ${samples_merged} ${stats} ${ref} ${align} | sed 's/Submitted batch job //')
    echo -e "align_stats\talign_stats.slurm\t${asid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt
fi

#             #
## CALL SNPS ##
#             #
# Calibrate DRAGEN STR model
stid=$(sbatch --account=yotegrp \
	      --job-name=str_table \
              -t 02-00:00:00 \
              --mem=5GB \
              --output=${std}/str_table_%A.out \
              --error=${std}/str_table_%A.err \
	      STRtable.slurm ${ref} ${envconda} | sed 's/Submitted batch job //')
echo -e "str_table\tSTRtable.slurm\t${stid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

# Call individual haplotypes
htcid=$(sbatch --account=yotegrp \
	       --job-name=GATK4hapcall \
	       -t 03-00:00:00 \
               --mem=10GB \
               --array=1-${ind_scaf}%500 \
	       --dependency=afterok:${asid} \
               --error=${std}/05_hapcall/05_hapcall_%A_%a.err \
               --output=${std}/05_hapcall/05_hapcall_%A_%a.out \
	       bam_to_gvcf.slurm ${temp} ${ref} ${align} ${vcfs} ${rmbam} | sed 's/Submitted batch job //')
echo -e "GATK4hapcall\tbam_to_gvcf.slurm\t${htcid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

# Joint genotype individuals by scaffold
jgtid=$(sbatch --account=yotegrp \
	       --job-name=GATK4genotype \
	       -t 10-00:00:00 \
               --mem=5GB \
               --array=1-${nscaffs} \
               --error=${std}/06_genotype/06_genotype_%A_%a.err \
               --output=${std}/06_genotype/06_genotype_%A_%a.out \
	       --dependency=afterok:${htcid} \
	       gvcf_to_vcf_scaff.slurm ${temp} ${vcfs} ${ref} ${samples_merged} | sed 's/Submitted batch job //')
echo -e "GATK4genotype\tgvcf_to_vcf_scaff.slurm\t${jgtid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

# Concatenate vcfs & filter SNPs
csid=$(sbatch --account=yotegrp \
	      --job-name=gvcf2snp \
              -t 03-00:00:00 \
	      --mem=25GB \
	      --dependency=afterok:${jgtid} \
	      --error=${std}/07_concat_vcfs_filter_SNPs_%A.err \
	      --output=${std}/07_concat_vcfs_filter_SNPs_%A.out \
	      scaff_vcf_to_snp_vcf.slurm ${vcfs} ${samples_merged} ${dpmin} ${dpmax} ${nindiv} ${rmivcf} ${ref} ${figs} | sed 's/Submitted batch job //')
echo -e "gvcf2snp\tscaff_vcf_to_snp_vcf.slurm\t${csid}" >> ${temp}/job_IDs_${SLURM_JOB_ID}.txt

# Report % successful jobs
cat ${temp}/job_IDs_${SLURM_JOB_ID}.txt
sucid=$(sbatch --account=yotegrp \
	       --job-name=jobrep \
	       -t 01-00:00:00 \
	       --mem=1MB \
	       --error=${std}/JOB_REPORT_%A.err \
	       --output=${std}/JOB_REPORT_%A.out \
	       --dependency=afterany:${csid} \
	       job_report.slurm ${std} ${temp} ${SLURM_JOB_ID} ${rmtemp}  | sed 's/Submitted batch job //')

(exit) && echo success
