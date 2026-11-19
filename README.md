# WGR DRAGEN-GATK SNP-Calling Pipeline
This pipeline performs whole-genome SNP calling of paired-end whole genome resequencing data following the [DRAGEN-GATK best practices protocol](https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939) with the exception of extra preprocessing utilizing HTStream. Preprocessing laregely follows the [recommendations of the UC Davis Bioinformatics Core](https://ucdavis-bioinformatics-training.github.io/2020-mRNA_Seq_Workshop/data_reduction/01-preproc_htstream_mm) with the addition of filtering reads with short (< 100bp) insert lengths, necessary for avoiding segmenation fault errors in DRAGMAP v1.2.1. This pipeline works with samples sequenced on a single lane or samples sequenced across multiple lanes. Restriction of SNP-calling to specific genomic intervals as well as masking of unwanted regions is also possible by specifying whitelist or blacklist bed files. This pipeline was developed on an HPC cluster utilizing [SLURM](https://slurm.schedmd.com/quickstart.html).

<p align="center">
<img width="1134" height="529" alt="Image" src="https://github.com/user-attachments/assets/c60afee1-395e-4688-89b8-1ce87b601022">
</p>

## Dependencies:
- [HTStream](https://s4hts.github.io/HTStream/#hts_QWindowTrim)
- [DRAGMAP v1.2.1](https://github.com/Illumina/DRAGMAP)
- [picard](https://github.com/broadinstitute/picard)
- [GATK4](https://github.com/broadinstitute/gatk)
- [samtools](https://github.com/samtools/samtools)
- [bedtools2](https://github.com/arq5x/bedtools2)
- [bcftools](https://github.com/samtools/bcftools)
- [R](https://www.r-project.org/)

NOTE: scripts use conda to load DRAGMAP & picard; HTStream, samtools, bedtools, and bcftools are all called using `module load` but can be incorporated into the conda environment specified

## Required Input Files:
- raw unmapped FASTQ files
- **unzipped** reference genome FASTA
- sample list(s): file containing sample prefixes (one per line); if individuals were sequenced across multiple lanes, create two separate prefix lists - one with and one without lane designations
    | Unmerged | Merged |
    |:----------:|:--------:|
    | <pre>sample1_L001<br>sample1_L002<br>sample2_L001<br>sample2_L002<br>sample3_L001<br>sample3_L002<br></pre> | <pre>sample1<br>sample2<br>sample3<br><br><br><br></pre> |
- read group information file: tab delimited file containing the following read group information for each pair of FASTQ files (one per line) in the following order - ID SM LB PL PU
    - ID = sample read group ID (e.g. sample prefix including lane designation; identical to SM if sequenced on a single lane)
    - SM = sample name (excluding lane designation)
    - LB = DNA library identifyier (only imporatant if multiple libraries were sequenced)
    - PU = platform uint; this can indicate if samples were ran on different lanes and/or different sequencing units
- OPTIONAL: bed file specifying genomic interval(s) to include OR bed file specifying genomic interval(s) to exclude/mask

## Instructions:
1. Download all scripts into a single directory
2. Open `clean_align_callSNPs.sbatch` in a text editor and set **ALL** required variables in the designated section (starting on line 53)
3. Execute the pipeline by running `sbatch clean_align_callSNPs.sbatch`

## Final outputs:
### Job information outputs (found in specified slum output directory)
- Success report for all jobs executed -- **JOB_REPORT_<*jobID*>.out**
- List of all job names and IDs -- **job_IDs_<*jobID*>.txt**
- Standard output and error for each job -- **<*job_name*>_<*jobID*>.out** & **<*job_name*>_<*jobID*>.err**
### DRAGMAP alignment statistics for each sample (found in 03_AlignStats/<*sampleID*>)
- Summary of alignment flags -- ***sampleID*.dragmap_flagstats.txt** (via `samtools flagstat`)
- Comprehensive alignment stats -- ***sampleID*.dragmap_stats.txt** (via `samtools stats`; can be plotted with `plot-bamstats`)
- Per-scaffold numer of mapped (column 1) and unmapped (column 2) reads -- ***sampleID*.dragmap_IDstats.tsv** (via `samtools idxstats`)
- Genome-wide mean depth of coverage -- ***sampleID*.dragmap_depth.txt** (via `samtools depth`)
- Per-scaffold summary alignment stats -- ***sampleID*.dragmap_coverage.tsv** (via `samtools coverage`)
- GATK alignment summary stats -- ***sampleID*.dragmap.gatk.stat.txt** (via `gatk CollectAlignmentSummaryMetrics`)
- Mean depth of coverage by 50 kb windows with 10% overlap -- ***sampleID*.dragmap_50kbwin_meancov.tsv** (via `bedtools coverage`)
### VCFs (found in 04_GATKvcfs)
- Unfiltered -- **gatk.all_indivs.vcf.gz**
- All raw SNPs -- **gatk.snp.vcf.gz**
- Insertions and deletions -- **gatk.indel.vcf.gz**
- Multialleleic SNPs -- **gatk.snp.multiallelic.vcf.gz**
- Soft filtered final set of SNPs with iterative filteres specified -- **gatk.snp.soft_filtered.vcf.gz**
- Hard filtered final set of SNPs -- **gatk.snp.hard_filtered.vcf.gz**
### VCF statistics (found in 04_GATKvcfs)
- Quality statistics for a subset of 1M raw SNPs -- **gatk.snps.sub1M.qual.tsv**
- Genotype read depth for a subset of 1M raw SNPs -- **gatk.snps.sub1M.GT_DP.tsv**
- Iterative per-scaffold SNP filters -- **gatk.snp.filtered.FILTER_DIST**
- Per-scaffold summarized distribution of SNP filters -- **gatk.snp.filtered.FILTER_DIST.tsv**
- Summary statistics for all filtered sites (via `bcftools stats`) -- **gatk.snp.filtered.stats**
- Per-sample summary stats pulled from `bcftools stats` -- **gatk.snp.filtered.per-ind_stats.tsv**
- Variant allele frequency for heterozygous calls a subset of 20K filtered SNPs -- **gatk.snp.hard_filtered.thin20Kbp_hetVAF.tsv.gz**
### Quality control plots (found in Figures)
- Quality metric distributions with default cutoff values -- **raw_SNP_QC_metric_distributions.svg**
- Read depth distribution for a subset of genotypes with optional cutoff values -- **raw_SNP_GT_read_depth_distribution.svg**
- Filtered SNP per-sample quality metrics -- **filtered_SNP_QC_plots.svg**
- SNP filtering proportions per-scaffold -- **per-scaffold_SNP_filtering.svg**

## Output file structure (within specified home directory)
```
slurmout/        01_HTS_PreProc/        02_AlignDRAGMAP/        03_AlignStats/        04_GATKvcfs/          Figures/      temp/
  01_htstream      sampleID1              sampleID1               sampleID1             dragstr_model
  02_DRAGMAP       sampleID2              sampleID2               sampleID2             genomicDB
  03_merge         ...                    ...                     ...                   indiv_scaff_gvcfs
  04_stats                                                                              scaff_vcfs
  05_hapcall
  06_genotype
```

## Pipeline components:
### 1. clean_align_callSNPs.sbatch
Master script for executing the pipeline.
### 2. HTS_preproc.slurm
Use HTStream to clean raw paired FASTQ files by screening Illumina PhiX library sequences, trimming adapers, quality trimming read ends, removing 'N's, and filtering reads smaller than 100bp.
### 3. hashDRAGMAP.slurm
Build the reference genome hash table for DRAGMAP.
### 4. alignDRAGMAP.slurm
Use DRAGMAP to align cleaned reads, use picard to replace read group information, and use GATK MarkDuplicates to mark and remove PCR duplicate reads.
### 5. samtools_merge.slurm
If samples were sequenced across multiple lanes, use samtools to merge bams by sample. 
### 6. genome_wins.slurm
Use bedtools to create a bed file of the reference genome broken into 50Kb windows with 10% overlap.
### 7. align_stats.slurm
Use samtools, bedtools, and GATK CollectAlignmentSummaryMetrics to alculate alignment statistics for each sample including mean depth of coverage per sample (whole genome, per scaffold, and per 50Kb window), # and % reads aligned, read and insert length mean and distribution, etc.
### 8. STRtable.slurm
Use GATK ComposeSTRTableFile to create a short tandem repeat (STR) location table of the reference genome used for DragSTR model auto-calibration.
### 9. bam_to_gvcf.slurm
Use GATK to calibrate the DragSTR model (CalibrateDragstrModel), call individual variants (HaplotypeCaller), and compress individual GVCFs (ReblockGVCF).
### 10. gvcf_to_vcf_scaff.slurm
Use GATK to import single-sample GVCFs into per-scaffold databases (GenomicsDBImport) and joint call variants (GenotypeGVCFs).
### 11. scaff_vcf_to_snp_vcf.slurm
Use bcftools to combine per-scaffold VCFs and filter variants. Minor allele frequency filters are internally calculated as $`5/(total * 2)`$ to ensure at least 3 individuals possess a given allele ([Rochette & Catchen 2017](http://dx.doi.org/10.1038/nprot.2017.123)). Other filters include removal of INDELs, mutialleleic SNPs, invariable sites, per-genotype read depth above or below a given value, QUAL value < 10.4139, per-SNP missingness > 20%, and other [hard filtering recommendations](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants) from the Broad Institute. These filters can all be adjusted by editing this script. Quality control figures are also generated using the following sub-scripts:
- **raw_SNP_QC_metric_dists.R** - Visualize quality metric distributions for raw SNPs
- **raw_SNP_GT_read_depth_dist.R** - Visualize per-genotype read depth distribution
- **filtered_SNP_QC.R** - Visualize per-sample quality metrics and per-scaffold filtering proportions
### 12. job_report.slurm
Report the percent of successful and failed jobs for each pipeline component.
