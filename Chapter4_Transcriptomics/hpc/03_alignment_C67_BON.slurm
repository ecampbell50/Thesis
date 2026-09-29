#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Methods: RNA-seq processing (Figures 4.2-4.4)
# PhageExpressionAtlas (Nextflow: Cutadapt, HISAT2, featureCounts) for C67 x BON.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/4_SuccessfulRuns/C67_BON_20251219/C67_BONNIE_20251219.sub
# ======================================================================
#SBATCH --job-name=C67_BON_20251219
#SBATCH --error=C67_BON_20251219_%j.err
#SBATCH --output=C67_BON_20251219_%j.out
#SBATCH --partition=k2-himem,k2-bioinf,k2-medpri
#SBATCH --time=07:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Load modules
module purge
module load nextflow/24.07.27/java-22.0.2

# Set JAVA_HOME (critical!)
export JAVA_HOME=/opt/apps/jdk/22.0.2

nextflow run pipeline \
	--reads '/mnt/scratch2/users/40204129/3_Transcriptomics/0_AzentaRAWdata/00_fastq/*-C67-BON-*_{R1,R2}_001.fastq.gz' \
	--hostGenome "/mnt/scratch2/users/40204129/3_Transcriptomics/8_AnnotatedRun/genomes_and_annotations/C67.fna" \
	--phageGenome "/mnt/scratch2/users/40204129/3_Transcriptomics/8_AnnotatedRun/genomes_and_annotations/BONNIE.fasta" \
	--hostGFF "/mnt/scratch2/users/40204129/3_Transcriptomics/8_AnnotatedRun/genomes_and_annotations/C67_annotated_noseq.gff3" \
	--phageGFF "/mnt/scratch2/users/40204129/3_Transcriptomics/8_AnnotatedRun/genomes_and_annotations/BONNIE_phold.gff3" \
	--outputDir "C67_BON_20251219" \
	--countFeature "CDS" \
	--featureIdentifier "ID" \
	--conda_path "/mnt/scratch2/users/40204129/conda/envs" \
	--pairedEnd=true \
	-resume

