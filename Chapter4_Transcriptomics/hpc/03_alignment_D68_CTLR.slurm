#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Methods: RNA-seq processing (Figures 4.2-4.4)
# PhageExpressionAtlas (Nextflow: Cutadapt, HISAT2, featureCounts) for D68 x CTLR.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/4_SuccessfulRuns/D68_CTLR_20260126/D68_CTLR_20260126.sub
# ======================================================================
#SBATCH --job-name=D68_CTLR_20260126
#SBATCH --error=D68_CTLR_20260126_%j.err
#SBATCH --output=D68_CTLR_20260126_%j.out
#SBATCH --partition=k2-bioinf,k2-medpri
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

###
# Originally did 7 hours but timed out
# Added an extra 3
###


# Load modules
module purge
module load jdk/22.0.2
module load nextflow/24.07.27/java-22.0.2

# Initialize conda FIRST
eval "$('/opt/apps/mamba/1.3.1/python3-3.10.5+gcc-9.3.0/mambaforge/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
export PATH="/opt/apps/mamba/1.3.1/python3-3.10.5+gcc-9.3.0/mambaforge/bin:$PATH"
export CONDA_PKGS_DIRS=/mnt/scratch2/users/40204129/conda/pkgs
export CONDA_ENVS_DIRS=/mnt/scratch2/users/40204129/conda/envs

# NOW set JAVA_HOME (AFTER conda, so it doesn't get overridden)
export JAVA_HOME=/opt/apps/jdk/22.0.2
export PATH=$JAVA_HOME/bin:$PATH

# Verify
echo "=== ENVIRONMENT DEBUG ==="
echo "JAVA_HOME: $JAVA_HOME"
echo "Java location: $(which java)"
echo "Java version: $(java -version 2>&1 | head -1)"
echo "========================="

# Run nextflow
nextflow run pipeline \
	--reads '/mnt/scratch2/users/40204129/3_Transcriptomics/0_AzentaRAWdata/00_fastq/*-D68-CTLR-*_{R1,R2}_001.fastq.gz' \
	--hostGenome "/mnt/scratch2/users/40204129/3_Transcriptomics/8_AnnotatedRun/genomes_and_annotations/D68.fna" \
	--phageGenome "/mnt/scratch2/users/40204129/3_Transcriptomics/8_AnnotatedRun/genomes_and_annotations/BOTHPHAGES.fasta" \
	--hostGFF "/mnt/scratch2/users/40204129/3_Transcriptomics/8_AnnotatedRun/genomes_and_annotations/D68_annotated_noseq.gff3" \
	--phageGFF "/mnt/scratch2/users/40204129/3_Transcriptomics/8_AnnotatedRun/genomes_and_annotations/BonnieClyde_dual.gff3" \
	--outputDir "D68_CTLR_20260126" \
	--countFeature "CDS" \
	--featureIdentifier "ID" \
	--conda_path "/mnt/scratch2/users/40204129/conda/envs" \
	--pairedEnd=true \
	-resume
