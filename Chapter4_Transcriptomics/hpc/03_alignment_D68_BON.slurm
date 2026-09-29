#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Methods: RNA-seq processing (Figures 4.2-4.4)
# PhageExpressionAtlas (Nextflow: Cutadapt, HISAT2, featureCounts) for D68 x BON.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/4_SuccessfulRuns/D68_BON_20251222/D68_BON_20251222.sub
# ======================================================================
#SBATCH --job-name=D68_BON_20251222
#SBATCH --error=D68_BON_20251222_%j.err
#SBATCH --output=D68_BON_20251222_%j.out
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


# Initialize conda (exactly what your .bashrc does)
eval "$('/opt/apps/mamba/1.3.1/python3-3.10.5+gcc-9.3.0/mambaforge/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
export PATH="/opt/apps/mamba/1.3.1/python3-3.10.5+gcc-9.3.0/mambaforge/bin:$PATH"
export CONDA_PKGS_DIRS=/mnt/scratch2/users/40204129/conda/pkgs
export CONDA_ENVS_DIRS=/mnt/scratch2/users/40204129/conda/envs

WORKING_DIR="/mnt/scratch2/users/40204129/3_Transcriptomics/4_SuccessfulRuns/D68_BON_20251222"
cd $WORKING_DIR

nextflow run pipeline \
	--reads '/mnt/scratch2/users/40204129/3_Transcriptomics/0_AzentaRAWdata/00_fastq/*-D68-BON-*_{R1,R2}_001.fastq.gz' \
	--hostGenome "${WORKING_DIR}/genomes_and_annotations/D68.fna" \
	--phageGenome "${WORKING_DIR}/genomes_and_annotations/BONNIE.fasta" \
	--hostGFF "${WORKING_DIR}/genomes_and_annotations/D68_annotated_noseq.gff3" \
	--phageGFF "${WORKING_DIR}/genomes_and_annotations/BONNIE_phold.gff3" \
	--outputDir "D68_BON_20251222" \
	--countFeature "CDS" \
	--featureIdentifier "ID" \
	--conda_path "/mnt/scratch2/users/40204129/conda/envs" \
	--pairedEnd=true \
	-resume

