#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Methods: RNA-seq processing
# FastQC on the raw RNA-seq reads.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/1_QC_Preprocessing/fastqc_EmmetTrans.slurm
# ======================================================================

#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --job-name=fastqc_EmmetTrans_20250919
#SBATCH --error=fastqc_EmmetTrans_20250919.error
#SBATCH --output=fastqc_EmmetTrans_20250919.txt
#SBATCH --partition=k2-medpri,k2-bioinf
#SBATCH --nodes=1
#SBATCH --time=05:00:00
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Activate the conda environment with FastQC installed
module load apps/fastqc/0.11.8/noarch

# The working directory to where your samples are
search_dir="/mnt/scratch2/users/40204129/3_Transcriptomics/0_AzentaRAWdata/00_fastq"

# Set the path to the directory where you want to store FastQC results
output_dir="/mnt/scratch2/users/40204129/3_Transcriptomics/1_QC_Preprocessing/1_FastQCOutput"

mkdir $output_dir

# Find all FASTQ files in subdirectories and loop through them

find $search_dir -type f -name "*R*_001.fastq.gz" | xargs fastqc -t 8 -o "$output_dir"

