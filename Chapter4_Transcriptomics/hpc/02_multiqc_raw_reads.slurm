#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Methods: RNA-seq processing
# MultiQC summary of the FastQC reports.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/1_QC_Preprocessing/multiqc_EmmetTrans_20250919.slurm
# ======================================================================

#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=multiqc_EmmetTrans_20250919
#SBATCH --error=multiqc_EmmetTrans_20250919.error
#SBATCH --output=multiqc_EmmetTrans_20250919.txt
#SBATCH --partition=k2-hipri,k2-bioinf
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mail-user=ecampbell50@qub.ac.uk@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load multiqc/1.2
module load Python-Packages/multiqc/1.23

# Change the working directory to where your FastQC output files are
cd /mnt/scratch2/users/40204129/3_Transcriptomics/1_QC_Preprocessing/1_FastQCOutput

# Set the path to the directory where you want to store FastQC results
output_dir="/mnt/scratch2/users/40204129/3_Transcriptomics/1_QC_Preprocessing/1_FastQCOutput"

# Find all FASTQ files in subdirectories and loop through them
multiqc --outdir $output_dir .

