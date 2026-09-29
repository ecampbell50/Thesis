#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Table 4.1 (Phold), Figure S4.1
# Phold structural annotation of phage Bonnie.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/2_PhageGenomes/2_Phold/PHOLD_BONNIE.sub
# ======================================================================

#SBATCH --job-name=PHOLDcpu_Bonnie_Test
#SBATCH --error=PHOLDcpu_Bonnie_Test.err
#SBATCH --output=PHOLDcpu_Bonnie_Test.log
#SBATCH --partition=k2-bioinf
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --mem=124G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

phold run \
	-i /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/2_PhageGenomes/1_Pharokka/BONNIE_pharokka/BONNIE.gbk \
	-o /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/2_PhageGenomes/2_Phold/BONNIE_phold_cpu \
	-d /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/2_PhageGenomes/phold_database_gpu \
	-t 16 \
	-f \
	--cpu \

