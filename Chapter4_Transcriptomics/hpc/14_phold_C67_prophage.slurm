#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Figure S4.3
# Phold on the C67 prophage Pharokka output.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/3_ProphageComparison/C67_prophage_pharokka_attempt2/RunPhold.sub
# ======================================================================

#SBATCH --job-name=Phold_C67_prophage
#SBATCH --error=Phold_C67_prophage.err
#SBATCH --output=Phold_C67_prophage.log
#SBATCH --partition=k2-sandbox,k2-bioinf,k2-hipri
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --mem=124G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL


#phold install -t 16
phold run -i pharokka.gbk  -o C67_prophage_phold -t 16 -p C67_prophage -f --cpu --ultra_sensitive

