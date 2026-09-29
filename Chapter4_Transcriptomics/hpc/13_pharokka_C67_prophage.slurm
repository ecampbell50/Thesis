#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Figure S4.3
# Pharokka annotation of the C67 prophage.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/3_ProphageComparison/RunPharokka.sub
# ======================================================================

#SBATCH --job-name=Pharokka_C67_prophage
#SBATCH --error=Pharokka_C67_prophage.err
#SBATCH --output=Pharokka_C67_prophage.log
#SBATCH --partition=k2-bioinf,k2-hipri
#SBATCH --time=1:59:59
#SBATCH --nodes=1
#SBATCH --mem=124G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

pharokka.py -i C67_prophage.fna  -o C67_prophage_pharokka_attempt2 -d ~/.local/bin/../databases/ -t 16 -f

