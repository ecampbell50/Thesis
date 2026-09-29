#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Figures S4.1, S4.2
# eggNOG-mapper on the Bonnie and Clyde proteins.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/2_PhageGenomes/5_eggNOG/eggNOG_phages.sub
# ======================================================================
#SBATCH --job-name=Clyde_eggnog
#SBATCH --cpus-per-task=16
#SBATCH --output=Clyde_eggnog.%j.out
#SBATCH --error=Clyde_eggnog.%j.err
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=100G
#SBATCH --time=00:45:00
#SBATCH --partition=k2-hipri,k2-bioinf

#module load apps/anaconda3
#source activate /mnt/scratch2/igfs-anaconda/conda-envs/eggnog

# Advanced emapper run with custom parameters
emapper.py \
    --cpu 16 \
    -i /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/2_PhageGenomes/1_Pharokka/CLYDE_pharokka/CLYDE.faa \
    -o Clyde_eggNOG_caudovirales \
    --seed_ortholog_evalue 0.001 \
    --seed_ortholog_score 60 \
    --tax_scope 28883 \
    --report_orthologs

