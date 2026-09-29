#!/bin/bash
# ======================================================================
# Thesis Chapter 4, Table 4.1 (eggNOG-mapper), Figures 4.4, S4.7, S4.8
# eggNOG-mapper on the D68 Bakta proteins.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/1_StrainGenomes/D68_eggNOG/eggNOG_root_D68.sub
# ======================================================================
#SBATCH --job-name=D68_root_eggnog
#SBATCH --cpus-per-task=16
#SBATCH --output=D68_root_eggnog.%j.out
#SBATCH --error=D68_root_eggnog.%j.err
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --partition=k2-hipri,k2-bioinf

#module load apps/anaconda3
#source activate /mnt/scratch/igfs-anaconda/conda-envs/eggnog

# Advanced emapper run with custom parameters
emapper.py \
    --cpu 16 \
    -i /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/1_StrainGenomes/D68_bakta/D68.faa \
    -o D68_root_eggNOG \
    --seed_ortholog_evalue 0.001 \
    --seed_ortholog_score 60 \
    --report_orthologs

