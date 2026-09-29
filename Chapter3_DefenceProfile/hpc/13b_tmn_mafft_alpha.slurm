#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.7 and Figure S3.3
# MAFFT alignment of the Tmn-alpha family (GeneRax input).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/Tmn_Gene_Work/MAFFT_Tmn_alpha.sub
# ======================================================================
#SBATCH --job-name=MAFFT_Tmn_alpha
#SBATCH --error=MAFFT_Tmn_alpha.err
#SBATCH --output=MAFFT_Tmn_alpha.out
#SBATCH --partition=k2-sandbox,k2-hipri,k2-bioinf
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --ntasks=4
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

INPUT_FILE="Tmn_alpha_genes.faa"

mafft --maxiterate 1000 --localpair --thread 4 "$INPUT_FILE" > "Tmn_alpha_MAFFT_aligned.faa"

