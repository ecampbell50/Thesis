#!/bin/bash
# ======================================================================
# Thesis Chapter 2, Figure 2.1 and Results (Supplementary S4-S10)
# Runs Ptolemaea on the 100 Saureus genomes.
# Original location (Kelvin2 HPC): /users/40204129/Chp4_Ptolemaea/Ptolemaea_Saureus/Ptolemaea_Saureus.sub
# ======================================================================

#SBATCH --job-name=Saureus_Ptolemaea
#SBATCH --error=Saureus_Ptolemaea_round3.err
#SBATCH --output=Saureus_Ptolemaea_round3.log
#SBATCH --partition=k2-bioinf,k2-hipri
#SBATCH --time=2:59:59
#SBATCH --nodes=1
#SBATCH --mem=124G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

bash scripts/Ptolemaea.sh .

