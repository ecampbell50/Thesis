#!/bin/bash
# ======================================================================
# Thesis Chapter 2, Figure 2.1 and Results (Supplementary S4-S10)
# Runs Ptolemaea on the 100 Paeruginosa genomes.
# Original location (Kelvin2 HPC): /users/40204129/Chp4_Ptolemaea/Ptolemaea_Paeruginosa/Ptolemaea_Paeruginosa.sub
# ======================================================================

#SBATCH --job-name=Paeruginosa_Ptolemaea
#SBATCH --error=Paeruginosa_Ptolemaea_round4.err
#SBATCH --output=Paeruginosa_Ptolemaea_round4.log
#SBATCH --partition=k2-bioinf,k2-hipri,k2-sandbox
#SBATCH --time=00:29:59
#SBATCH --nodes=1
#SBATCH --mem=124G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

bash scripts/Ptolemaea.sh .

