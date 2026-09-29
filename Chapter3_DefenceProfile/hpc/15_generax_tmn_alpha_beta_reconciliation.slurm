#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.7, S3.3, S3.4
# GeneRax reconciliation of Tmn-alpha and Tmn-beta against the species tree.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/Phylogeny_alignments/generax_split.slurm
# ======================================================================
#SBATCH --job-name=GeneRax_Tmn_split
#SBATCH --error=generax_split.err
#SBATCH --output=generax_split.out
#SBATCH --partition=k2-hipri,k2-bioinf
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# cd /users/40204129/sharedscratch/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/Phylogeny_alignments

# rm -rf generax_Tmn_split_output

./GeneRax/build/bin/generax \
    -f Tmn_split_families.txt \
    -s rooted_tree.nwk \
    -r UndatedDTL \
    --prefix generax_Tmn_split_output

