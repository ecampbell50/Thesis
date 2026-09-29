#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.7
# IQ-TREE gene tree of the Tmn alignment (-boot 1000).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/Tmn_Gene_Work/IQTree_Tmn.sub
# ======================================================================
#SBATCH --job-name=IQtree_Tmn
#SBATCH --error=IQtree_Tmn.err
#SBATCH --output=IQtree_Tmn.out
#SBATCH --partition=k2-hipri
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=8
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Run IQ-TREE on single gene alignment
iqtree2 -s Tmn_MAFFT_aligned.faa \
        -m MFP \
        -B 1000 \
        -alrt 1000 \
        -nt 8

echo "Tree building complete!"
