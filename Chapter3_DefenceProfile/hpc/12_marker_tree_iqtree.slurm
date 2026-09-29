#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.7 (marker-gene species tree)
# IQ-TREE on the catsequences supermatrix (species tree for GeneRax).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/Phylogeny_alignments/IQTree.sub
# ======================================================================
#SBATCH --job-name=IQtree
#SBATCH --error=IQtree.err
#SBATCH --output=IQtree.out
#SBATCH --partition=k2-medpri,k2-bioinf
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=8
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/iqtree

# Run IQ-TREE with partitions and ultrafast bootstrap
iqtree2 -s allseqs.fas \
        -p allseqs.partitions.nex \
        -m MFP \
        -B 1000 \
        -alrt 1000 \
        -nt 8

# -s: alignment file
# -p: partition file (each COG as separate partition)
# -m MFP: ModelFinder Plus (auto-selects best model per partition)
# -B 1000: 1000 ultrafast bootstrap replicates
# -alrt 1000: 1000 SH-aLRT replicates (additional support)
# -nt 8: use 8 threads
