#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.6
# FastTree (GTR, 1000 bootstraps) on the Roary core-gene alignment.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/output/06_pangenome/Roary_Phylogeny_wSpenue_21Jan2026/Ssuis_PhyloTree_20260126.sub
# ======================================================================

#SBATCH --job-name=PhyloTree_FastTree
#SBATCH --error=PhyloTree_FastTree.err
#SBATCH --output=PhyloTree_FastTree.log
#SBATCH --partition=k2-bioinf,k2-medpri
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

FastTree -nt -gtr -boot 1000 core_gene_alignment.aln > SsuisPhylo_SpneuOutgroup_119cg.nwk

