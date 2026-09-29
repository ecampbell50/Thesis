#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.6
# Roary pangenome of 2,119 S. suis + 5 S. pneumoniae, core-gene alignment (-e).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/output/06_pangenome/PhyloPangenome_wSpenu_21Jan2026.sub
# ======================================================================

#SBATCH --job-name=Roary_phylogeny
#SBATCH --error=Roary_phylogeny_21Jan2026.err
#SBATCH --output=Roary_phylogeny_21Jan2026.log
#SBATCH --partition=k2-bioinf,k2-lowpri
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=30
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Run Roary
roary -p 30 \
      -e \
      -v \
      -f Roary_Phylogeny_wSpenue_21Jan2026 \
      -g 70000 \
      all_gffs/*.gff

