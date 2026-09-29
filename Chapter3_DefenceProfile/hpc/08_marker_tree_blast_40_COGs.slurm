#!/bin/bash
# ======================================================================
# Thesis Chapter 3, Figure 3.7 (marker-gene species tree)
# BLASTs the 40 single-copy marker genes against each genome.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/COG_to_Seq_BLASTs.sub
# ======================================================================
#SBATCH --job-name=COG_blast
#SBATCH --error=COG_blast.err
#SBATCH --output=COG_blast.out
#SBATCH --partition=k2-hipri,k2-bioinf
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=8
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

module load apps/ncbiblast

for i in COGs/*.faa; do
    COG_NAME=$(basename "$i" .faa)
    echo "Processing ${COG_NAME}..."
    blastp -query "$i" \
           -db AllTmnContainingGenomes_andoutgroups.faa \
           -out "COG_to_Seq_BLASTS/${COG_NAME}_blast.txt" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
           -evalue 1e-5 \
           -max_target_seqs 500 \
           -num_threads 8

    echo "${COG_NAME} completed..."
done

echo "BLAST complete for all COGs"
