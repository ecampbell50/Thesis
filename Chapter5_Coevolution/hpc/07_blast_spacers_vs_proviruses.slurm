#!/bin/bash
# ======================================================================
# Thesis Chapter 5, Figures 5.10-5.13
# BLASTn of all spacers against all proviruses (protospacer hits).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/Chapter3_Coevolution/spacerProphageBlast.sub
# ======================================================================
#SBATCH --job-name=spacer_blast
#SBATCH --partition=k2-bioinf,k2-lowpri
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --output=spacer_blast.out
#SBATCH --error=spacer_blast.err
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=END,FAIL

cd /mnt/scratch2/users/40204129/Chapter3_Coevolution/09_minced
BASE=/mnt/scratch2/users/40204129/Chapter3_Coevolution
# provirus db, with each provirus's FILENAME as its blast id (so hits map to prophage_id)
: > all_proviruses.fna
for f in "$BASE"/05_proviruses/*.fna; do
  id=$(basename "$f")
  awk -v id="$id" '/^>/{print ">"id; next}{print}' "$f" >> all_proviruses.fna
done
makeblastdb -in all_proviruses.fna -dbtype nucl -out provirus_db >/dev/null

# match spacers -> proviruses, then keep only protospacer-grade hits
blastn -task blastn-short -query all_spacers.fasta -db provirus_db -dust no -evalue 1 -num_threads 8 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore' \
  > spacer_vs_provirus.raw.tsv
awk -F'\t' '($3>=95)&&($5<=1)&&($4/$7>=0.90)' spacer_vs_provirus.raw.tsv > spacer_vs_provirus.tsv
echo "protospacer-grade hits: $(wc -l < spacer_vs_provirus.tsv)"
