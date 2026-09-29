#!/bin/bash
# ======================================================================
# Thesis Chapter 5, Methods: CRISPR spacers (Figures 5.2, 5.10-5.13)
# MinCED CRISPR array and spacer calling on all genomes.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/Chapter3_Coevolution/minced_17jun26.sub
# ======================================================================
#SBATCH --job-name=crispr_spacers
#SBATCH --partition=k2-bioinf,k2-medpri
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=16G
#SBATCH --output=minced_17jun26.out
#SBATCH --error=minced_17jun26.err
#SBATCH --mail-user=ecampbell50@qub.ac.uk
#SBATCH --mail-type=END,FAIL

cd /mnt/scratch2/users/40204129/Chapter3_Coevolution/09_minced
mkdir -p arrays
: > all_spacers.fasta
for f in ../01_genomes/*.fna; do
  gid=$(basename "$f" .fna)
  minced -gffFull -spacers "$f" "arrays/${gid}.crisprs" "arrays/${gid}.gff" >/dev/null 2>&1
  sp="${gid}_spacers.fa"                                   # minced drops this in the cwd
  [[ -s "$sp" ]] || sp=$(ls "arrays/${gid}".*spacers*.fa 2>/dev/null | head -1)  # fallback
  [[ -s "$sp" ]] && sed "s/^>/>${gid}|/" "$sp" >> all_spacers.fasta
done
echo "spacers: $(grep -c '^>' all_spacers.fasta) from $(grep '^>' all_spacers.fasta | cut -d'|' -f1 | sort -u | wc -l) genomes"
