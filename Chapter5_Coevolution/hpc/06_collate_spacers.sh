# ======================================================================
# Thesis Chapter 5, Methods: CRISPR spacers
# Collates per-genome spacers into all_spacers.fasta.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/Chapter3_Coevolution/09_minced/collate_arrays.sh
# ======================================================================
cd /mnt/scratch2/users/40204129/Chapter3_Coevolution/09_minced
: > all_spacers.fasta
for sp in arrays/*_spacers.fa; do
  [[ -s "$sp" ]] || continue                 # skip empties (genomes with no arrays)
  gid=$(basename "$sp" _spacers.fa)          # e.g. 1005042.3
  sed "s/^>/>${gid}|/" "$sp" >> all_spacers.fasta
done
# sanity check:
echo "spacers: $(grep -c '^>' all_spacers.fasta) from $(grep '^>' all_spacers.fasta | cut -d'|' -f1 | sort -u | wc -l) genomes"
head -4 all_spacers.fasta
