# ======================================================================
# Thesis Chapter 3, Figure 3.7 (marker-gene species tree)
# Prefixes protein IDs with genome IDs before the marker-gene BLAST.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/prokka_faas/add_genomeID.sh
# ======================================================================
for i in *genomic.faa
do
	ID="$(echo $i | cut -d'_' -f1,2)"
	sed "s/^>/>$ID@/g" $i > "${ID}_editedheaders.faa"
done


