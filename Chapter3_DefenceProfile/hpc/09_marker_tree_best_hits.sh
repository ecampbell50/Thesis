# ======================================================================
# Thesis Chapter 3, Figure 3.7 (marker-gene species tree)
# Keeps the best hit per genome for each COG.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/11_DefenseProfileofSsuis/1_SsuisGenomes/ssuis_defence_project/GeneReconcilliation_of_tmn/COG_to_Seq_BLASTS/GetTopHits.sh
# ======================================================================
for i in COG*_blast.txt;
do
	COG="${i/_blast.txt/}"
	cat "${COG}_blast.txt" | tr '\t' ' ' | sed 's/@/ /' | sort -t' ' -k2,2 -k13,13rn | awk 'BEGIN{PREV=""};{if(PREV!=$2 && $13>60) {PREV=$2; print($2"@"$3)}}' > "${COG}_bestgenomehits.txt"
done
