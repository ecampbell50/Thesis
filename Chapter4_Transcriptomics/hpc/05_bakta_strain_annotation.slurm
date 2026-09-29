# ======================================================================
# Thesis Chapter 4, Table 4.1 (Bakta)
# Bakta annotation of C67, D32 and D68 (switch GENOME/STRAIN lines per strain).
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/3_Transcriptomics/7_GenomeAnnotations/1_StrainGenomes/BAKTA_tester.sh
# ======================================================================
export MPLCONFIGDIR=/mnt/scratch2/users/40204129/tmp_matplotlib
mkdir -p $MPLCONFIGDIR
export PYTHONNOUSERSITE=1

#GENOME="C67_HBHC_prefixed.fasta"
#GENOME="D32_LBLC_prefixed.fasta"
GENOME="D68_LBHC_prefixed.fasta"

STRAIN="D68"

bakta \
	--db "/mnt/scratch2/igfs-databases/Bakta/bakta_1.11.0_db_6.0/db" \
	-p "${STRAIN}" \
	--output "${STRAIN}_bakta" \
	--force \
	--genus Streptococcus \
	--species suis \
	--gram + \
	--locus-tag "${STRAIN}_prot" \
	--verbose \
	--threads 8 \
	"$GENOME"

