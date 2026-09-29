# ======================================================================
# Thesis Chapter 5, Methods: prophage annotation
# Collects every 'Provirus' row from the geNomad summaries into one table.
# Original location (Kelvin2 HPC): /mnt/scratch2/users/40204129/Chapter3_Coevolution/04_genomad/GetProviruses.sh
# ======================================================================
for genomad_dir in *_genomad; do
    [[ "$genomad_dir" == *test* ]] && continue
    GENOME="${genomad_dir/_genomad/}"
    awk -F'\t' -v genome="$GENOME" \
        '$3 == "Provirus" {print genome "	" $0}' \
        "$genomad_dir"/*summary/*virus_summary.tsv
done
