# ======================================================================
# Thesis Chapter 2, Methods: genome retrieval (Supplementary S1)
# Same as 02 for S. aureus only, completeness relaxed to >=98.85%.
# Original location (Kelvin2 HPC): /users/40204129/Chp4_Ptolemaea/02.2_DefineQuality.sh
# ======================================================================
for f in metadata/Staphylococcus_aureus.tsv; do
  tag=$(basename "$f" .tsv)
  tail -n +2 "$f" \
    | awk -F'\t' '($3>=98.85) && ($4<=2)' \
    | sort -t$'\t' -k7,7r -k5,5nr \
    | head -100 \
    | cut -f1 \
    > "accessions/${tag}.txt"
  echo "$tag: selected $(wc -l < accessions/${tag}.txt)"
done
