# ======================================================================
# Thesis Chapter 2, Methods: genome retrieval (Supplementary S1)
# CheckM filter (>=99% complete, <=2% contamination) and picks 100 genomes per species.
# Original location (Kelvin2 HPC): /users/40204129/Chp4_Ptolemaea/02.1_DefineQuality.sh
# ======================================================================
for f in metadata/*.tsv; do
  tag=$(basename "$f" .tsv)
  tail -n +2 "$f" \
    | awk -F'\t' '($3>=99) && ($4<=2)' \
    | sort -t$'\t' -k7,7r -k5,5nr \
    | head -100 \
    | cut -f1 \
    > "accessions/${tag}.txt"
  echo "$tag: selected $(wc -l < accessions/${tag}.txt)"
done
