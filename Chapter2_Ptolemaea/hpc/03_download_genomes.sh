# ======================================================================
# Thesis Chapter 2, Methods: genome retrieval
# Downloads the selected assemblies.
# Original location (Kelvin2 HPC): /users/40204129/Chp4_Ptolemaea/03.1_DownloadGenomes.sh
# ======================================================================
for f in accessions/*.txt; do
  tag=$(basename "$f" .txt)
  datasets download genome accession --inputfile "$f" \
      --include genome --filename "genomes/${tag}.zip"
  unzip -q -o "genomes/${tag}.zip" -d "genomes/${tag}"
done

# Build the accession table for Supplementary Table S1
{ echo -e "species\taccession";
  for f in accessions/*.txt; do
    tag=$(basename "$f" .txt)
    awk -v s="$tag" '{print s"\t"$0}' "$f"
  done
} > SupplementaryTableS1_accessions.tsv
