# ======================================================================
# Thesis Chapter 2, Methods: genome retrieval (Supplementary S2)
# Pulls metadata for all complete RefSeq genomes of the 7 species (NCBI datasets).
# Original location (Kelvin2 HPC): /users/40204129/Chp4_Ptolemaea/01_PullMetadata.sh
# ======================================================================
SPECIES=("Escherichia coli" "Klebsiella pneumoniae" "Acinetobacter baumannii" \
         "Pseudomonas aeruginosa" "Enterococcus faecium" "Staphylococcus aureus" \
         "Enterobacter")

mkdir -p metadata accessions genomes

for sp in "${SPECIES[@]}"; do
  tag=$(echo "$sp" | tr ' ' '_')
  datasets summary genome taxon "$sp" \
      --assembly-source RefSeq \
      --assembly-level complete \
      --assembly-version latest \
      --as-json-lines \
  | dataformat tsv genome \
      --fields accession,organism-name,checkm-completeness,checkm-contamination,assmstats-contig-n50,assmstats-number-of-contigs,assminfo-release-date \
      > "metadata/${tag}.tsv"
  echo "$sp: $(( $(wc -l < metadata/${tag}.tsv) - 1 )) complete RefSeq genomes available"
done
