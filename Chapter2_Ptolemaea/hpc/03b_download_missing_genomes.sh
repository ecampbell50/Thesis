# ======================================================================
# Thesis Chapter 2, Methods: genome retrieval
# Re-downloads the few assemblies that failed in 03.
# Original location (Kelvin2 HPC): /users/40204129/Chp4_Ptolemaea/03.2_DownloadMissingGenomes.sh
# ======================================================================
# Efaecium — 2 missing
printf 'GCF_907165365.1\nGCF_965641755.1\n' > /tmp/missing_Efaecium.txt
datasets download genome accession \
    --inputfile /tmp/missing_Efaecium.txt \
    --include genome \
    --filename /tmp/Efaecium_missing.zip
unzip -q -o /tmp/Efaecium_missing.zip -d /tmp/Efaecium_missing
find /tmp/Efaecium_missing/ncbi_dataset/data -name "*.fna" \
    -exec cp {} Ptolemaea_Efaecium/genomes/ \;

# Kpneumoniae — 1 missing
printf 'GCF_978018475.1\n' > /tmp/missing_Kpneumoniae.txt
datasets download genome accession \
    --inputfile /tmp/missing_Kpneumoniae.txt \
    --include genome \
    --filename /tmp/Kpneumoniae_missing.zip
unzip -q -o /tmp/Kpneumoniae_missing.zip -d /tmp/Kpneumoniae_missing
find /tmp/Kpneumoniae_missing/ncbi_dataset/data -name "*.fna" \
    -exec cp {} Ptolemaea_Kpneumoniae/genomes/ \;

# Saureus — 3 missing
printf 'GCF_900324225.1\nGCF_900324365.1\nGCF_900324385.1\n' > /tmp/missing_Saureus.txt
datasets download genome accession \
    --inputfile /tmp/missing_Saureus.txt \
    --include genome \
    --filename /tmp/Saureus_missing.zip
unzip -q -o /tmp/Saureus_missing.zip -d /tmp/Saureus_missing
find /tmp/Saureus_missing/ncbi_dataset/data -name "*.fna" \
    -exec cp {} Ptolemaea_Saureus/genomes/ \;
