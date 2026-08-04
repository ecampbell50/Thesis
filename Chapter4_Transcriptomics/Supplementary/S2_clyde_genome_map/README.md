# Supplementary S2 — Clyde genome map

**Script:** `script/phage_genome_map.R`  (generates **both** Bonnie and Clyde maps)
**Figure:** `figure/clyde_genome_map.pdf`
**Data:** `data/Clyde_allannos.gff3`  (genome length 34,734 bp)

Circular genome map: forward-strand genes on the outer ring, reverse on the
inner ring, coloured by functional category inferred from the GFF3 `product`
annotation.

> Note: same script as S1. The original reads the GFF from
> `~/Desktop/files_from_hpc/files_to_local/Clyde_allannos.gff3`; a copy is in
> `data/`. (A standalone older version, `clyde_genome_map.R`, also exists in the
> working tree but `phage_genome_map.R` is the current one.)
