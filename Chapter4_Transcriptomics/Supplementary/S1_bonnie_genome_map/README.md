# Supplementary S1 — Bonnie genome map

**Script:** `script/phage_genome_map.R`  (generates **both** Bonnie and Clyde maps)
**Figure:** `figure/bonnie_genome_map.pdf`
**Data:** `data/Bonnie_allannos.gff3`

Circular genome map: forward-strand genes on the outer ring, reverse on the
inner ring, drawn as arrows and coloured by functional category inferred from
the GFF3 `product` annotation. Genome length auto-read from the
`##sequence-region` header.

> Note: the original script reads the GFF from
> `~/Desktop/files_from_hpc/files_to_local/Bonnie_allannos.gff3`. A copy is
> included here in `data/`; update the path in the script if re-running from
> this folder.
