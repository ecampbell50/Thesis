# Chapter 5 revision analyses (2026-07-30)

Scripts written while addressing supervisor comments on the coevolution chapter.

| Script | Purpose | Addresses |
|---|---|---|
| `select_medoids.py` | Pick the medoid representative of each of the 15 largest prophage communities (highest mean within-community Sourmash Jaccard). Produces the VIRIDIC input list. | comment 346 |
| `blast_pc74_vs_ice99.sh` | Test whether the PRO_com_74 representative is contained within the ICE_community_99 representative, following VIRIDIC's 50.1% similarity score. Result: 85.1% of PRO_com_74 (40,273 / 47,313 bp) aligns at 89-95% identity across ~44 kb of the 105 kb ICE, in reverse orientation. | comments 144, 191, F12 |
| `procom4_treedist.py` | Patristic distance summaries: clonality of bacterial community 2, species-wide scale, and PRO_com_4 carriage per bacterial community. Needs `treedist_all` output. | comments 211-214 |

VIRIDIC itself was run on the web service (<https://rhea.icbm.uni-oldenburg.de/VIRIDIC/>) using
`../Top15Medoids_VIRIDIC_input.fasta`; outputs are `../VIRIDIC_{heatmap.pdf,sim-dist_table.tsv,cluster_table.tsv}`.
The standalone tarball is no longer distributed and there is no Docker image.

`treedist` (Creevey, doi:10.5281/zenodo.1244019) is compiled in `../treedist/`.
