# Chapter 5 — Prophage–host co-evolution

Run from this folder. `scripts/prep_*.py` and `scripts/analyse_prophage_defence.py` build the tables in `data/` that the rest read.

| Thesis | Code |
|---|---|
| Fig 5.4, 5.5 | Netpass: [github.com/ecampbell50/Netpass](https://github.com/ecampbell50/Netpass) |
| Fig 5.6 | `scripts_revisions/select_medoids.py`, then the VIRIDIC web server (`viridic/`) |
| Fig 5.7 | `scripts_revisions/blast_pc74_vs_ice99.sh`, `plot_ice99_genome_map.py` |
| Fig 5.8 | `scripts/plot_bipartite_network.R` |
| Fig 5.9 | `scripts/make_itol_community_datasets.R`, drawn in iTOL (`data/itol/`) |
| Fig 5.10 | `CrisprSpacer_Analysis.ipynb` |
| Fig 5.11–5.13 | `scripts/make_itol_crispr_linkages.py`, drawn in iTOL |
| Fig 5.14, Fig S5.4 | `scripts/plot_community_defence_split_heatmap.py` |
| Fig 5.15 | `scripts/analyse_modularity.R`, `scripts/analyse_nestedness.R` |
| Fig 5.16, S5.3, S5.6, S5.7 | `BipartiteNetwork_Analysis.ipynb` (tables from `scripts/analyse_confounders.py`) |
| Fig S5.1 | `scripts_revisions/plot_pc74_vs_ice139.py` |
| Fig S5.2 | sourmash commands in `pipeline/Methodology.md` |
| Fig S5.5 | `scripts/crosswalk_old_vs_new_communities.R`, drawn in iTOL |
| CRISPR targeting stats | `scripts/analyse_crispr_targeting.py` |
| Patristic distances | `scripts/community_phylo_distance.R`, `scripts_revisions/procom4_treedist.py` |

`hpc/` has the Kelvin2 jobs: geNomad, MinCED and the spacer BLAST.
