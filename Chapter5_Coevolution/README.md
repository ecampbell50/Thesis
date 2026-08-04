# Chapter 5 — Prophage–host co-evolutionary dynamics in *Streptococcus suis*

```
Chapter5_Coevolution/
├── data/                          all inputs (16 files + community_questions/, confounders/)
├── BipartiteNetwork_Analysis.ipynb  Python: rarefaction, promiscuity, confounders, length bimodality
├── CrisprSpacer_Analysis.ipynb      Python: CRISPR spacer → prophage network
├── scripts/                       10 scripts, one per analysis in the chapter
├── scripts_revisions/             ICE maps, medoid selection, patristic distances
├── pipeline/                      Methodology.md, crispr_spacer_pipeline.sh
├── figures/  tables/              generated outputs (figures gitignored)
└── README.md
```

### Scope

**Only code used by Chapter 5 is included.** The source project contains a number of
exploratory analyses that were not carried into the chapter — a Pagel test,
serotype–defence and prophage–defence associations, defence-burden confounders,
defence location, Scoary module analysis, coexistence and tripartite networks. All were
checked against the chapter text (0 mentions each) and **deliberately excluded**. If you
are looking for them, they remain in
`Desktop/Thesis/1_CoevolutionofSsuisProphages/`.

| Script | Chapter section |
|---|---|
| `prep_coevolution_inputs.py`, `prep_community_question_data.py` | builds the incidence matrices and per-question tables |
| `coevolution_common.R` | single source of truth for the BAC × PRO community matrix |
| `plot_bipartite_network.R` | Fig — community-level infection network |
| `analyse_modularity.R`, `analyse_nestedness.R` | Methods §Network modularity and nestedness |
| `plot_community_defence_split_heatmap.py` | Methods §Prophage community – defence system associations |
| `analyse_crispr_targeting.py` | Methods §CRISPR targeting analysis |
| `community_phylo_distance.R` | Methods §Core-genome phylogeny (patristic distances) |
| `make_itol_community_datasets.R` | iToL annotation files for the phylogeny figure |

## Figure → code map

**✅ reproducible · ⚠️ partial · 🌐 external service · ❌ no code found**

| Figure | File | Produced by | |
|---|---|---|---|
| Fig — PBIN schematic | `PBIN_infographic.png` | author-made, from published sources | — ¹ |
| Fig — CRISPR overview | `CRISPR_overview.png` | author-made, from published sources | — ¹ |
| Methods — weighted overlap | `WeightedOverlap_Schematic.png` | author-made, adapted from Edge_overlap Fig 1 | — ¹ |
| Fig — Netpass thresholds (4 panels) | `communities_vs_weight_{BAC,PRO}`, `elbow_detection_{49_BAC,44_PRO}` | **`github.com/ecampbell50/Netpass`** (`Main.py`) | ✅ |
| Fig — VIRIDIC heatmap | `VIRIDIC_heatmap.png` | VIRIDIC **web service** ² | 🌐 |
| Fig — ICE99/PC74 map | `ICE99_PC74_genome_map.pdf` | `scripts_revisions/plot_ice99_genome_map.py` | ✅ |
| Fig — bipartite network | `bipartite_full_network.png` | `scripts/plot_bipartite_network.R` | ✅ |
| Fig — phylogeny + communities | `Phylo_BC_PC.png` | iToL ³ | ✅ |
| Fig — spacer/prophage network | `CRISPR_bipartite.png` | `CrisprSpacer_Analysis.ipynb` | ✅ |
| Fig — BC3 CRISPR linkages | `BC3_CRISPR_linkages.png` | `scripts/make_itol_crispr_linkages.py 3` → iToL ⁴ | ✅ |
| Fig — BC5 CRISPR linkages | `BC5_CRISPR_linkages.png` | `scripts/make_itol_crispr_linkages.py 5` → iToL ⁴ | ✅ |
| Fig — pruned BC2–5 linkages | `pruned_BC2-5_linkages.png` | `scripts/make_itol_crispr_linkages.py 2 5` → iToL, pruned tree ⁴ | ✅ |
| Fig — defence split heatmap | `community_defence_split_heatmap_10SIG.png` | `scripts/plot_community_defence_split_heatmap.py` | ✅ |
| Fig — modularity & nestedness | `Q_NODF_min5.png` | `scripts/analyse_modularity.R` + `analyse_nestedness.R` ⁵ | ✅ |
| Fig — rarefaction | `Rarefraction.png` | `BipartiteNetwork_Analysis.ipynb` ⁶ | ✅ |
| S — ICE139/PC74 map | `ICE139_PC74_genome_map.pdf` | `scripts_revisions/plot_pc74_vs_ice139.py` | ✅ |
| S — sourmash matrices | `BAC_k21s1000.matrix.pdf`, `PRO_k15s500.matrix.pdf` | `sourmash plot` ⁷ | ✅ |
| S — promiscuity vs size | `Q1_promiscuity_vs_size.png` | `BipartiteNetwork_Analysis.ipynb` | ✅ |
| S — full binomial heatmap | `community_defence_split_heatmap_ALL.png` | `scripts/plot_community_defence_split_heatmap.py` | ✅ |
| S — BAC recalculation | `BC_recalculation_supp.png` | iToL, two rings from `Dataset_Colourstrip_BACcom_top15{,_noBAC2}.txt`; crosswalk by `scripts/crosswalk_old_vs_new_communities.R` | ✅ |
| S — assembly confounder | `assembly_quality_confounder.png` | `BipartiteNetwork_Analysis.ipynb` | ✅ |
| S — prophage length bimodality | `prophage_length_bimodality.png` | `BipartiteNetwork_Analysis.ipynb` | ✅ |

¹ Schematics/infographics, not data figures — no code is expected.

² Run at <https://rhea.icbm.uni-oldenburg.de/VIRIDIC/> on `Top15Medoids_VIRIDIC_input.fasta`.
Input, similarity/cluster tables and the heatmap are all preserved. **Note:** the published
heatmap labels the ICE community as `PRO_99`; it is annotated as "the ICE community" in
`VIRIDIC_top15_representatives.txt`.

³ Tree from the Roary core-gene alignment + FastTree (Chapter 3), rendered in **iToL v6**
with the top-15 bacterial-community colourstrip and prophage-community external-shape
datasets. Both are generated by `scripts/make_itol_community_datasets.R`; the treefile and
both datasets are in `data/itol/`.

⁴ iToL renderings driven by `DATASET_CONNECTION` linkage files. The original linkage files
were produced ad hoc and no generating script survived, so
`scripts/make_itol_crispr_linkages.py` reconstructs them from `Protospacer_hits.csv` +
`genome_community_mapping_49_BAC.csv`. The derivation was reverse-engineered from the
originals and **verified to reproduce them exactly** (communities 2, 3 and 5: 12,105 /
35,844 / 5,149 connections, row-for-row identical):

- source genome = spacer ID before `|`; target genome = provirus ID before `_`
- line width = `1.0 + 0.45 × (hits between that genome pair)`
- colour = grey `#cccccc` within community · red `#e31a1c` inside→outside ·
  blue `#1f78b4` outside→inside

Load the treefile from `data/itol/` plus the generated `Linkage_community_*.txt` in iToL.
The originals are also kept in `data/itol/per_community/` for reference.

⁵ Two-panel composite assembled in a vector editor. Panel A from `fig_modules_min5.png`,
panel B from `fig_nested_min5.png`. **Every statistic in the caption verified exactly**
against `tables/network_{modularity,nestedness}_results.csv` — Q=0.725, 12 modules,
null 0.637±0.014, z=+6.4; NODF=10.0, null 7.46±0.68, z=+3.8.

⁶ Composite of `Rarefaction_by_genome_size.png` + `Rarefaction_projected_forward.png`.
(The thesis filename is spelled `Rarefraction`.)

⁷ `Methodology.md` records the `sourmash compare` commands exactly, but not the
`sourmash plot` call that produced the two matrix PDFs.

## Naming caveat

Thesis figure files carry a `Draft3_` prefix added when copying into the thesis
`Figures/` folder. The scripts write the un-prefixed base name.

## Pipeline provenance (`pipeline/`)

`Methodology.md` is a 276-line step-by-step record of the upstream HPC work with real
SBATCH blocks — sourmash sketching/comparison (k=21 s=1000 for genomes, k=15 s=500 for
proviruses), geNomad `end-to-end --conservative --sensitivity 5`, provirus extraction,
MINCED v0.4.2 spacer detection, and the spacer→protospacer BLAST. `crispr_spacer_pipeline.sh`
is the runnable SLURM version of the CRISPR step.

## Inputs (`data/`)

| File | Size | Used by |
|---|---|---|
| `Bipartite_BAC_PRO_genome_edgetable.csv` | 51 M | bipartite network, rarefaction |
| `spacer_vs_provirus.raw.csv` | 48 M | CRISPR spacer analysis |
| `Protospacer_hits.csv` | 7.8 M | CRISPR targeting |
| `SupplementaryTable1-BVBRC_29Nov23.csv` | 2.4 M | genome metadata |
| `all_proviruses_QC.tsv` | 804 K | provirus QC |
| `all_spacers.fasta` | 688 K | spacer sequences |
| `genome_community_mapping_{44_PRO,49_BAC}.csv` | 200 K | Netpass community assignments |
| `community_questions/q1–q5*.csv` | — | per-question summaries |
| `confounders/{genome_master_table,prophage_length_table}.csv` | 532 K | confounder analyses |
