# Chapter 3 — Anti-phage defence profile of *Streptococcus suis*

Code and data for every figure, table and statistic in Chapter 3.

```
Chapter3_DefenceProfile/
├── data/                        all inputs (nothing outside this folder is read)
├── Ssuis_DefenceAnalysis.ipynb  Python: Figures 2–5, Table 1, all Results numbers
├── defence_network.R            R: Figure 1 + Figure S2 (protein-similarity network)
├── supporting_scripts/          methods-supporting scripts (medoids, tmn, iToL)
├── figures/                     generated figures
└── tables/                      generated statistic tables
```

## Running

Both scripts run **from this directory** so the relative `data/` paths resolve.

```bash
jupyter nbconvert --to notebook --execute --inplace Ssuis_DefenceAnalysis.ipynb
```
```bash
Rscript defence_network.R
```

Python 3.9+ with `pandas numpy matplotlib scipy scikit-learn statsmodels`
(notebook runs in ~30 s). R with `igraph ggraph tidyverse ggrepel`, installed on
first run if missing.

## Figure → code map

| Figure | Output file | Produced by |
|---|---|---|
| **Fig 1** — defence landscape | `figures/network/<system>/` | `defence_network.R` ¹ |
| **Fig 2** — richness vs gene count | `figures/OutcomeVsCount.pdf` | notebook §4 |
| **Fig 3** — prevalence heatmap + histogram | `figures/HeatmapAndHistogram_noUbiq.pdf` | notebook §6 (markers from §5) |
| **Fig 4** — PCA of defence profiles | `figures/PCA_Draft3.png` | notebook §7 |
| **Fig 5** — co-occurrence heatmap | `figures/defence_cooccurrence_heatmap.png` | notebook §8 |
| **Fig 6** — core-genome phylogeny | — | **external** ² |
| **Fig 7** — *Tmn* reconciliation | — | **external** ³ |
| **Fig S2** — full defence network | `figures/network/all/` | `defence_network.R` |
| **Fig S5 / S6** — *Tmn* β/α reconciliation | — | **external** ³ |
| **Table 1** — subtypes >99% | printed by notebook §3 | notebook §3 |

¹ The script produces the **components** — one network plus per-community panels for
each defence type with ≥100 genes. The published Figure 1 is a composite assembled
from these in a vector editor, so the final layout is not reproduced by code.

² Tree inferred with Roary (core gene alignment) + FastTree, then **rendered in iToL**.
Annotation files are in `supporting_scripts/itol/`; the rendering itself has no code.

³ Produced with GeneRax / ThirdKind / AlphaFold. External tools, no plotting code.

## Inputs (`data/`)

| File | Used by | Description |
|---|---|---|
| `ssuis_defence_CLEANED.csv` | notebook | per-genome defence gene-count matrix — **the main input** |
| `Ssuis_BVBRC_29Nov23.csv` | notebook | genome metadata (size, CDS, status); BV-BRC retrieved 29 Nov 2023 (NB: >2,119 IDs, some not downloaded when retrieving) |
| `ConsensusSerotypes_23Jul24.csv` | notebook | per-genome serotype |
| `serotype_colours.csv` | notebook | serotype palette for Figs 3 & 4 |
| `SmallCluster_Tmngenes.txt` | notebook | small CLANS cluster = *tmn*-beta |
| `LargeCluster_Tmngenes.txt` | notebook | large CLANS cluster = *tmn*-alpha |
| `Ssuis_rawptolemaeaoutput.csv` | R | per-protein annotations (97,099 proteins) |
| `network_edges.tsv` | R | filtered all-vs-all MMseqs2 hits → edges |
| `node_sizes.tsv` | R | representative → cluster size (node area) |
| `cluster_membership.tsv` | R | representative → member protein |

### Important: `ssuis_defence_CLEANED.csv` has a single, unsplit `tmn` column

This is correct and deliberate. The notebook performs the *tmn*-α/β split itself in
§2b, using `SmallCluster_Tmngenes.txt`, and §2a asserts that exactly one `tmn` column
is present. **A pre-split matrix will make the notebook fail.** Any file elsewhere in
the project named `*withTmnSplit*` or `*CLEAN_FINAL*` is a stale downstream artefact,
not an input.

## Two notes on the numbers

**Subtype counts, 148 vs 149.** The chapter text reports 74 types / 144 subtypes /
148 systems (outcome split 43/19/10/67/9). The notebook prints 149 / 145 / 68 because
it counts *after* the *tmn*-α/β split. Both are correct: the text describes the
annotation output, and the chapter only introduces the α/β split later.

**Protein counts in the network, 29,535 vs 29,521.** `defence_network.R` filters with a
case-insensitive **substring** match, so `"RM"` also captures composite annotations such
as `(p::Other|d::RM)` and `(p::RM|d::PrrC)` — 5 extra clusters, 14 proteins. This is
intentional (a cluster is shown if *any* member matches) and is why the RM panel reports
29,535 proteins, of which 27,200 are exact-type RM and 2,301 are DMS_other that
co-cluster with RM.
