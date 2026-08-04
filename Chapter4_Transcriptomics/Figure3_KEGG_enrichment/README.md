# Figure 3 — Top-3 KEGG pathway enrichment (temporal, directional)

**Script:** `script/temporal_directional_enrichment_mainKEGG.R`
**Upstream DE script:** `script/DESeq2_vs_control.R` — generated the `data/DESeq2_vs_control/` CSVs this enrichment reads (design `~ group`, low-count filter `rowSums ≥ 10`, each phage vs time-matched CTLR per timepoint, `lfcShrink` ashr).
**Figure:** `figure/temporal_directional_enrichment_mainKEGG.png`
**Underlying values:** `figure/temporal_KEGG_mainText_directional_table.csv`

Three main-text KEGG pathways — **Ribosome**, **Galactose metabolism**,
**Phosphotransferase system (PTS)** — one panel each with a real 0–100% y-axis.
Band = % of the pathway's tested genes that are DE vs uninfected control
(padj < 0.05, |log2FC| ≥ 1) at each timepoint; above 0 = induced, below =
repressed. Colour = strain hue (darker induced, lighter repressed); columns =
strain × phage. A colour legend is included.

This is the pared-down main-text version of the full COG+KEGG figure
(`temporal_directional_enrichment.R`), which remains available as a broader
supplementary if needed.

## Data
- `data/functional_annotation_{STRAIN}_with_names.csv` — eggNOG functional
  annotation; `KEGG_Pathway_names` column maps genes → pathways.
- `data/DESeq2_vs_control/DESeq2_{STRAIN}_{PHAGE}_vs_CTLR_T{t}min.csv` — host
  PHAGE-vs-CTLR DESeq2 results (same files as Figure 2).

Pathway selection is hardcoded in the script (`SEL_KEGG`); change that vector
to add/remove pathways (e.g. to bring back Starch and sucrose metabolism).
