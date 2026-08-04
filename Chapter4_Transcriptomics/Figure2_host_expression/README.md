# Figure 2 — Host gene expression trajectories

**Script:** `script/fig2_host_trajectories_FINALISED.R`
**Upstream DE script:** `script/DESeq2_vs_control.R` — generated the `data/DESeq2_vs_control/` CSVs (set STRAIN/PHAGE at top; design `~ group`, low-count filter `rowSums ≥ 10`, each phage vs time-matched CTLR per timepoint, `lfcShrink` ashr).
**Figure:** `figure/fig2_host_trajectories.pdf`

Strip plot: one dot per gene per timepoint (y = log2FC, x = timepoint). Only
genes with padj < 0.05 shown. Colour: grey (|LFC| < 1), up/down colour
(|LFC| ≥ 1); defence/prophage genes always coloured (locus-tag lists hardcoded
in the script). Layout: 3 strains (rows) × Bonnie/Clyde (cols).

## Data
- `data/DESeq2_vs_control/DESeq2_{STRAIN}_{PHAGE}_vs_CTLR_T{2,10,20,30,50}min.csv`
  — host PHAGE-vs-CTLR DESeq2 results per timepoint
  (cols: gene, baseMean, log2FoldChange, lfcSE, pvalue, padj).

T=0 excluded (same samples across conditions). No genes removed beyond the padj
filter applied at plot time.
