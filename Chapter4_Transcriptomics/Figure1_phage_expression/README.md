# Figure 1 — Phage temporal expression

**Script:** `script/fig1_phage_expression_FINALISED.R`
**Figure:** `figure/fig1_phage_expression.pdf`

Panel A: CSR ternary of Early/Middle/Late gene-class proportions per strain ×
phage (class counts hardcoded in script; source = PEA `ClassThreshold`).
Panels B–G: ridgeline/joy plots of phage gene expression across the timecourse,
one panel per strain × phage in genomic order.

## Data
- `data/raw_counts/{STRAIN}_{PHAGE}_full_raw_counts.tsv` — per-condition count
  matrices (18 count cols: 6 timepoints × 3 reps + Entity + Symbol).
- `data/PEA/{STRAIN}_{PHAGE}_fractional_expression.tsv` — `ClassThreshold`
  column (Early/Middle/Late/None) used for the ternary.

## Method (see script header for full detail)
T=0 dropped; phage rows only (prefix `PENJXGPI_CDS_` Bonnie / `ZUDWSPYW_CDS_`
Clyde); rowSum<10 excluded; VST per dataset; reps averaged; per-gene 0–1 scaling;
cubic-spline smoothing.
