# Supplementary S4 — Prophage region coverage

**Script:** `script/prophage_coverage_plot.R`
**Figures:** `figure/fig_prophage_coverage_{C67,D32,D68}.pdf` (+ `.png`)

Read-coverage across each strain's resident prophage region over the infection
timecourse, normalised to depth per million bases mapped (RPM-like), replicates
averaged with a ribbon. Prophage boundaries (shaded box) are hardcoded in the
script: C67 55,949–69,028; D32 & D68 5,291–18,370.

## Data
`data/coverage_prophage/` — per-sample region depth TSVs
(`{tp}-{strain}-{condition}-{rep}.tsv`) plus `normalization_factors.tsv`
(total bases mapped per sample, used for normalisation).

> Original location: `data/coverage/prophage/`.
