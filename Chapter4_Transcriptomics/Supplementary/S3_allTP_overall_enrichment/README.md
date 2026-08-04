# Supplementary S3 — All-timepoint overall enrichment dotplot

**Script:** `script/allTP_overall_enrichment_figure.R`
**Upstream scripts:** `script/COG-analysis-allTP.R`, `script/KEGG-pathway-analysis-allTP.R`, `script/GO-Enrichment-allTP.R`
**Figure:** `figure/allTP_overall_enrichment_figure.png`

Combined COG / KEGG / GO functional enrichment for every strain × phage, using
the union of DEGs across all post-infection timepoints (one Fisher/hypergeometric
test per term). Three stacked dotplot panels; y = term, x = condition, colour =
p-value, size = DEG count. Only terms passing BH FDR < 0.05 in ≥1 condition shown.

## Data
`data/{STRAIN}/{PHAGE}/{STRAIN}_{PHAGE}_allTP_{COG,KEGG,GO}_results/…csv` — the
per-condition enrichment results read by the figure script (18 CSVs).

### Inputs to the enrichment tests (`data/_inputs/`)
The upstream enrichment scripts (`COG/KEGG/GO-…-allTP.R`) consume these. Bundled
here so the enrichment can be re-run from the package alone:
- `_inputs/functional_annotation_{STRAIN}_with_names.csv` (3 files) — eggNOG
  functional annotation; `KEGG_Pathway_names` / COG / GO columns map genes → terms.
- `_inputs/DESeq2_vs_control/DESeq2_{STRAIN}_{PHAGE}_vs_CTLR_T{t}min.csv` (36 files)
  — host PHAGE-vs-CTLR DESeq2 results, all timepoints (identical to the Figure 2/3
  copies).

Note: the scripts currently point at absolute paths (`…/2_Analysis/Maria_Pipeline/…`
and `…/2_Analysis/results/DESeq2_vs_control/…`); repoint them at `data/_inputs/`
to run from within the package.

## Regeneration history
This figure was originally produced **before the C67 functional annotations were
available**, leaving C67 incomplete. On **2026-06-08** the C67 upstream
enrichments were re-run with the now-available
`functional_annotation_C67_with_names.csv`:

```
Rscript COG-analysis-allTP.R C67 {BON,CLY}
Rscript KEGG-pathway-analysis-allTP.R C67 {BON,CLY}
Rscript GO-Enrichment-allTP.R C67 {BON,CLY}
Rscript allTP_overall_enrichment_figure.R
```

C67 KEGG enrichment grew from 2.2 KB → 11.7 KB and C67 is now fully populated in
all three panels. (Also: line 34 of `allTP_overall_enrichment_figure.R` was made
conditional — `if (!requireNamespace("ggh4x")) install.packages(...)` — so it no
longer errors under `Rscript`.)
