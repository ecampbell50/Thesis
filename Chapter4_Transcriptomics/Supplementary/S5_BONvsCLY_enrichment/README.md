# Supplementary S5 — Bonnie-vs-Clyde differential expression enrichment

**Script:** `script/BONvsCLY_directional_enrichment.R`  (enrichment + figure)
**Upstream:** `script/DESeq2_BONvsCLY.R`  (produces the BONvsCLY DESeq2 CSVs)
**Figure:** `figure/BONvsCLY_directional_enrichment_figure.png`

Functional enrichment (COG/KEGG/GO) of the **direct** Bonnie-vs-Clyde host DE
contrast, split by direction. Union of significant genes across timepoints
(padj < 0.05, |LFC| ≥ 1) split into Bonnie-up (LFC ≥ +1) vs Clyde-up (LFC ≤ −1);
one hypergeometric enrichment per direction. x = direction, facets = strain,
panels = COG/KEGG/GO.

> **Interpretation caveat:** Clyde is the reference, so this is a phage-vs-phage
> contrast with NO uninfected baseline. "Bonnie-up" = higher under Bonnie *than
> Clyde*, not induced vs resting.

## Data
- `data/deseq2_BONvsCLY/DESeq2_{STRAIN}_BONvsCLY_T{2,10,20,30,50}min.csv`
- `data/functional_annotation_{STRAIN}_with_names.csv` (TERM2GENE source)
