# Chapter 4 — Host transcriptional response to phage infection

> **Thesis repo note.** This is the curated figure package, organised one directory per
> figure. Two additions and three caveats:
>
> **Added in `scripts_extra/`** — two thesis figures that were not in the original package:
> `temporal_directional_enrichment.R` (all-timepoint enrichment figure, recovered from
> `Maria_Pipeline/scripts/`) and `DESeq2_timecourse_overview.R` (`grid_pca_3x3.png`,
> from `2_Analysis/scripts/R/`).
>
> **⚠️ Scripts are archival, not runnable in place.** As noted below, they use paths
> relative to the original analysis root, not to this folder. Reproducing a figure means
> restoring those paths, not just running the script here.
>
> **Two figures have no code.** `S1_PhagesvsProphagesClinker.png` was made with clinker
> v0.0.32 interactively; `Galactose_ClydeEffect.png` has no generating script anywhere in
> the project.
>
> **Supplementary numbering differs.** Thesis Supplementary Figure S2 (prophage coverage)
> lives here as `Supplementary/S4_prophage_coverage/`.
>
> **Upstream read processing is deliberately absent** — it ran through the published
> **PhageExpressionAtlas** pipeline (FastQC v0.12.1, Cutadapt v4.9, HISAT2 v2.2.1,
> SAMtools v1.6, featureCounts v2.0.1, MultiQC v1.32), which has its own documentation.
> Annotation tools (Bakta, Pharokka, Phold, eggNOG-mapper, geNomad, Ptolemaea) have
> versions and parameters in the chapter's Table `tab:Tools`, but no run scripts survived.

---

## Original package documentation

Self-contained, traceable bundle of every main and supplementary figure: each
folder holds the **figure**, the **script** that produced it, and **all input
data** it consumes.

Built by copying from the working tree on **2026-06-08**; originals are
untouched. Scripts use project-relative paths (`results/...`, `data/...`) and
were written to run from the analysis root
(`Project_Chp3_Transcriptomics/2_Analysis/`); the copies here are for archival
and traceability. Paths to the original inputs are listed per figure below.

## Experimental design
- **Strains (host):** C67, D32, D68
- **Conditions:** BON (Bonnie phage), CLY (Clyde phage), CTLR (no phage)
- **Timepoints:** 0, 2, 10, 20, 30, 50 min post-infection · 3 replicates
- **DEG definition:** padj < 0.05 AND |log2FC| ≥ 1 (vs uninfected control)

## Main figures

| Figure | Folder | Script | Original figure | Key inputs (original location) |
|---|---|---|---|---|
| **Fig 1** — Phage temporal expression (CSR ternary + ridgelines) | `Figure1_phage_expression/` | `fig1_phage_expression_FINALISED.R` | `results/figures/fig1_phage_expression.pdf` | `data/raw_counts/{S}_{P}_full_raw_counts.tsv`; `results/PEA/{S}_{P}_fractional_expression.tsv` |
| **Fig 2** — Host gene DE trajectories | `Figure2_host_expression/` | `fig2_host_trajectories_FINALISED.R` | `results/figures/fig2_host_trajectories.pdf` | `results/DESeq2_vs_control/DESeq2_{S}_{P}_vs_CTLR_T{t}min.csv` |
| **Fig 3** — Top-3 KEGG pathway enrichment (Ribosome, Galactose metabolism, PTS) | `Figure3_KEGG_enrichment/` | `temporal_directional_enrichment_mainKEGG.R` | `results/temporal_enrichment_vsCTLR/temporal_directional_enrichment_mainKEGG.png` | `Maria_Pipeline/functional_annotation_{S}_with_names.csv`; `results/DESeq2_vs_control/...vs_CTLR...csv` |

## Supplementary figures

| Suppl. | Folder | Script | Figure | Key inputs |
|---|---|---|---|---|
| **S1** — Bonnie genome map | `Supplementary/S1_bonnie_genome_map/` | `phage_genome_map.R` | `bonnie_genome_map.pdf` | `Bonnie_allannos.gff3` |
| **S2** — Clyde genome map | `Supplementary/S2_clyde_genome_map/` | `phage_genome_map.R` | `clyde_genome_map.pdf` | `Clyde_allannos.gff3` |
| **S3** — All-timepoint overall enrichment dotplot | `Supplementary/S3_allTP_overall_enrichment/` | `allTP_overall_enrichment_figure.R` (+ upstream COG/KEGG/GO-allTP scripts) | `allTP_overall_enrichment_figure.png` | per-condition COG/KEGG/GO `*_allTP_*_results` CSVs |
| **S4** — Prophage region coverage | `Supplementary/S4_prophage_coverage/` | `prophage_coverage_plot.R` | `fig_prophage_coverage_{C67,D32,D68}.pdf/png` | `data/coverage/prophage/*.tsv` (+ `normalization_factors.tsv`) |
| **S5** — Bonnie-vs-Clyde differential expression enrichment | `Supplementary/S5_BONvsCLY_enrichment/` | `BONvsCLY_directional_enrichment.R` (+ `DESeq2_BONvsCLY.R`) | `BONvsCLY_directional_enrichment_figure.png` | `results/deseq2_BONvsCLY/DESeq2_{S}_BONvsCLY_T{t}min.csv`; `functional_annotation_{S}_with_names.csv` |

`{S}` = strain (C67/D32/D68), `{P}` = phage (BON/CLY), `{t}` = timepoint.

## Notes
- **S3 (allTP dotplot) was regenerated on 2026-06-08** after the C67 functional
  annotations became available. The C67 upstream COG/KEGG/GO-allTP enrichments
  were re-run (C67 KEGG enrichment grew from 2.2 KB → 11.7 KB), then the figure
  rebuilt — C67 is now fully populated in all three panels. The refreshed inputs
  and figure are in `S3_allTP_overall_enrichment/`.
