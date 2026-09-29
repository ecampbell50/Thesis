# Chapter 4 — Transcriptomics

Strains C67, D32, D68; phages Bonnie (BON) and Clyde (CLY); uninfected control (CTLR).
Timepoints 0, 2, 10, 20, 30, 50 min, 3 replicates. DEGs: padj < 0.05 and |log2FC| ≥ 1 vs control.

Scripts were run from my original analysis folder, so input paths may need changing.

| Thesis | Code |
|---|---|
| Fig 4.2 | `Figure1_phage_expression/` |
| Fig 4.3 | `Figure2_host_expression/` (DESeq2: `script/DESeq2_vs_control.R`) |
| Fig 4.4 | `Figure3_KEGG_enrichment/` |
| Fig S4.1, S4.2 | `Supplementary/S1_bonnie_genome_map/`, `S2_clyde_genome_map/` |
| Fig S4.3 | `hpc/13`–`16` (clinker) |
| Fig S4.4–S4.6 | `hpc/04`, then `Supplementary/S4_prophage_coverage/` |
| Fig S4.7, S4.8 | `scripts_extra/temporal_directional_enrichment.R` |
| Fig S4.12 | `scripts_extra/DESeq2_timecourse_overview.R` |
| Tables S4.1, S4.3 | `hpc/10`, `hpc/11` |

`hpc/` has the Kelvin2 jobs: read QC, PhageExpressionAtlas runs, and genome annotation.
