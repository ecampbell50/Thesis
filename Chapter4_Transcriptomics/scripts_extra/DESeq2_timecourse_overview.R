# DESeq2_timecourse_overview.R
# Runs all 9 strain × phage combinations and saves three 3×3 comparison grids:
#   grid_volcano_3x3.png   — faceted volcano plots
#   grid_pca_3x3.png       — PCA plots
#   grid_heatmap_3x3.png   — sample distance heatmaps
#
# Layout of each grid:
#              BON (Bonnie)   CLY (Clyde)   CTLR (Control)
#   C67  │  C67-BON       │  C67-CLY    │  C67-CTLR
#   D32  │  D32-BON       │  D32-CLY    │  D32-CTLR
#   D68  │  D68-BON       │  D68-CLY    │  D68-CTLR

# ── Package installation (run once, then comment out) ─────────────────────────
if (!("DESeq2"    %in% installed.packages())) BiocManager::install("DESeq2",  update = FALSE)
if (!("ashr"      %in% installed.packages())) BiocManager::install("ashr",    update = FALSE)
if (!("ggrepel"   %in% installed.packages())) install.packages("ggrepel")
if (!("RColorBrewer" %in% installed.packages())) install.packages("RColorBrewer")
if (!("patchwork" %in% installed.packages())) install.packages("patchwork")

library(DESeq2)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# Thresholds for labelling extreme genes on volcano plots
LFC_THRESHOLD  <- 3
PADJ_THRESHOLD <- 1e-10

# Counts files (TSV, gene IDs as row names)
COUNTS_FILES <- list(
  C67 = "C67_DESeq2data.tsv",
  D32 = "D32_DESeq2data.tsv",
  D68 = "D68_DESeq2data.tsv"
)

# Output filenames (saved in working directory)
OUT_VOLCANO <- "grid_volcano_3x3.png"
OUT_PCA     <- "grid_pca_3x3.png"
OUT_HEATMAP <- "grid_heatmap_3x3.png"

# Samples to exclude per combo — paste sample names as strings, e.g.:
#   C67_CLY = c("0_C67_CLY_R1", "0_C67_CLY_R3")
# Set to character(0) if no exclusions for that run.
EXCLUDE_SAMPLES <- list(
  C67_BON  = character(0),
  C67_CLY  = character(0),
  C67_CTLR = character(0),
  D32_BON  = character(0),
  D32_CLY  = character(0),
  D32_CTLR = character(0),
  D68_BON  = character(0),
  D68_CLY  = character(0),
  D68_CTLR = character(0)
)

# ══════════════════════════════════════════════════════════════════════════════
# DEFENCE GENE DICTIONARIES
# (named vector: locus_tag = "SystemName_SxGy")
# ══════════════════════════════════════════════════════════════════════════════

DEFENCE_GENES <- list(

  C67 = c(
    # ── Original ────────────────────────────────────────────────────────────
    "C67_prot_00121" = "PDC-M22_S1G1",
    "C67_prot_00122" = "PDC-M22_S1G2",
    "C67_prot_00133" = "RM_Type_IV_S1G1",
    "C67_prot_00134" = "RM_Type_IV_S1G2",
    "C67_prot_00284" = "RM_Type_I_S1G1",
    "C67_prot_00285" = "RM_Type_I_S1G2",
    "C67_prot_00286" = "RM_Type_I_S1G3",
    "C67_prot_00287" = "RM_Type_I_S1G4",
    "C67_prot_00297" = "RM_Type_III_S1G1",
    "C67_prot_00298" = "RM_Type_III_S1G2",
    "C67_prot_01021" = "RM_Type_I_S2G1",
    "C67_prot_01022" = "RM_Type_I_S2G2",
    "C67_prot_01023" = "RM_Type_I_S2G3",
    "C67_prot_01024" = "RM_Type_I_S2G4",
    "C67_prot_01211" = "PDC-S30_S1G1",
    "C67_prot_01330" = "RosmerTA_S1G1",
    "C67_prot_01331" = "RosmerTA_S1G2",
    "C67_prot_01418" = "Ogmios_S1G1",
    "C67_prot_01435" = "tmn_S1G1",
    "C67_prot_01476" = "PD-T4-6_S1G1",
    # ── New from Ptolemaea ───────────────────────────────────────────────────
    "C67_prot_00060"  = "PDC-S07_S1G1",
    "C67_prot_00603"  = "PDC-S07_S2G1",
    "C67_prot_01356"  = "PDC-S07_S3G1",
    "C67_prot_00253"  = "Gabija_S1G1",
    "C67_prot_00763"  = "CAS_Class1-IV_S1G1",
    "C67_prot_01000"  = "DRT_other_S1G1",
    "C67_prot_01545"  = "DRT_other_S2G1",
    "C67_prot_01598"  = "Theoris_S1G1",
    "C67_prot_00210"  = "DMS_other_S2G1",
    "C67_prot_00597"  = "DMS_other_S3G1",
    "C67_prot_00619"  = "DMS_other_S4G1",
    "C67_prot_01170"  = "DMS_other_S5G1",
    "C67_prot_01206"  = "DMS_other_S6G1",
    "C67_prot_01226"  = "DMS_other_S7G1",
    "C67_prot_01617"  = "DMS_other_S8G1",
    "C67_prot_01819"  = "RM_Type_I_S3G1"
  ),

  D32 = c(
    # ── CAS_Class1-IV ───────────────────────────────────────────────────────
    "D32_prot_00742" = "CAS_Class1-IV_S1G1",
    # ── DMS_other ───────────────────────────────────────────────────────────
    "D32_prot_00009" = "DMS_other_S1G1",
    "D32_prot_00085" = "DMS_other_S2G1",
    "D32_prot_00247" = "DMS_other_S3G1",
    "D32_prot_00257" = "DMS_other_S3G2",
    "D32_prot_00379" = "DMS_other_S6G1",
    "D32_prot_00401" = "DMS_other_S7G1",
    "D32_prot_00920" = "DMS_other_S8G1",
    "D32_prot_00940" = "DMS_other_S9G1",
    "D32_prot_01852" = "DMS_other_S4G1",
    "D32_prot_01911" = "DMS_other_S5G1",
    # ── DRT_other ───────────────────────────────────────────────────────────
    "D32_prot_01263" = "DRT_other_S1G1",
    "D32_prot_01431" = "DRT_other_S2G1",
    # ── Gabija ──────────────────────────────────────────────────────────────
    "D32_prot_00042" = "Gabija_S1G1",
    # ── Ogmios ──────────────────────────────────────────────────────────────
    "D32_prot_01386" = "Ogmios_S1G1",     # ⚠ MAPPING-status only
    # ── PD-T4-6 ─────────────────────────────────────────────────────────────
    "D32_prot_01328" = "PD-T4-6_S1G1",
    # ── PDC-M22 ─────────────────────────────────────────────────────────────
    "D32_prot_00172" = "PDC-M22_S1G1",
    "D32_prot_00173" = "PDC-M22_S1G2",
    # ── PDC-S07 ─────────────────────────────────────────────────────────────
    "D32_prot_00292" = "PDC-S07_S1G1",
    "D32_prot_00395" = "PDC-S07_S2G1",
    "D32_prot_01175" = "PDC-S07_S3G1",
    # ── PDC-S30 ─────────────────────────────────────────────────────────────
    "D32_prot_00925" = "PDC-S30_S1G1",
    # ── RM_Type_I ───────────────────────────────────────────────────────────
    "D32_prot_00008" = "RM_Type_I_S1G1",
    "D32_prot_00010" = "RM_Type_I_S1G2",
    "D32_prot_00011" = "RM_Type_I_S1G3",
    "D32_prot_00882" = "RM_Type_I_S2G1",
    "D32_prot_00883" = "RM_Type_I_S2G2",
    "D32_prot_00884" = "RM_Type_I_S2G3",
    "D32_prot_01452" = "RM_Type_I_S3G1",
    "D32_prot_01453" = "RM_Type_I_S3G2",
    "D32_prot_01454" = "RM_Type_I_S3G3",
    "D32_prot_01455" = "RM_Type_I_S3G4",
    # ── RM_Type_III ─────────────────────────────────────────────────────────
    "D32_prot_01970" = "RM_Type_III_S1G1",
    "D32_prot_01971" = "RM_Type_III_S1G2",
    # ── RM_Type_IV ──────────────────────────────────────────────────────────
    "D32_prot_00160" = "RM_Type_IV_S1G1",
    "D32_prot_00161" = "RM_Type_IV_S1G2",
    # ── RosmerTA ────────────────────────────────────────────────────────────
    "D32_prot_01200" = "RosmerTA_S1G1",
    "D32_prot_01201" = "RosmerTA_S1G2",
    # ── Theoris ─────────────────────────────────────────────────────────────
    "D32_prot_01867" = "Theoris_S1G1",
    # ── tmn ─────────────────────────────────────────────────────────────────
    "D32_prot_01369" = "tmn_S1G1"         # ⚠ MAPPING-status only
  ),

  D68 = c(
    # ── Abi2DF ──────────────────────────────────────────────────────────────
    "D68_prot_01256" = "Abi2DF_S1G1",
    "D68_prot_01327" = "Abi2DF_S2G1",
    # ── CAS_Class1-IV ───────────────────────────────────────────────────────
    "D68_prot_00729" = "CAS_Class1-IV_S1G1",
    # ── DMS_other ───────────────────────────────────────────────────────────
    "D68_prot_00232" = "DMS_other_S1G1",
    "D68_prot_00254" = "DMS_other_S2G1",
    "D68_prot_00291" = "DMS_other_S4G1",
    "D68_prot_00367" = "DMS_other_S5G1",
    "D68_prot_00528" = "DMS_other_S7G1",
    "D68_prot_00548" = "DMS_other_S8G1",
    "D68_prot_00919" = "DMS_other_S9G1",
    "D68_prot_00929" = "DMS_other_S9G2",
    "D68_prot_01833" = "DMS_other_S3G1",
    "D68_prot_01915" = "DMS_other_S6G1",
    # ── DRT_other ───────────────────────────────────────────────────────────
    "D68_prot_01201" = "DRT_other_S2G1",
    "D68_prot_01696" = "DRT_other_S1G1",
    # ── Gabija ──────────────────────────────────────────────────────────────
    "D68_prot_00324" = "Gabija_S1G1",
    # ── Ogmios ──────────────────────────────────────────────────────────────
    "D68_prot_01102" = "Ogmios_S1G1",     # ⚠ MAPPING-status only
    # ── PD-T4-6 ─────────────────────────────────────────────────────────────
    "D68_prot_01160" = "PD-T4-6_S1G1",
    # ── PDC-S07 ─────────────────────────────────────────────────────────────
    "D68_prot_00237" = "PDC-S07_S1G1",
    "D68_prot_00238" = "PDC-S07_S1G2",
    "D68_prot_01046" = "PDC-S07_S2G1",
    # ── PDC-S30 ─────────────────────────────────────────────────────────────
    "D68_prot_00533" = "PDC-S30_S1G1",
    # ── RM_Type_I ───────────────────────────────────────────────────────────
    "D68_prot_00289" = "RM_Type_I_S1G1",
    "D68_prot_00292" = "RM_Type_I_S1G2",
    "D68_prot_00293" = "RM_Type_I_S1G3",
    "D68_prot_00490" = "RM_Type_I_S2G1",
    "D68_prot_00491" = "RM_Type_I_S2G2",
    "D68_prot_00492" = "RM_Type_I_S2G3",
    "D68_prot_01222" = "RM_Type_I_S3G1",
    "D68_prot_01223" = "RM_Type_I_S3G2",
    "D68_prot_01224" = "RM_Type_I_S3G3",
    "D68_prot_01225" = "RM_Type_I_S3G4",
    # ── RM_Type_III ─────────────────────────────────────────────────────────
    "D68_prot_01973" = "RM_Type_III_S1G1",
    "D68_prot_01974" = "RM_Type_III_S1G2",
    # ── RM_Type_IV ──────────────────────────────────────────────────────────
    "D68_prot_00442" = "RM_Type_IV_S1G1",
    "D68_prot_00443" = "RM_Type_IV_S1G2",
    # ── RosmerTA ────────────────────────────────────────────────────────────
    "D68_prot_01020" = "RosmerTA_S1G1",
    "D68_prot_01021" = "RosmerTA_S1G2",
    # ── Theoris ─────────────────────────────────────────────────────────────
    "D68_prot_01871" = "Theoris_S1G1",
    # ── tmn ─────────────────────────────────────────────────────────────────
    "D68_prot_01119" = "tmn_S1G1"         # ⚠ MAPPING-status only
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# MAIN FUNCTION — runs one strain/phage combo and returns the three plots
# ══════════════════════════════════════════════════════════════════════════════

run_combo <- function(strain, phage) {
  cat("\n\n══════════════════════════════════════\n")
  cat("  Running:", strain, "—", phage, "\n")
  cat("══════════════════════════════════════\n")

  defence_genes <- DEFENCE_GENES[[strain]]
  exclude_key   <- paste0(strain, "_", phage)
  exclude_samps <- EXCLUDE_SAMPLES[[exclude_key]]

  # ── Load & filter counts ──────────────────────────────────────────────────
  counts_raw           <- read.delim(COUNTS_FILES[[strain]], row.names = 1)
  colnames(counts_raw) <- sub("^X", "", colnames(counts_raw))

  counts_sub <- counts_raw[, grepl(phage, colnames(counts_raw)), drop = FALSE]
  counts_sub <- counts_sub[, !colnames(counts_sub) %in% exclude_samps, drop = FALSE]

  if (ncol(counts_sub) == 0) {
    stop("No columns matched phage '", phage, "' for strain '", strain,
         "'. Check COUNTS_FILES and PHAGE spelling.")
  }

  # ── Timepoints ────────────────────────────────────────────────────────────
  all_times     <- sort(unique(as.numeric(sub("_.*", "", colnames(counts_sub)))))
  all_times_chr <- as.character(all_times)
  timepoints    <- all_times_chr[all_times_chr != "0"]

  # ── DESeq2 ───────────────────────────────────────────────────────────────
  colData_sub <- data.frame(
    row.names = colnames(counts_sub),
    time = factor(sub("_.*", "", colnames(counts_sub)), levels = all_times_chr)
  )

  dds      <- DESeqDataSetFromMatrix(countData = counts_sub,
                                     colData   = colData_sub,
                                     design    = ~ time)
  keep     <- rowSums(counts(dds)) >= 10
  dds      <- dds[keep, ]
  dds$time <- relevel(dds$time, ref = "0")
  dds      <- DESeq(dds)

  # ── LFC shrinkage for each timepoint vs T=0 ───────────────────────────────
  results_list <- list()
  for (tp in timepoints) {
    results_list[[tp]] <- lfcShrink(dds,
                                    coef = paste0("time_", tp, "_vs_0"),
                                    type = "ashr")
  }

  # ── VST for PCA and heatmap ───────────────────────────────────────────────
  vsd <- vst(dds, blind = TRUE)

  # ──────────────────────────────────────────────────────────────────────────
  # PCA
  # ──────────────────────────────────────────────────────────────────────────
  time_colours <- c("0"  = "#2ca02c", "2"  = "#bcbd22",
                    "10" = "#ff7f0e", "20" = "#d62728",
                    "30" = "#9467bd", "50" = "#1f77b4")

  pca_data    <- plotPCA(vsd, intgroup = "time", returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, colour = time, label = name)) +
    geom_point(size = 3) +
    geom_text_repel(size = 2.5, max.overlaps = Inf, box.padding = 0.4,
                    segment.size = 0.3) +
    scale_colour_manual(values = time_colours) +
    labs(x      = paste0("PC1: ", percent_var[1], "%"),
         y      = paste0("PC2: ", percent_var[2], "%"),
         colour = "Time (min)",
         title  = paste0(strain, " — ", phage)) +
    theme_bw(base_size = 10) +
    theme(legend.position  = "bottom",
          legend.key.size  = unit(0.4, "cm"),
          legend.text      = element_text(size = 8),
          plot.title       = element_text(size = 11, face = "bold"),
          aspect.ratio     = 1)

  # ──────────────────────────────────────────────────────────────────────────
  # Sample distance heatmap
  # pheatmap(..., silent=TRUE) draws nothing but returns a list with $gtable
  # ──────────────────────────────────────────────────────────────────────────
  sampleDists      <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  heatmap_colours  <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

  hm <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows = sampleDists,
                 clustering_distance_cols = sampleDists,
                 col      = heatmap_colours,
                 main     = paste0(strain, " — ", phage),
                 fontsize = 7,
                 silent   = TRUE)

  # ──────────────────────────────────────────────────────────────────────────
  # Faceted volcano
  # ──────────────────────────────────────────────────────────────────────────
  all_results <- bind_rows(
    lapply(timepoints, function(tp) {
      as.data.frame(results_list[[tp]]) %>%
        rownames_to_column("gene") %>%
        filter(!is.na(padj)) %>%
        mutate(timepoint = tp)
    })
  ) %>%
    mutate(
      timepoint     = factor(timepoint, levels = timepoints),
      sig           = padj < 0.05 & abs(log2FoldChange) > 1,
      is_defence    = gene %in% names(defence_genes),
      defence_label = ifelse(is_defence, defence_genes[gene], NA_character_),
      is_extreme    = sig & !is_defence &
                      (abs(log2FoldChange) > LFC_THRESHOLD | padj < PADJ_THRESHOLD),
      category = case_when(
        is_defence ~ "defence",
        sig        ~ "significant",
        TRUE       ~ "not_significant"
      )
    )

  annot_right <- data.frame(
    timepoint = factor(timepoints, levels = timepoints),
    label     = paste0("Higher at T=", timepoints, "min")
  )
  annot_left <- data.frame(
    timepoint = factor(timepoints, levels = timepoints),
    label     = "Higher at T=0min"
  )

  volcano_facet <- ggplot(all_results, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = filter(all_results, category == "not_significant"),
               colour = "grey70", alpha = 0.4, size = 0.5) +
    geom_point(data = filter(all_results, category == "significant"),
               colour = "red", alpha = 0.6, size = 0.7) +
    geom_point(data = filter(all_results, is_defence),
               colour = "dodgerblue", size = 1.8, shape = 18) +
    geom_text_repel(data = filter(all_results, is_extreme),
                    aes(label = gene),
                    colour = "black", size = 1.8,
                    max.overlaps = 10, box.padding = 0.2,
                    fontface = "italic") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "black",
               linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black",
               linewidth = 0.3) +
    geom_text(data = annot_left,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.5, inherit.aes = FALSE,
              colour = "grey40", size = 2, fontface = "italic") +
    geom_text(data = annot_right,
              aes(x = Inf, y = Inf, label = label),
              hjust = 1.05, vjust = 1.5, inherit.aes = FALSE,
              colour = "grey40", size = 2, fontface = "italic") +
    facet_wrap(~ timepoint, nrow = 2,
               labeller = labeller(
                 timepoint = function(x) paste0("T=", x, "min vs T=0")
               )) +
    labs(title    = paste0(strain, " — ", phage),
         subtitle = paste0("red=sig.DE (padj<0.05, |LFC|>1)  |  ",
                           "blue=defence genes (n=", length(defence_genes), ")"),
         x = "Log2 Fold Change",
         y = "-log10(adjusted p-value)") +
    theme_bw(base_size = 9) +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 7,  colour = "grey40"),
          aspect.ratio  = 1,
          strip.text    = element_text(size = 8))

  return(list(pca     = pca_plot,
              heatmap = hm$gtable,
              volcano = volcano_facet))
}

# ══════════════════════════════════════════════════════════════════════════════
# RUN ALL 9 COMBINATIONS
# ══════════════════════════════════════════════════════════════════════════════
# Order defines the grid layout (left→right, top→bottom):
#   [1] C67-BON   [2] C67-CLY   [3] C67-CTLR
#   [4] D32-BON   [5] D32-CLY   [6] D32-CTLR
#   [7] D68-BON   [8] D68-CLY   [9] D68-CTLR

combos <- data.frame(
  strain = c("C67","C67","C67", "D32","D32","D32", "D68","D68","D68"),
  phage  = c("BON","CLY","CTLR","BON","CLY","CTLR","BON","CLY","CTLR"),
  stringsAsFactors = FALSE
)

results_all <- vector("list", nrow(combos))
for (i in seq_len(nrow(combos))) {
  results_all[[i]] <- run_combo(combos$strain[i], combos$phage[i])
}

cat("\n\nAll 9 runs complete. Assembling grids...\n")

# ══════════════════════════════════════════════════════════════════════════════
# ASSEMBLE AND SAVE 3×3 GRIDS
# ══════════════════════════════════════════════════════════════════════════════

# Shared annotation for all three grids
grid_subtitle <- "Columns: BON (Bonnie)  |  CLY (Clyde)  |  CTLR (Control)     Rows: C67  |  D32  |  D68"

# ── PCA grid ─────────────────────────────────────────────────────────────────
cat("Saving PCA grid...\n")
pca_plots <- lapply(results_all, function(x) x$pca)
pca_grid  <- wrap_plots(pca_plots, ncol = 3) +
  plot_annotation(
    title    = "PCA — Sample Overview (all 9 runs)",
    subtitle = grid_subtitle,
    theme    = theme(
      plot.title    = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 14, colour = "grey40")
    )
  )

ggsave(OUT_PCA, pca_grid,
       width = 27, height = 27, dpi = 300, limitsize = FALSE)
cat("  Saved:", OUT_PCA, "\n")

# ── Heatmap grid ──────────────────────────────────────────────────────────────
# wrap_elements() lets patchwork handle non-ggplot grobs (pheatmap gtables)
cat("Saving heatmap grid...\n")
hm_grobs  <- lapply(results_all, function(x) wrap_elements(x$heatmap))
hm_grid   <- wrap_plots(hm_grobs, ncol = 3) +
  plot_annotation(
    title    = "Sample Distance Heatmaps — all 9 runs",
    subtitle = grid_subtitle,
    theme    = theme(
      plot.title    = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 14, colour = "grey40")
    )
  )

ggsave(OUT_HEATMAP, hm_grid,
       width = 27, height = 27, dpi = 300, limitsize = FALSE)
cat("  Saved:", OUT_HEATMAP, "\n")

# ── Volcano grid ──────────────────────────────────────────────────────────────
# Each cell contains a 5-panel faceted volcano, so needs a larger canvas
cat("Saving volcano grid...\n")
vol_plots <- lapply(results_all, function(x) x$volcano)
vol_grid  <- wrap_plots(vol_plots, ncol = 3) +
  plot_annotation(
    title    = "Volcano Plots — Time course vs T=0 (all 9 runs)",
    subtitle = paste0(grid_subtitle,
                      "\nblue diamond = defence gene  |  red = sig. DE (padj<0.05, |LFC|>1)"),
    theme    = theme(
      plot.title    = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 13, colour = "grey40")
    )
  )

ggsave(OUT_VOLCANO, vol_grid,
       width = 42, height = 36, dpi = 300, limitsize = FALSE)
cat("  Saved:", OUT_VOLCANO, "\n")

# ── Open all three in the default image viewer ────────────────────────────────
cat("\nOpening output files...\n")
browseURL(OUT_PCA)
browseURL(OUT_HEATMAP)
browseURL(OUT_VOLCANO)

cat("\nDone! Files saved in:\n", getwd(), "\n")
