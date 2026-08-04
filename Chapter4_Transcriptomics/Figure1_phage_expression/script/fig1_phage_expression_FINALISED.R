# fig1_phage_expression.R  ── Figure 1: Phage gene expression overview
#
# Panel A  — CSR ternary: Early / Middle / Late gene class proportions
#             for each strain × phage combination
# Panels B–G — Joy / ridgeline plots of phage gene expression across the
#               infection timecourse (one panel per strain × phage, genomic order)
#
# Layout (2 rows × 3 cols):
#   Row 1: C67 Bonnie  |  D32 Bonnie  | D68 Clyde
#   Row 2: C67 Clyde   |  D32 Clyde   | D68 Clyde
#
# ── RUN FROM ──────────────────────────────────────────────────────────────────
#   setwd("...2_Analysis/")   ← project root
#   source("scripts/R/fig1_phage_expression.R")
#
# ── DATA IN ───────────────────────────────────────────────────────────────────
# Panel A (ternary):
#   · Gene class counts HARDCODED below.
#     Source: PhageExpressionAtlas ClassThreshold column in
#             results/PEA/{STRAIN}_{PHAGE}_fractional_expression.tsv
#     Totals: Bonnie = 72 genes, Clyde = 65 genes
#
# Panels B–G (ridgelines):
#   · data/raw_counts/{STRAIN}_{PHAGE}_full_raw_counts.tsv
#     18 count columns: {0,2,10,20,30,50}_{STRAIN}_{PHAGE}_R{1,2,3}
#     + Entity and Symbol columns (dropped before analysis)
#   · results/PEA/{STRAIN}_{PHAGE}_fractional_expression.tsv
#     Column used: ClassThreshold  (Early / Middle / Late / None)
#
# ── FILTERING ─────────────────────────────────────────────────────────────────
#   · T=0 columns dropped  (pre-infection; phage not yet injected)
#   · Only phage gene rows kept  (row name prefix: PENJXGPI_CDS_ or ZUDWSPYW_CDS_)
#   · Genes with rowSums(counts) < 10 excluded before VST
#
# ── STATS ─────────────────────────────────────────────────────────────────────
#   · VST: DESeq2::varianceStabilizingTransformation(blind = TRUE)
#           applied independently per strain × phage dataset
#   · Replicates averaged: mean(R1, R2, R3) per timepoint after VST
#   · Per-gene scaling: (x - min) / (max - min + 1e-9)  → 0–1 waveform height
#   · Smooth curves: cubic spline with 300 output points
#   · Ternary: pct = count / total × 100  (no statistical test)
#
# ── OUTPUT ────────────────────────────────────────────────────────────────────
#   results/figures/fig1_phage_expression.pdf   (10 × 16 in)
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Packages ───────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("DESeq2"))
  if (!p %in% installed.packages()) BiocManager::install(p, update = FALSE)
for (p in c("ggtern", "ggridges", "tidyverse", "cowplot", "magick"))
  if (!p %in% installed.packages()) install.packages(p)

library(DESeq2);   library(ggtern);  library(ggridges)
library(tidyverse); library(cowplot); library(magick)

# ── 1. Paths ──────────────────────────────────────────────────────────────────
# Self-contained: resolve relative to this script's location in paper_package
# (…/Figure1_phage_expression/{script,data,figure}). No external paths.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT         <- normalizePath(file.path(.script_dir, ".."))
COUNTS_DIR   <- file.path(ROOT, "data", "raw_counts")
PEA_DIR      <- file.path(ROOT, "data", "PEA")
OUT_DIR      <- file.path(ROOT, "figure")
CLASS_METHOD <- "ClassThreshold"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# PANEL A — TERNARY
# ══════════════════════════════════════════════════════════════════════════════
# Gene class counts from PhageExpressionAtlas (ClassThreshold).
# Taken from results/PEA/{STRAIN}_{PHAGE}_fractional_expression.tsv
phage_data <- tribble(
  ~label,          ~Phage,   ~Strain,  ~Early, ~Middle, ~Late, ~Total,
  "C67 · Bonnie",  "Bonnie", "C67",       37,      26,     9,     72,
  "D32 · Bonnie",  "Bonnie", "D32",       28,      44,     0,     72,
  "D68 · Bonnie",  "Bonnie", "D68",       35,      31,     6,     72,
  "C67 · Clyde",   "Clyde",  "C67",       35,       4,    26,     65,
  "D32 · Clyde",   "Clyde",  "D32",       33,       4,    28,     65,
  "D68 · Clyde",   "Clyde",  "D68",       16,      30,    19,     65
) %>%
  mutate(pE = Early  / Total * 100,
         pM = Middle / Total * 100,
         pL = Late   / Total * 100)

cat("── TERNARY INPUT DATA ────────────────────────────────────────────────────\n")
print(phage_data %>% select(label, Early, Middle, Late, Total, pE, pM, pL))

strain_cols  <- c(C67 = "#FCE800", D32 = "#8AE788", D68 = "#F52DAD")
phage_shapes <- c(Bonnie = 21, Clyde = 24)

ternary_plt <- ggtern(
  phage_data,
  aes(x = pM, y = pE, z = pL, fill = Strain, shape = Phage)
) +
  geom_point(size = 5, stroke = 1.2, colour = "white") +
  geom_text(aes(label = label, colour = Strain),
            size = 3, fontface = "bold",
            vjust = -1.0, hjust = 0.5, lineheight = 0.9) +
  Tlab("Early") + Llab("Middle") + Rlab("Late") +
  scale_colour_manual(values = strain_cols, guide = "none") +
  scale_fill_manual(values   = strain_cols, name = "Strain") +
  scale_shape_manual(values  = phage_shapes, name = "Phage") +
  labs(title    = "CSR representation of S. suis phages",
       subtitle = "Early = Competitor  ·  Middle = Stress-tolerator  ·  Late = Ruderal") +
  theme_bw(base_size = 11) + theme_showarrows() +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5),
        plot.subtitle   = element_text(size = 8, colour = "grey40", hjust = 0.5),
        legend.position = "right")

# render to a temp PNG, reload as a raster image for cowplot.
tmp_tern <- tempfile(fileext = ".png")
ggsave(tmp_tern, ternary_plt, width = 7, height = 6, dpi = 200, bg = "white")
tern_panel <- ggdraw() + draw_image(tmp_tern)

# ══════════════════════════════════════════════════════════════════════════════
# PANELS B–G — RIDGELINES
# Colour scheme and white background from phage_ridgeline_indv.R
# ══════════════════════════════════════════════════════════════════════════════
class_cols <- c(Early   = "#99ACFF",   # soft blue
                Middle  = "#F4D06F",   # soft yellow
                Late    = "#9A4C95",   # purple
                Unknown = "grey70")

# Loop order: fill row-by-row so plot_grid produces 2 rows × 3 cols
# Row 1 (BON): C67 | D32 | D68
# Row 2 (CLY): C67 | D32 | D68
COMBOS <- data.frame(
  STRAIN = rep(c("C67", "D32", "D68"), times = 2),
  PHAGE  = rep(c("BON", "CLY"), each = 3),
  stringsAsFactors = FALSE
)

ridge_plots <- lapply(seq_len(nrow(COMBOS)), function(i) {

  STRAIN  <- COMBOS$STRAIN[i]
  PHAGE   <- COMBOS$PHAGE[i]
  pname   <- if (PHAGE == "BON") "Bonnie" else "Clyde"
  pfx     <- if (PHAGE == "BON") "PENJXGPI_CDS_" else "ZUDWSPYW_CDS_"

  # ── Read raw counts ────────────────────────────────────────────────────────
  f   <- file.path(COUNTS_DIR, paste0(STRAIN, "_", PHAGE, "_full_raw_counts.tsv"))
  cat("\nReading:", f, "\n")
  raw <- read.delim(f, row.names = 1, check.names = FALSE)
  mat <- raw[, !colnames(raw) %in% c("Entity", "Symbol")]

  # ── Filter ─────────────────────────────────────────────────────────────────
  mat <- mat[grepl(pfx, rownames(mat)),    ]   # phage genes only
  mat <- mat[, !grepl("^0_", colnames(mat))]   # drop T=0
  mat <- mat[rowSums(mat) >= 10,           ]   # low-count filter
  cat(sprintf("  Phage genes retained: %d\n", nrow(mat)))

  # ── VST ────────────────────────────────────────────────────────────────────
  cd  <- data.frame(row.names = colnames(mat),
                    time = factor(sub("_R[0-9]+$", "", colnames(mat))))
  dds <- DESeqDataSetFromMatrix(mat, cd, ~ time)
  vst <- assay(varianceStabilizingTransformation(dds, blind = TRUE))

  # Replicate average → ordered timepoints
  tp_labs  <- sub("_R[0-9]+$", "", colnames(vst))
  vst_mean <- sapply(unique(tp_labs), function(tp)
    rowMeans(vst[, tp_labs == tp, drop = FALSE]))
  vst_mean <- vst_mean[, order(as.numeric(colnames(vst_mean))), drop = FALSE]

  # Per-gene min-max scale → 0–1
  vst_sc <- t(apply(vst_mean, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-9)))

  # ── Gene class from PEA ────────────────────────────────────────────────────
  pea <- read.delim(
    file.path(PEA_DIR, paste0(STRAIN, "_", PHAGE, "_fractional_expression.tsv")),
    row.names = 1, check.names = FALSE
  )
  class_df <- data.frame(
    gene_id = rownames(pea)[grepl(pfx, rownames(pea))],
    class   = pea[grepl(pfx, rownames(pea)), CLASS_METHOD]
  ) %>%
    mutate(class = factor(
      ifelse(class %in% c("Early", "Middle", "Late"), class, "Unknown"),
      levels = c("Early", "Middle", "Late", "Unknown")
    ))

  # ── Spline interpolation ───────────────────────────────────────────────────
  tps    <- as.numeric(colnames(vst_sc))
  tp_seq <- seq(min(tps), max(tps), length.out = 300)

  df_sm <- vst_sc %>%
    as.data.frame() %>% rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "timepoint", values_to = "expr") %>%
    mutate(timepoint = as.numeric(timepoint)) %>%
    group_by(gene_id) %>%
    group_modify(~ {
      s <- spline(.x$timepoint, .x$expr, xout = tp_seq)
      data.frame(timepoint = s$x, expr = pmax(s$y, 0))
    }) %>%
    ungroup() %>%
    left_join(class_df, by = "gene_id")

  # Genomic order; fct_rev → CDS_0001 at top of plot
  gene_ord  <- paste0(pfx, sprintf("%04d",
    sort(as.integer(sub(pfx, "", unique(df_sm$gene_id))))))
  df_sm$gene_id <- factor(df_sm$gene_id, levels = rev(gene_ord))

  # ── Plot (white background, grey text) ────-----------------------------------
  ggplot(df_sm, aes(x = timepoint, y = gene_id,
                    height = expr, fill = class, colour = class)) +
    geom_ridgeline(alpha = 0.75, linewidth = 0.35, scale = 3, min_height = 0) +
    scale_fill_manual(values   = class_cols, name = CLASS_METHOD, drop = FALSE) +
    scale_colour_manual(values = class_cols, drop = FALSE, guide = "none") +
    scale_x_continuous(breaks = tps, labels = paste0(tps, " min"), expand = c(0.01, 0)) +
    labs(title = paste0(STRAIN, " · ", pname), x = "Time post-infection", y = NULL) +
    theme_void(base_size = 9) +
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      axis.text.x       = element_text(colour = "grey40", size = 7, margin = margin(t = 4)),
      axis.text.y       = element_text(colour = "grey50", size = 4, hjust = 1, margin = margin(r = 2)),
      axis.ticks.x      = element_line(colour = "grey40", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      plot.title        = element_text(colour = "black", size = 10, face = "bold",
                                        hjust = 0.5, margin = margin(b = 4)),
      legend.position   = "none",
      plot.margin       = margin(10, 10, 10, 10)
    )
})
names(ridge_plots) <- paste0(COMBOS$STRAIN, "_", COMBOS$PHAGE)

# ── Verify gene class counts match ternary hardcoded values ───────────────────
cat("\n── GENE CLASS COUNTS IN PEA FILES (check vs ternary hardcoded values) ──\n")
for (nm in names(ridge_plots)) {
  parts  <- strsplit(nm, "_")[[1]]
  S <- parts[1]; P <- parts[2]
  pfx <- if (P == "BON") "PENJXGPI_CDS_" else "ZUDWSPYW_CDS_"
  pea <- read.delim(
    file.path(PEA_DIR, paste0(S, "_", P, "_fractional_expression.tsv")),
    row.names = 1, check.names = FALSE
  )
  cls <- pea[grepl(pfx, rownames(pea)), CLASS_METHOD]
  cat(sprintf("  %-12s  %s\n", nm,
              paste(names(table(cls)), table(cls), sep = "=", collapse = "  ")))
}

# ══════════════════════════════════════════════════════════════════════════════
# COMBINE AND SAVE
# ══════════════════════════════════════════════════════════════════════════════
# Ridge grid: 3 rows × 2 cols
# COMBOS is ordered (C67-BON, C67-CLY, D32-BON, D32-CLY, D68-BON, D68-CLY)
ridge_grid <- plot_grid(
  plotlist    = ridge_plots,
  nrow        = 2,
  ncol        = 3,
  labels      = LETTERS[2:7],
  label_size  = 12
)

full_fig <- plot_grid(
  tern_panel, ridge_grid,
  nrow        = 2,
  rel_heights = c(1, 2.2),
  labels      = c("A", ""),
  label_size  = 14
)

out_file <- file.path(OUT_DIR, "fig1_phage_expression.pdf")
ggsave(out_file, full_fig, width = 10, height = 16, dpi = 300)
cat("\nSaved:", out_file, "\n")
browseURL(out_file)
