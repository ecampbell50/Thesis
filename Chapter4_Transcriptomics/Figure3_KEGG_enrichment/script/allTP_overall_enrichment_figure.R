# ══════════════════════════════════════════════════════════════════════════════
# allTP_overall_enrichment_figure.R
#
# WHAT THIS SCRIPT DOES
# ─────────────────────
# Builds one combined figure showing functional enrichment (COG / KEGG / GO) for
# every strain × phage combination, using the "all-timepoints unioned DEGs"
# enrichment results that were generated upstream by:
#   - COG-analysis-allTP.R
#   - KEGG-pathway-analysis-allTP.R
#   - GO-Enrichment-allTP.R
#
# Each of those upstream scripts:
#   1. Reads the DESeq2 result for each post-infection timepoint (T2,10,20,30,50)
#   2. Takes the UNION of DEGs (padj < 0.05 AND |log2FC| >= 1) across timepoints
#   3. Runs a SINGLE Fisher's exact test (hypergeometric) per functional term
#      asking: "are there more DEGs in this term than expected by chance?"
#   4. Writes results to a CSV (one row per term: pvalue, padj, gene Count, etc.)
#
# This script then loads those CSVs for all 6 conditions (3 strains × 2 phages)
# and renders them as three stacked dotplot panels.
#
# WHAT THE FIGURE SHOWS
# ─────────────────────
#   y-axis  = functional term (top N by significance)
#   x-axis  = condition (strain × phage)
#   colour  = p-value (Fisher's exact, gradient — redder = more significant)
#   size    = number of DEGs falling into that term ("Count")
#   shown   = ONLY terms passing FDR (Benjamini-Hochberg padj < 0.1) in at
#             least one condition.
# ══════════════════════════════════════════════════════════════════════════════

if (!requireNamespace("ggh4x", quietly = TRUE)) install.packages("ggh4x")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  library(scales)
  # ggh4x: per-strain strip background colours
  library(ggh4x)
})

# ── Paths & global thresholds ─────────────────────────────────────────────────

# Self-contained: resolve relative to this script's location in paper_package
# (…/S3_allTP_overall_enrichment/{script,data,figure}). Enrichment result CSVs in data/.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT <- normalizePath(file.path(.script_dir, ".."))
BASE <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "figure", "allTP_overall_enrichment_figure.png")

STRAINS <- c("C67", "D32", "D68")
PHAGES  <- c("BON", "CLY")

# Significance filter: Benjamini-Hochberg adjusted p (FDR) per condition.
PADJ_FILTER <- 0.1

# Maximum terms per panel, ranked by best padj across the 6 conditions
TOPN_COG  <- 15
TOPN_KEGG <- 15
TOPN_GO   <- 15

# x-axis tick label colours (BON / CLY)
PHAGE_COLOURS <- c(
  BON = "#000000",
  CLY = "#000000"
)

# Strain strip background colours
STRAIN_COLOURS <- c(
  C67 = "#aed1ea",
  D32 = "#93f0a3",
  D68 = "#f98ab3"
)

# ──────────────────────────────────────────────────────────────────────────────
# Step 1: LOAD enrichment results from disk
# ──────────────────────────────────────────────────────────────────────────────
# Result filenames differ by annotation type, hence the switch().
#
# GO_ALL_results.csv combines BP/MF/CC and can list the same Description under
# several IDs; collapse to one row per Description (min pvalue). The same
# deduplication is applied to all three term types.

load_results <- function(strain, phage, type) {
  dir <- file.path(BASE, strain, phage,
                   paste0(strain, "_", phage, "_allTP_", type, "_results"))
  f <- switch(type,
    COG  = file.path(dir, paste0(strain, "_", phage, "_allTP_COG_enrichment_results.csv")),
    KEGG = file.path(dir, "KEGG_allTP_enrichment_results.csv"),
    GO   = file.path(dir, "GO_ALL_results.csv")
  )
  if (!file.exists(f)) return(NULL)

  d <- read_csv(f, show_col_types = FALSE)
  if (nrow(d) == 0) return(NULL)

  d %>%
    # Make a tidy term label (trimmed + line-wrapped for the plot y-axis)
    mutate(
      Strain    = strain,
      Phage     = phage,
      Condition = paste0(strain, "\n", phage),
      term_raw  = str_trim(Description),
      term      = str_wrap(term_raw, width = 40)
    ) %>%
    # Keep only the columns we need downstream
    dplyr::select(term, term_raw, pvalue, p.adjust, Count, Strain, Phage, Condition) %>%
    # Deduplicate: one row per (term × condition), keep the strongest signal
    group_by(term, Condition) %>%
    slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()
}

# Pull all 6 conditions into one tidy long-format dataframe
collect_all <- function(type) {
  bind_rows(lapply(STRAINS, function(s) {
    lapply(PHAGES, function(p) load_results(s, p, type)) %>% bind_rows()
  }))
}

# ──────────────────────────────────────────────────────────────────────────────
# Step 2: BUILD a panel for one term type (COG / KEGG / GO)
# ──────────────────────────────────────────────────────────────────────────────

make_panel <- function(type, topN, title_label) {

  dat <- collect_all(type)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)

  # ── Term selection ────────────────────────────────────────────────────────
  # Keep terms with padj < PADJ_FILTER in at least one of the 6 conditions,
  # so condition-specific responses are retained.
  sig_terms <- dat %>%
    filter(!is.na(p.adjust), p.adjust < PADJ_FILTER) %>%
    pull(term) %>% unique()

  if (length(sig_terms) == 0) {
    message(type, ": no terms pass padj < ", PADJ_FILTER,
            " — nothing to plot for this panel.")
    return(NULL)
  }

  # Rank by best (minimum) padj across conditions and take the top N
  term_order <- dat %>%
    filter(term %in% sig_terms) %>%
    group_by(term) %>%
    summarise(best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_padj) %>%
    slice_head(n = topN) %>%
    pull(term)

  # ── Build the plotting dataframe ──────────────────────────────────────────
  # Full term × strain × phage grid; missing cells are NA and are not drawn.
  # x-axis = phage, faceted by strain.
  full_grid <- expand.grid(
    term   = term_order,
    Strain = STRAINS,
    Phage  = PHAGES,
    stringsAsFactors = FALSE
  )

  dat_plot <- full_grid %>%
    left_join(
      dat %>% filter(term %in% term_order),
      by = c("term", "Strain", "Phage")
    ) %>%
    mutate(
      term   = factor(term, levels = rev(term_order)),    # rev = top of plot
      Strain = factor(Strain, levels = STRAINS),
      Phage  = factor(Phage,  levels = PHAGES),
      # -log10(p) for colour; +1e-300 avoids log(0), capped at 15
      neg_log_p = pmin(-log10(pvalue + 1e-300), 15)
    )

  # Only draw cells passing FDR
  d_draw <- dat_plot %>%
    filter(!is.na(p.adjust), p.adjust < PADJ_FILTER)

  # Colour legend labelled at p = 0.05, 0.01, 0.001, 1e-5, 1e-10
  p_breaks <- c(0.05, 0.01, 0.001, 1e-5, 1e-10)
  p_breaks <- p_breaks[p_breaks >= 10^(-15)]      # within capped range
  break_pos <- -log10(p_breaks)
  break_lbl <- c("0.05", "0.01", "0.001",
                 expression(10^-5), expression(10^-10))[seq_along(p_breaks)]

  x_text_cols <- PHAGE_COLOURS[PHAGES]

  # Strip colours only for strains present in this panel
  present_strains <- levels(droplevels(d_draw$Strain))
  strip_rects <- lapply(present_strains, function(s) {
    element_rect(fill = STRAIN_COLOURS[s], colour = NA)
  })

  ggplot(d_draw, aes(x = Phage, y = term)) +
    geom_point(
      aes(size = Count, fill = neg_log_p),
      shape = 21, colour = "black", stroke = 0.5, alpha = 0.95
    ) +

    facet_grid2(. ~ Strain,
                strip = strip_themed(background_x = strip_rects)) +

    scale_fill_gradientn(
      colours  = c("#4575B4", "#74ADD1", "#FEE090", "#F46D43", "#A50026"),
      limits   = c(-log10(PADJ_FILTER), 15),    # start at the FDR threshold
      breaks   = break_pos,
      labels   = break_lbl,
      name     = "Fisher's\nexact p",
      oob      = scales::squish,
      na.value = NA
    ) +

    scale_size_continuous(
      name   = "DEGs in\nterm",
      range  = c(2, 12),
      breaks = c(2, 5, 10, 25, 50, 100)
    ) +

    scale_x_discrete(drop = FALSE) +
    labs(title = title_label, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major   = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size = 9),
      axis.text.x        = element_text(size = 9.5, face = "bold",
                                        colour = x_text_cols),
      strip.text.x       = element_text(size = 10, face = "bold",
                                        margin = margin(t = 3, b = 3)),
      panel.spacing.x    = unit(0.6, "lines"),
      plot.title         = element_text(size = 11, face = "bold",
                                        margin = margin(b = 4)),
      legend.position    = "right",
      legend.key.size    = unit(0.5, "cm"),
      legend.text        = element_text(size = 8),
      legend.title       = element_text(size = 8.5),
      plot.margin        = margin(4, 4, 4, 4)
    ) +
    guides(
      fill = guide_colourbar(order = 1, barwidth = 0.7, barheight = 5,
                             title.position = "top", reverse = TRUE),
      size = guide_legend(order = 2,
                          override.aes = list(shape = 21,
                                              fill = "grey70",
                                              colour = "black"))
    )
}

# ──────────────────────────────────────────────────────────────────────────────
# Step 3: ASSEMBLE three panels into the final figure
# ──────────────────────────────────────────────────────────────────────────────

p_cog  <- make_panel("COG",  TOPN_COG,  "A  COG functional categories")
p_kegg <- make_panel("KEGG", TOPN_KEGG, "B  KEGG pathways")
p_go   <- make_panel("GO",   TOPN_GO,   "C  Gene Ontology (BP/MF/CC combined)")

panels <- Filter(Negate(is.null), list(p_cog, p_kegg, p_go))
if (length(panels) == 0) stop("No panels produced.")

combined <- wrap_plots(panels, ncol = 1) +
  plot_annotation(
    title    = "Functional enrichment across all post-infection timepoints (T2–T50)",
    subtitle = paste0("Union of DEGs (|log2FC| >= 1, padj < 0.05) per strain x phage. ",
                      "Only terms with FDR < ", PADJ_FILTER, " shown."),
    theme = theme(
      plot.title    = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  )

# Dynamic figure height — keep dots from getting squashed when many terms shown
total_rows <- sum(sapply(panels, function(p) length(levels(p$data$term))))
fig_h <- max(12, total_rows * 0.40 + 5)

ggsave(OUT, combined, width = 11, height = fig_h, dpi = 300, bg = "white")
cat("Saved:", OUT, "\n")
