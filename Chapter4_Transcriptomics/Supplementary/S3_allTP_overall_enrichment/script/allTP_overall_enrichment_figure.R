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
#   shown   = ONLY terms passing FDR (Benjamini-Hochberg padj < 0.05) in at
#             least one condition. Less ambiguous than the previous "nominal p"
#             filter — every dot you see is statistically defensible.
# ══════════════════════════════════════════════════════════════════════════════

install.packages("ggh4x")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  library(scales)
  # ggh4x extends ggplot2's faceting — we use it to colour each strain's
  # strip background independently (vanilla ggplot2 only supports one colour
  # for all strips). Install once with: install.packages("ggh4x")
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

# Significance filter applied to the figure.
# We filter on Benjamini-Hochberg adjusted p (FDR), NOT nominal pvalue.
# Rationale: with hundreds of terms tested per condition, a fair chunk would
# hit nominal p<0.1 by chance. FDR adjustment controls the expected fraction
# of false positives. Anything shown here passes a defensible threshold.
PADJ_FILTER <- 0.1

# Maximum number of terms to display per panel (ranked by best padj across all
# 6 conditions). Lower this if the figure gets cramped.
TOPN_COG  <- 15
TOPN_KEGG <- 15
TOPN_GO   <- 15

# ── Phage colour overlay ──────────────────────────────────────────────────────
# These hex codes are used to colour the x-axis tick labels (BON / CLY) inside
# each strain sub-panel. Edit these to match the phage colours used in your
# other figures so the whole thesis chapter is visually consistent.
PHAGE_COLOURS <- c(
  BON = "#000000",   # ← edit to your Bonnie colour
  CLY = "#000000"    # ← edit to your Clyde colour
)

# Strain strip background colours — these now drive the coloured header above
# each strain's sub-panel (requires ggh4x). Edit the hex codes to match your
# strain colour scheme from other figures.
STRAIN_COLOURS <- c(
  C67 = "#aed1ea",
  D32 = "#93f0a3",
  D68 = "#f98ab3"
)

# ──────────────────────────────────────────────────────────────────────────────
# Step 1: LOAD enrichment results from disk
# ──────────────────────────────────────────────────────────────────────────────
# For each strain × phage, point to the right CSV. The upstream scripts write
# them to different filenames (no consistency), hence the switch().
#
# DEDUPLICATION: GO_ALL_results.csv contains entries from BP/MF/CC ontologies
# and sometimes lists the same Description multiple times under different IDs.
# Without dedup, those plot as overlapping circles at the same y-coordinate
# (e.g. "structural constituent of ribosome" appears 5× in some files).
# We collapse to one row per Description by keeping the most significant (min
# pvalue). Same dedup is applied to all three term types for consistency.

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
  # A term is shown if it is FDR-significant (padj < PADJ_FILTER) in at least
  # ONE of the 6 conditions. We don't require every condition to be FDR-sig —
  # that would hide condition-specific responses (a term enriched only in
  # C67-BON is interesting even if it's not enriched anywhere else).
  sig_terms <- dat %>%
    filter(!is.na(p.adjust), p.adjust < PADJ_FILTER) %>%
    pull(term) %>% unique()

  if (length(sig_terms) == 0) {
    message(type, ": no terms pass padj < ", PADJ_FILTER,
            " — nothing to plot for this panel.")
    return(NULL)
  }

  # Of those, rank by BEST padj across conditions and take the top N.
  # Using min(padj) makes the ranking reflect "strongest signal anywhere",
  # which is usually what you want for a summary figure.
  term_order <- dat %>%
    filter(term %in% sig_terms) %>%
    group_by(term) %>%
    summarise(best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_padj) %>%
    slice_head(n = topN) %>%
    pull(term)

  # ── Build the plotting dataframe ──────────────────────────────────────────
  # We want EVERY (term × strain × phage) combination represented — even those
  # that didn't reach significance — so the dotplot reads as a grid.
  # Missing cells are NA and won't draw a dot.
  #
  # New layout: x-axis is now Phage (BON/CLY) and we facet horizontally by
  # Strain (C67/D32/D68). This makes within-strain phage comparisons easier.
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
      # Compute -log10(p) for the colour scale. The +1e-300 is to avoid log(0)
      # if any pvalue is exactly zero. The pmin() cap prevents one extreme
      # outlier from dominating the colour scale.
      neg_log_p = pmin(-log10(pvalue + 1e-300), 15)
    )

  # ── Layer: only draw dots for cells passing FDR ───────────────────────────
  # Cells that didn't pass FDR are NA — they leave empty space in the grid.
  d_draw <- dat_plot %>%
    filter(!is.na(p.adjust), p.adjust < PADJ_FILTER)

  # ── Build the plot ────────────────────────────────────────────────────────
  # Colour scale: -log10(p), but we relabel the legend at meaningful p-value
  # thresholds (p=0.05, 0.01, 0.001, 1e-5, 1e-10) so the colour is interpretable
  # without doing log math in your head.
  p_breaks <- c(0.05, 0.01, 0.001, 1e-5, 1e-10)
  p_breaks <- p_breaks[p_breaks >= 10^(-15)]      # within capped range
  break_pos <- -log10(p_breaks)
  break_lbl <- c("0.05", "0.01", "0.001",
                 expression(10^-5), expression(10^-10))[seq_along(p_breaks)]

  # X-axis text colours: order matches the PHAGES vector (BON then CLY).
  # ggplot applies the vector left-to-right across the discrete x-axis ticks.
  x_text_cols <- PHAGE_COLOURS[PHAGES]

  # Build strip colour list dynamically from whichever strains actually have
  # FDR-sig terms in this panel. Positional assignment (always 3 rects) would
  # give the wrong colour when fewer than 3 strains are present (e.g. GO panel).
  present_strains <- levels(droplevels(d_draw$Strain))
  strip_rects <- lapply(present_strains, function(s) {
    element_rect(fill = STRAIN_COLOURS[s], colour = NA)
  })

  ggplot(d_draw, aes(x = Phage, y = term)) +
    # Single layer, consistent style: all dots are FDR-significant by definition.
    geom_point(
      aes(size = Count, fill = neg_log_p),
      shape = 21, colour = "black", stroke = 0.5, alpha = 0.95
    ) +

    # ── Facet horizontally by strain ──────────────────────────────────────
    facet_grid2(. ~ Strain,
                strip = strip_themed(background_x = strip_rects)) +

    # Colour: -log10(p), labelled with actual p-values
    scale_fill_gradientn(
      colours  = c("#4575B4", "#74ADD1", "#FEE090", "#F46D43", "#A50026"),
      limits   = c(-log10(PADJ_FILTER), 15),    # start at the FDR threshold
      breaks   = break_pos,
      labels   = break_lbl,
      name     = "Fisher's\nexact p",
      oob      = scales::squish,
      na.value = NA
    ) +

    # Size: wider range for better visual discrimination than before.
    # range = (2, 12) means smallest dot = 2pt, largest = 12pt.
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
      # X-axis labels coloured per phage (BON / CLY) — edit PHAGE_COLOURS above
      axis.text.x        = element_text(size = 9.5, face = "bold",
                                        colour = x_text_cols),
      # Strain facet header text — background fill comes from strip_themed()
      # in the facet_grid2 call above, so we only style the text here
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
