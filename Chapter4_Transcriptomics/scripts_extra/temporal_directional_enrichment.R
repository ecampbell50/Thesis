# COG and KEGG terms: % of genes up/down vs control over time (Fig S4.7, S4.8).

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler", update = FALSE)
for (pkg in c("ggh4x", "patchwork")) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

suppressPackageStartupMessages({
  library(clusterProfiler); library(dplyr); library(readr); library(tidyr)
  library(stringr); library(ggplot2); library(patchwork); library(ggh4x); library(scales)
})

# ── 1. Paths & settings ───────────────────────────────────────────────────────
ANALYSIS  <- "/Users/emmet/Desktop/Thesis/Project_Chp3_Transcriptomics/2_Analysis"
DESEQ_DIR <- file.path(ANALYSIS, "results", "DESeq2_vs_control")  # vs-control CSVs
ANNO_DIR  <- file.path(ANALYSIS, "Maria_Pipeline")
OUT_DIR   <- file.path(ANALYSIS, "results", "temporal_enrichment_vsCTLR")
FIG_OUT   <- file.path(OUT_DIR, "temporal_directional_enrichment.png")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STRAINS     <- c("C67", "D32", "D68")
PHAGES      <- c("BON", "CLY")
PHAGE_LABEL <- c(BON = "Bonnie", CLY = "Clyde")
TIMES       <- c(2, 10, 20, 30, 50)

SIG_PADJ <- 0.05   # "perturbed" = padj < this ...
SIG_LFC  <- 1      # ... AND |log2FC| >= this (same as your DEG definition)

PADJ_FILTER      <- 0.1   # term shown only if union-ORA FDR < this in >=1 condition
TOPN             <- c(COG = 12, KEGG = 15)
MIN_PATHWAY_SIZE <- 1     # skip terms with < this many tested genes (noisy %)

ROW_HALF <- 0.42                       # max half-band height (row units)
STRAIN_COLOURS <- c(C67 = "#aed1ea", D32 = "#93f0a3", D68 = "#f98ab3")

# ── Ridge colours ─────────────────────────────────────────────────────────────
# COLOUR_BY = "direction" : simple red (induced) / blue (repressed), one legend.
# COLOUR_BY = "phage"     : per strain × phage colour families — keeps each
#   strain in one hue but separates Bonnie vs Clyde by shade. Within a phage,
#   the DARKER shade = induced (up), LIGHTER = repressed (down). Edit freely.
COLOUR_BY <- "phage"     # "phage" or "direction"

DIR_COLS <- c(up = "#D73027", down = "#4575B4")   # used when COLOUR_BY = "direction"

# c(up = darker, down = lighter) per strain × phage. C67 = two blues, etc.
PHAGE_COLS <- list(
  C67 = list(BON = c(up = "#08306b", down = "#c6dbef"),   # navy / light blue
             CLY = c(up = "#08306b", down = "#c6dbef")),   # medium / pale blue
  D32 = list(BON = c(up = "#586B24", down = "#C0D684"),    # dark / light green
             CLY = c(up = "#586B24", down = "#C0D684")),
  D68 = list(BON = c(up = "#3D0B37", down = "#DDA6C6"),    # dark / light pink
             CLY = c(up = "#3D0B37", down = "#DDA6C6")))

# c(up = darker, down = lighter) per strain × phage. C67 = two blues, etc.
#PHAGE_COLS <- list(
#  C67 = list(BON = c(up = "#08306b", down = "#9ecae1"),   # navy / light blue
#             CLY = c(up = "#2171b5", down = "#c6dbef")),   # medium / pale blue
#  D32 = list(BON = c(up = "#00441b", down = "#a1d99b"),    # dark / light green
#             CLY = c(up = "#41ab5d", down = "#d9f0d3")),
#  D68 = list(BON = c(up = "#7a0177", down = "#fbb4b9"),    # dark / light pink
#             CLY = c(up = "#dd3497", down = "#fde0dd")))

# ── 2. Helpers ────────────────────────────────────────────────────────────────
clean_term_string <- function(x) { x <- str_trim(x); x[x == ""] <- NA_character_; x }

# TERM2GENE (term name → gene) for COG or KEGG, from the eggNOG annotation file.
build_terms <- function(anno, type) {
  col <- if (type == "COG") "COG_names" else "KEGG_Pathway_names"
  anno %>% dplyr::select(GeneID, !!sym(col)) %>%
    rename(val = !!sym(col)) %>%
    filter(!is.na(val), val != "") %>%
    separate_rows(val, sep = "; ") %>%
    mutate(val = clean_term_string(val)) %>%
    filter(!is.na(val)) %>%
    transmute(term = val, gene = GeneID) %>% distinct()
}

# Read all vs-control timepoints for one strain × phage into a long table.
read_vsctrl <- function(strain, phage) {
  rows <- list()
  for (tp in TIMES) {
    f <- file.path(DESEQ_DIR, sprintf("DESeq2_%s_%s_vs_CTLR_T%dmin.csv", strain, phage, tp))
    if (!file.exists(f)) { message("  [missing] ", basename(f)); next }
    d <- read_csv(f, show_col_types = FALSE); colnames(d)[1] <- "gene"
    rows[[as.character(tp)]] <- d %>%
      filter(!str_starts(gene, "PEN"), !str_starts(gene, "ZUD")) %>%
      transmute(gene, time = tp, log2FoldChange, padj)
  }
  bind_rows(rows)
}

# ── 3. Compute per-term, per-timepoint directional % (+ select enriched terms) ─
cat("── Computing temporal directional enrichment (vs control) ────────────────\n")
pert_all <- list(); sel_terms <- list()

for (type in c("COG", "KEGG")) {
  enriched <- list(); pert_rows <- list()

  for (strain in STRAINS) {
    af <- file.path(ANNO_DIR, sprintf("functional_annotation_%s_with_names.csv", strain))
    if (!file.exists(af)) stop("Missing annotation: ", af)
    anno <- read_csv(af, show_col_types = FALSE) %>%
      filter(!str_starts(GeneID, "PEN"), !str_starts(GeneID, "ZUD"))
    t2g <- build_terms(anno, type)
    if (nrow(t2g) == 0) next
    membership <- split(t2g$gene, t2g$term)

    for (phage in PHAGES) {
      res <- read_vsctrl(strain, phage)
      if (nrow(res) == 0) next
      tested <- unique(res$gene)
      bg     <- intersect(tested, unique(t2g$gene))

      # (a) ORA on the union DEG set → which terms are worth showing
      sig <- res %>%
        filter(!is.na(padj), padj < SIG_PADJ,
               !is.na(log2FoldChange), abs(log2FoldChange) >= SIG_LFC) %>%
        pull(gene) %>% unique() %>% intersect(bg)
      if (length(sig) > 0) {
        e <- enricher(sig, universe = bg, TERM2GENE = t2g,
                      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
        if (!is.null(e)) {
          edf <- as.data.frame(e) %>% filter(!is.na(p.adjust), p.adjust < PADJ_FILTER)
          if (nrow(edf) > 0) enriched[[paste(strain, phage)]] <- edf %>% dplyr::select(ID, p.adjust)
        }
      }

      # (b) per-term, per-timepoint directional % of tested genes
      for (term in names(membership)) {
        pg <- intersect(membership[[term]], tested)
        n_tested <- length(pg)
        if (n_tested < MIN_PATHWAY_SIZE) next
        sub <- res %>% filter(gene %in% pg)
        for (tp in TIMES) {
          st <- sub %>% filter(time == tp, !is.na(padj), !is.na(log2FoldChange))
          up <- sum(st$padj < SIG_PADJ & st$log2FoldChange >=  SIG_LFC)
          dn <- sum(st$padj < SIG_PADJ & st$log2FoldChange <= -SIG_LFC)
          pert_rows[[paste(strain, phage, term, tp)]] <- tibble(
            term_type = type, strain = strain, phage = phage, term = term, time = tp,
            n_tested = n_tested, up_pct = 100 * up / n_tested, down_pct = 100 * dn / n_tested)
        }
      }
    }
  }

  pert <- bind_rows(pert_rows)
  pert_all[[type]] <- pert

  # Common term set across strains: enriched terms ranked by best (min) FDR.
  if (length(enriched) == 0 || nrow(pert) == 0) { sel_terms[[type]] <- character(0); next }
  term_padj <- bind_rows(enriched) %>% group_by(ID) %>%
    summarise(best = min(p.adjust), .groups = "drop")
  ord <- term_padj %>% filter(ID %in% pert$term) %>%
    arrange(best) %>% slice_head(n = TOPN[[type]]) %>% pull(ID)
  sel_terms[[type]] <- ord
  write_csv(pert %>% filter(term %in% ord),
            file.path(OUT_DIR, paste0("temporal_", type, "_directional_table.csv")))
  cat(sprintf("  %s: %d terms selected\n", type, length(ord)))
}

# ── 4. Build one ridgeline panel for a term type ──────────────────────────────
make_ridge <- function(type, title) {
  terms <- sel_terms[[type]]
  if (length(terms) == 0) { message(type, ": no terms — skipping."); return(NULL) }
  pert <- pert_all[[type]] %>% filter(term %in% terms)

  # Gene-count per term (mean tested genes across strains) for the y-axis labels.
  term_n <- pert %>% group_by(term) %>%
    summarise(n = round(mean(n_tested)), .groups = "drop")
  lab_map <- setNames(sprintf("%s  (n≈%d)", term_n$term, term_n$n), term_n$term)

  pert$term <- factor(pert$term, levels = rev(terms))   # first selected = top row
  pert$y    <- as.integer(pert$term)

  maxpct <- max(c(pert$up_pct, pert$down_pct), na.rm = TRUE); if (maxpct <= 0) maxpct <- 1
  sc <- ROW_HALF / maxpct
  pert <- pert %>% mutate(up_h = up_pct * sc, dn_h = down_pct * sc)

  # Long ribbon table; keep raw strain/phage/dir codes for colour keying.
  rib <- bind_rows(
    pert %>% transmute(strain, phage, term, time, y, dir = "up",
                       ymin = y, ymax = y + up_h, edge = y + up_h),
    pert %>% transmute(strain, phage, term, time, y, dir = "down",
                       ymin = y - dn_h, ymax = y, edge = y - dn_h)
  )

  # Colour key + values depend on COLOUR_BY.
  if (COLOUR_BY == "phage") {
    rib$fill_key <- paste(rib$strain, rib$phage, rib$dir, sep = "|")
    fill_values <- unlist(lapply(STRAINS, function(s) lapply(PHAGES, function(p)
      setNames(c(PHAGE_COLS[[s]][[p]][["up"]], PHAGE_COLS[[s]][[p]][["down"]]),
               c(paste(s, p, "up", sep = "|"), paste(s, p, "down", sep = "|"))))))
    fill_values <- unlist(fill_values); fill_guide <- "none"
  } else {
    rib$fill_key <- ifelse(rib$dir == "up", "Induced (up vs ctrl)", "Repressed (down vs ctrl)")
    fill_values <- setNames(c(DIR_COLS[["up"]], DIR_COLS[["down"]]),
                            c("Induced (up vs ctrl)", "Repressed (down vs ctrl)"))
    fill_guide <- guide_legend(title = NULL)
  }

  rib <- rib %>%
    mutate(strain = factor(strain, levels = STRAINS),
           phage  = factor(PHAGE_LABEL[phage], levels = unname(PHAGE_LABEL[PHAGES]))) %>%
    arrange(time)

  ylabs <- levels(pert$term)
  top_y <- length(ylabs)

  # ±% scale key + sign labels, drawn only in the first facet (top row).
  key_pct <- pretty(c(0, maxpct), n = 2); key_pct <- key_pct[key_pct > 0 & key_pct <= maxpct]
  key_df <- data.frame(
    strain = factor(STRAINS[1], levels = STRAINS),
    phage  = factor(PHAGE_LABEL[[PHAGES[1]]], levels = unname(PHAGE_LABEL[PHAGES])),
    pct = key_pct, y0 = top_y)
  sign_df <- data.frame(
    strain = factor(STRAINS[1], levels = STRAINS),
    phage  = factor(PHAGE_LABEL[[PHAGES[1]]], levels = unname(PHAGE_LABEL[PHAGES])),
    lab = c("+% induced", "−% repressed"),
    y   = c(top_y + ROW_HALF * 0.55, top_y - ROW_HALF * 0.55))

  ggplot(rib, aes(x = time)) +
    # baselines + faint "ceiling" envelope (max-band height = maxpct%)
    geom_hline(yintercept = seq_along(ylabs), colour = "grey80", linewidth = 0.3) +
    geom_hline(yintercept = seq_along(ylabs) + ROW_HALF, colour = "grey90",
               linewidth = 0.2, linetype = "dotted") +
    geom_hline(yintercept = seq_along(ylabs) - ROW_HALF, colour = "grey90",
               linewidth = 0.2, linetype = "dotted") +
    geom_ribbon(aes(ymin = ymin, ymax = ymax,
                    group = interaction(term, dir, strain, phage), fill = fill_key),
                alpha = 0.92) +
    # mark the actual sampled timepoints on each band edge
    geom_point(aes(y = edge, colour = fill_key,
                   group = interaction(term, dir, strain, phage)),
               size = 0.7, show.legend = FALSE) +
    # ±% scale key (first facet only)
    geom_segment(data = key_df,
                 aes(x = min(TIMES), xend = min(TIMES),
                     y = y0, yend = y0 + ROW_HALF * (pct / maxpct)),
                 inherit.aes = FALSE, colour = "grey45", linewidth = 0.3) +
    geom_text(data = key_df,
              aes(x = min(TIMES), y = y0 + ROW_HALF * (pct / maxpct),
                  label = paste0(pct, "%")),
              inherit.aes = FALSE, hjust = 1.15, size = 2.1, colour = "grey45") +
    geom_text(data = sign_df, aes(x = max(TIMES), y = y, label = lab),
              inherit.aes = FALSE, hjust = 1, size = 2.2, colour = "grey45",
              fontface = "italic") +
    facet_nested(. ~ strain + phage,
                 strip = strip_nested(
                   background_x = elem_list_rect(
                     fill = c(STRAIN_COLOURS[STRAINS], rep("grey95", length(STRAINS) * length(PHAGES)))))) +
    scale_fill_manual(values = fill_values, guide = fill_guide, name = NULL) +
    scale_colour_manual(values = fill_values, guide = "none") +
    scale_y_continuous(breaks = seq_along(ylabs), labels = lab_map[ylabs],
                       expand = expansion(add = 0.6)) +
    scale_x_continuous(breaks = TIMES) +
    labs(title = paste0(title, "    (tallest band = ", round(maxpct), "% of a term's genes)"),
         x = "Time post-infection (min)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      strip.text = element_text(face = "bold", size = 9, margin = margin(2, 2, 2, 2)),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 7),
      panel.spacing.x = unit(0.25, "lines"),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom")
}

# ── 5. Assemble & save ────────────────────────────────────────────────────────
cat("── Rendering figure ──────────────────────────────────────────────────────\n")
p_cog  <- make_ridge("COG",  "A  COG functional categories")
p_kegg <- make_ridge("KEGG", "B  KEGG pathways")
panels <- Filter(Negate(is.null), list(p_cog, p_kegg))
if (length(panels) == 0) stop("No panels produced.")

combined <- wrap_plots(panels, ncol = 1, guides = "collect") +
  plot_annotation(
    title    = "Temporal, directional functional response to phage infection (vs control)",
    subtitle = paste0("Band = % of a term's tested genes that are DE (padj<", SIG_PADJ,
                      ", |log2FC|>=", SIG_LFC, ") at each timepoint.  ",
                      "Above baseline = induced, below = repressed.  ",
                      if (COLOUR_BY == "phage")
                        "Colour = strain hue; within each phage darker = induced, lighter = repressed.  "
                      else "",
                      "Terms = union-ORA enriched (FDR<", PADJ_FILTER, "), >=", MIN_PATHWAY_SIZE, " genes."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(size = 9, colour = "grey40"))) &
  theme(legend.position = "bottom")

n_rows <- sum(sapply(panels, function(p) length(levels(p$data$term))))
fig_h  <- max(8, n_rows * 0.55 + 3)

ggsave(FIG_OUT, combined, width = 13, height = fig_h, dpi = 300, bg = "white")
cat("\nSaved figure:", FIG_OUT, "\n")
cat("Tables in:", OUT_DIR, "\n")
