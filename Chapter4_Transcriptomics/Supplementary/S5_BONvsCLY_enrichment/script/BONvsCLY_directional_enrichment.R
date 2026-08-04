# ══════════════════════════════════════════════════════════════════════════════
# BONvsCLY_directional_enrichment.R
#
# WHAT THIS SCRIPT DOES
# ─────────────────────
# Functional enrichment (COG / KEGG / GO) of the *direct* Bonnie-vs-Clyde host
# DE results, split by DIRECTION of change — exactly the analysis Chris suggested:
# "separate the genes that are up in Bonnie from those up in Clyde and look at the
#  pathways each indicates."
#
# It is a single self-contained merge of two existing scripts:
#   - enrichment computation  ← *-analysis-strept.R  (clusterProfiler::enricher,
#                                TERM2GENE from functional_annotation_{STRAIN}_with_names.csv)
#   - combined figure         ← allTP_overall_enrichment_figure.R (stacked dotplots)
#
# PIPELINE (per strain)
#   1. Read the Bonnie-vs-Clyde DESeq2 CSVs for each timepoint (T2,10,20,30,50):
#        results/deseq2_BONvsCLY/DESeq2_{STRAIN}_BONvsCLY_T{tp}min.csv
#      (Clyde is the reference, so +log2FC = higher under Bonnie.)
#   2. Take the UNION across timepoints of significant genes (padj<0.05 & |LFC|>=1)
#      and split into two sets:
#        Bonnie-up : log2FC >= +1   (higher under Bonnie than Clyde)
#        Clyde-up  : log2FC <= -1   (higher under Clyde than Bonnie)
#   3. Run ONE Fisher/hypergeometric enrichment (enricher) per direction, against
#      the same genome background (all tested genes ∩ annotated genes).
#   4. Write per-condition enrichment CSVs, then render one figure:
#        x-axis  = direction (Bonnie-up / Clyde-up)
#        facets  = strain (C67 / D32 / D68)
#        panels  = COG / KEGG / GO (stacked)
#
# IMPORTANT INTERPRETATION CAVEAT
#   This is a phage-vs-phage contrast with NO uninfected baseline, so "Bonnie-up"
#   means *higher under Bonnie than under Clyde* — NOT up vs resting. Enriched
#   terms describe how the two infection programmes DIFFER, not absolute induction.
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Packages ───────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler", update = FALSE)
for (pkg in c("ggh4x", "patchwork", "forcats")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  library(scales)
  library(ggh4x)   # per-strain coloured facet strips
})

# ── 1. Paths & global thresholds ──────────────────────────────────────────────
# Self-contained: resolve relative to this script's location in paper_package
# (…/S5_BONvsCLY_enrichment/{script,data,figure}). No external paths.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT      <- normalizePath(file.path(.script_dir, ".."))
DATA      <- file.path(ROOT, "data")
DESEQ_DIR <- file.path(DATA, "deseq2_BONvsCLY")
ANNO_DIR  <- DATA                                  # functional_annotation_*_with_names.csv live here
OUT_DIR   <- file.path(ROOT, "figure")
FIG_OUT   <- file.path(OUT_DIR, "BONvsCLY_directional_enrichment_figure.png")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STRAINS <- c("C67", "D32", "D68")
TIMES   <- c(2, 10, 20, 30, 50)

# DEG selection (matches the strept / allTP pipeline)
PADJ_THR <- 0.05
LFC_THR  <- 1

# Directions of the Bonnie-vs-Clyde contrast (Clyde = reference, +LFC = Bonnie)
#   internal key  → display label (used on the x-axis and in filenames)
DIRECTIONS  <- c("BONup", "CLYup")
DIR_LABELS  <- c(BONup = "Bonnie-up", CLYup = "Clyde-up")
DIR_COLOURS <- c(BONup = "#1b7837", CLYup = "#762a83")  # axis-tick label colours

# Figure significance filter — FDR (Benjamini-Hochberg padj), not nominal p.
PADJ_FILTER <- 0.1
TOPN_COG  <- 15
TOPN_KEGG <- 15
TOPN_GO   <- 15

# Strain facet-strip background colours (consistent with your other figures)
STRAIN_COLOURS <- c(C67 = "#aed1ea", D32 = "#93f0a3", D68 = "#f98ab3")

# ── 2. Helpers ────────────────────────────────────────────────────────────────
clean_term_string <- function(x) { x <- str_trim(x); x[x == ""] <- NA_character_; x }

# Build the union DEG sets (per direction) + background for one strain.
load_strain_genes <- function(strain) {
  tested <- character(0); bon_up <- character(0); cly_up <- character(0)
  for (tp in TIMES) {
    f <- file.path(DESEQ_DIR, sprintf("DESeq2_%s_BONvsCLY_T%dmin.csv", strain, tp))
    if (!file.exists(f)) { message("  [MISSING] ", f); next }
    df <- read_csv(f, show_col_types = FALSE)
    if ("...1" %in% colnames(df)) colnames(df)[colnames(df) == "...1"] <- "gene"
    # Host genes only (defensive: drop any phage IDs if present)
    df <- df %>% filter(!str_starts(gene, "PEN"), !str_starts(gene, "ZUD"))
    tested <- union(tested, unique(df$gene))
    sig <- df %>% filter(!is.na(padj), padj < PADJ_THR, !is.na(log2FoldChange))
    bon_up <- union(bon_up, sig %>% filter(log2FoldChange >=  LFC_THR) %>% pull(gene))
    cly_up <- union(cly_up, sig %>% filter(log2FoldChange <= -LFC_THR) %>% pull(gene))
  }
  list(tested = tested, BONup = bon_up, CLYup = cly_up)
}

# Build TERM2GENE (and TERM2NAME for GO) for one annotation type.
build_terms <- function(anno, type) {
  if (type == "COG") {
    t2g <- anno %>% dplyr::select(GeneID, COG_names) %>%
      filter(!is.na(COG_names), COG_names != "") %>%
      separate_rows(COG_names, sep = "; ") %>%
      mutate(COG_names = clean_term_string(COG_names)) %>%
      filter(!is.na(COG_names)) %>%
      transmute(term = COG_names, gene = GeneID) %>% distinct()
    return(list(t2g = t2g, t2n = NULL))
  }
  if (type == "KEGG") {
    t2g <- anno %>% dplyr::select(GeneID, KEGG_Pathway_names) %>%
      filter(!is.na(KEGG_Pathway_names), KEGG_Pathway_names != "") %>%
      separate_rows(KEGG_Pathway_names, sep = "; ") %>%
      mutate(KEGG_Pathway_names = clean_term_string(KEGG_Pathway_names)) %>%
      filter(!is.na(KEGG_Pathway_names)) %>%
      transmute(term = KEGG_Pathway_names, gene = GeneID) %>% distinct()
    return(list(t2g = t2g, t2n = NULL))
  }
  if (type == "GO") {  # all ontologies combined (BP+MF+CC), like GO_ALL_results
    t2g <- anno %>% dplyr::select(GeneID, GOs) %>%
      filter(!is.na(GOs), GOs != "") %>%
      separate_rows(GOs, sep = ";") %>%
      mutate(GOs = clean_term_string(GOs)) %>%
      filter(!is.na(GOs)) %>%
      transmute(term = GOs, gene = GeneID) %>% distinct()
    t2n <- anno %>% dplyr::select(GOs, GO_names) %>%
      filter(!is.na(GOs), GOs != "", !is.na(GO_names), GO_names != "") %>%
      separate_rows(GOs, sep = ";") %>%
      separate_rows(GO_names, sep = "; ") %>%
      mutate(GOs = clean_term_string(GOs), GO_names = clean_term_string(GO_names)) %>%
      filter(!is.na(GOs), !is.na(GO_names)) %>%
      transmute(term = GOs, name = GO_names) %>% distinct(term, .keep_all = TRUE)
    return(list(t2g = t2g, t2n = t2n))
  }
  stop("Unknown type: ", type)
}

# ── 3. Load annotations & gene sets once per strain ───────────────────────────
cat("── Loading annotations and Bonnie-vs-Clyde DEG sets ──────────────────────\n")
annos <- list(); strain_genes <- list()
for (s in STRAINS) {
  af <- file.path(ANNO_DIR, paste0("functional_annotation_", s, "_with_names.csv"))
  if (!file.exists(af)) stop("Missing annotation file: ", af)
  a <- read_csv(af, show_col_types = FALSE) %>%
    filter(!str_starts(GeneID, "PEN"), !str_starts(GeneID, "ZUD"))
  annos[[s]] <- a
  g <- load_strain_genes(s)
  strain_genes[[s]] <- g
  cat(sprintf("  %s — tested: %d | Bonnie-up: %d | Clyde-up: %d\n",
              s, length(g$tested), length(g$BONup), length(g$CLYup)))
}

# ── 4. Run enrichment for every (strain × direction × type) ───────────────────
# Returns one tidy long-format dataframe per term type, ready for the figure.
collect_type <- function(type) {
  rows <- list()
  for (s in STRAINS) {
    terms <- build_terms(annos[[s]], type)
    if (nrow(terms$t2g) == 0) next
    g  <- strain_genes[[s]]
    bg <- intersect(g$tested, unique(annos[[s]]$GeneID))
    for (d in DIRECTIONS) {
      sig <- intersect(g[[d]], bg)
      if (length(sig) == 0) { message("  ", type, " ", s, " ", d, ": no genes"); next }
      e <- enricher(gene = sig, universe = bg,
                    TERM2GENE = terms$t2g, TERM2NAME = terms$t2n,
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
      if (is.null(e)) next
      edf <- as.data.frame(e)
      if (nrow(edf) == 0) next
      if (!"Count" %in% colnames(edf)) edf$Count <- lengths(strsplit(as.character(edf$geneID), "/"))

      # Persist the full per-condition table (mirrors the strept pipeline outputs)
      write_csv(edf, file.path(OUT_DIR,
                 sprintf("%s_%s_%s_enrichment_results.csv", s, d, type)))

      rows[[paste(s, d)]] <- edf %>%
        mutate(
          term_raw  = ifelse(is.na(Description) | Description == "", ID, str_trim(Description)),
          term      = str_wrap(term_raw, width = 40),
          Strain    = s,
          Direction = d,
          Condition = paste0(s, "\n", DIR_LABELS[[d]])
        ) %>%
        dplyr::select(term, term_raw, pvalue, p.adjust, Count, Strain, Direction, Condition)
    }
  }
  if (length(rows) == 0) return(NULL)
  bind_rows(rows) %>%
    group_by(term, Condition) %>%
    slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()
}

# ── 5. Build one figure panel for a term type ─────────────────────────────────
make_panel <- function(type, topN, title_label) {
  dat <- collect_type(type)
  if (is.null(dat) || nrow(dat) == 0) {
    message(type, ": no enrichment results — skipping panel."); return(NULL)
  }

  # Show a term if FDR-significant in >=1 condition; rank by best padj; take topN.
  sig_terms <- dat %>% filter(!is.na(p.adjust), p.adjust < PADJ_FILTER) %>% pull(term) %>% unique()
  if (length(sig_terms) == 0) {
    message(type, ": no terms pass padj < ", PADJ_FILTER, " — skipping panel."); return(NULL)
  }
  term_order <- dat %>% filter(term %in% sig_terms) %>%
    group_by(term) %>% summarise(best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_padj) %>% slice_head(n = topN) %>% pull(term)

  # Full grid so the dotplot reads as a matrix (missing = blank).
  full_grid <- expand.grid(term = term_order, Strain = STRAINS,
                           Direction = DIRECTIONS, stringsAsFactors = FALSE)
  dat_plot <- full_grid %>%
    left_join(dat %>% filter(term %in% term_order),
              by = c("term", "Strain", "Direction")) %>%
    mutate(
      term      = factor(term, levels = rev(term_order)),
      Strain    = factor(Strain, levels = STRAINS),
      Direction = factor(Direction, levels = DIRECTIONS),
      neg_log_p = pmin(-log10(pvalue + 1e-300), 15)
    )

  d_draw <- dat_plot %>% filter(!is.na(p.adjust), p.adjust < PADJ_FILTER)
  if (nrow(d_draw) == 0) { message(type, ": nothing to draw."); return(NULL) }

  p_breaks  <- c(0.05, 0.01, 0.001, 1e-5, 1e-10)
  break_pos <- -log10(p_breaks)
  break_lbl <- c("0.05", "0.01", "0.001", expression(10^-5), expression(10^-10))

  x_text_cols <- DIR_COLOURS[DIRECTIONS]

  present_strains <- levels(droplevels(d_draw$Strain))
  strip_rects <- lapply(present_strains, function(s)
    element_rect(fill = STRAIN_COLOURS[[s]], colour = NA))

  ggplot(d_draw, aes(x = Direction, y = term)) +
    geom_point(aes(size = Count, fill = neg_log_p),
               shape = 21, colour = "black", stroke = 0.5, alpha = 0.95) +
    facet_grid2(. ~ Strain, strip = strip_themed(background_x = strip_rects)) +
    scale_fill_gradientn(
      colours  = c("#4575B4", "#74ADD1", "#FEE090", "#F46D43", "#A50026"),
      limits   = c(-log10(PADJ_FILTER), 15),
      breaks   = break_pos, labels = break_lbl,
      name     = "Fisher's\nexact p", oob = scales::squish, na.value = NA) +
    scale_size_continuous(name = "DEGs in\nterm", range = c(2, 12),
                          breaks = c(2, 5, 10, 25, 50, 100)) +
    scale_x_discrete(drop = FALSE, labels = DIR_LABELS[DIRECTIONS]) +
    labs(title = title_label, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 9),
      axis.text.x      = element_text(size = 9.5, face = "bold", colour = x_text_cols),
      strip.text.x     = element_text(size = 10, face = "bold", margin = margin(t = 3, b = 3)),
      panel.spacing.x  = unit(0.6, "lines"),
      plot.title       = element_text(size = 11, face = "bold", margin = margin(b = 4)),
      legend.position  = "right",
      legend.key.size  = unit(0.5, "cm"),
      legend.text      = element_text(size = 8),
      legend.title     = element_text(size = 8.5),
      plot.margin      = margin(4, 4, 4, 4)) +
    guides(
      fill = guide_colourbar(order = 1, barwidth = 0.7, barheight = 5,
                             title.position = "top", reverse = TRUE),
      size = guide_legend(order = 2,
                          override.aes = list(shape = 21, fill = "grey70", colour = "black")))
}

# ── 6. Assemble the three panels ──────────────────────────────────────────────
cat("── Running enrichment & building figure ──────────────────────────────────\n")
p_cog  <- make_panel("COG",  TOPN_COG,  "A  COG functional categories")
p_kegg <- make_panel("KEGG", TOPN_KEGG, "B  KEGG pathways")
p_go   <- make_panel("GO",   TOPN_GO,   "C  Gene Ontology (BP/MF/CC combined)")

panels <- Filter(Negate(is.null), list(p_cog, p_kegg, p_go))
if (length(panels) == 0) stop("No panels produced — check inputs.")

combined <- wrap_plots(panels, ncol = 1) +
  plot_annotation(
    title    = "Bonnie vs Clyde: directional functional enrichment (host genes)",
    subtitle = paste0("Union of DEGs across T2-T50 (|log2FC| >= ", LFC_THR,
                      ", padj < ", PADJ_THR, "), split by direction.  ",
                      "Bonnie-up = higher under Bonnie than Clyde; Clyde-up = vice versa.  ",
                      "Only terms with FDR < ", PADJ_FILTER, " shown."),
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey40")))

total_rows <- sum(sapply(panels, function(p) length(levels(p$data$term))))
fig_h <- max(12, total_rows * 0.40 + 5)

ggsave(FIG_OUT, combined, width = 11, height = fig_h, dpi = 300, bg = "white")
cat("\nSaved figure:", FIG_OUT, "\n")
cat("Per-condition enrichment CSVs in:", OUT_DIR, "\n")
