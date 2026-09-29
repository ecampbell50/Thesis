# Ribosome, galactose metabolism and PTS: % of pathway genes up/down vs control over time (Fig 4.4).

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler", update = FALSE)
for (pkg in c("ggh4x", "patchwork")) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)

suppressPackageStartupMessages({
  library(clusterProfiler); library(dplyr); library(readr); library(tidyr)
  library(stringr); library(ggplot2); library(patchwork); library(ggh4x); library(scales)
})

# ── 1. Paths & settings ───────────────────────────────────────────────────────
# Self-contained: all paths resolve relative to this script's location inside the
# paper_package (…/Figure3_KEGG_enrichment/{script,data,figure}). No external paths.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT      <- normalizePath(file.path(.script_dir, ".."))   # Figure3_KEGG_enrichment/
DATA      <- file.path(ROOT, "data")
DESEQ_DIR <- file.path(DATA, "DESeq2_vs_control")  # vs-control CSVs
ANNO_DIR  <- DATA                                  # functional_annotation_*_with_names.csv live here
OUT_DIR   <- file.path(ROOT, "figure")
FIG_OUT   <- file.path(OUT_DIR, "temporal_directional_enrichment_mainKEGG.png")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

STRAINS     <- c("C67", "D32", "D68")
PHAGES      <- c("BON", "CLY")
PHAGE_LABEL <- c(BON = "Bonnie", CLY = "Clyde")
TIMES       <- c(2, 10, 20, 30, 50)

SIG_PADJ <- 0.05   # "perturbed" = padj < this ...
SIG_LFC  <- 1      # ... AND |log2FC| >= this (same as your DEG definition)

MIN_PATHWAY_SIZE <- 1     # skip terms with < this many tested genes (noisy %)

# The three KEGG pathways shown in the main text (top → bottom in the figure).
# Names must match the eggNOG KEGG_Pathway_names strings exactly.
SEL_KEGG <- c("Ribosome", "Galactose metabolism", "Phosphotransferase system (PTS)")

STRAIN_COLOURS <- c(C67 = "#aed1ea", D32 = "#93f0a3", D68 = "#f98ab3")

# Per strain: darker shade = induced (up vs control), lighter = repressed (down).
# (Bonnie & Clyde share a strain's hue; phage is separated by the facet strips.)
PHAGE_COLS <- list(
  C67 = c(up = "#08306b", down = "#c6dbef"),   # navy / light blue
  D32 = c(up = "#586B24", down = "#C0D684"),   # dark / light green
  D68 = c(up = "#3D0B37", down = "#DDA6C6"))   # dark / light pink

# ── 2. Helpers ────────────────────────────────────────────────────────────────
clean_term_string <- function(x) { x <- str_trim(x); x[x == ""] <- NA_character_; x }

# TERM2GENE (KEGG pathway name → gene) from the eggNOG annotation file.
build_terms <- function(anno) {
  anno %>% dplyr::select(GeneID, KEGG_Pathway_names) %>%
    rename(val = KEGG_Pathway_names) %>%
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

# ── 3. Compute per-pathway, per-timepoint directional % (selected pathways only) ─
cat("── Computing temporal directional enrichment for main-text KEGG pathways ──\n")
pert_rows <- list()

for (strain in STRAINS) {
  af <- file.path(ANNO_DIR, sprintf("functional_annotation_%s_with_names.csv", strain))
  if (!file.exists(af)) stop("Missing annotation: ", af)
  anno <- read_csv(af, show_col_types = FALSE) %>%
    filter(!str_starts(GeneID, "PEN"), !str_starts(GeneID, "ZUD"))
  t2g <- build_terms(anno)
  membership <- split(t2g$gene, t2g$term)

  for (phage in PHAGES) {
    res <- read_vsctrl(strain, phage)
    if (nrow(res) == 0) next
    tested <- unique(res$gene)

    for (term in SEL_KEGG) {
      pg <- intersect(membership[[term]], tested)
      n_tested <- length(pg)
      if (n_tested < MIN_PATHWAY_SIZE) {
        message(sprintf("  [%s/%s] '%s' has %d tested genes — skipped", strain, phage, term, n_tested))
        next
      }
      sub <- res %>% filter(gene %in% pg)
      for (tp in TIMES) {
        st <- sub %>% filter(time == tp, !is.na(padj), !is.na(log2FoldChange))
        up <- sum(st$padj < SIG_PADJ & st$log2FoldChange >=  SIG_LFC)
        dn <- sum(st$padj < SIG_PADJ & st$log2FoldChange <= -SIG_LFC)
        pert_rows[[paste(strain, phage, term, tp)]] <- tibble(
          strain = strain, phage = phage, term = term, time = tp,
          n_tested = n_tested, up_pct = 100 * up / n_tested, down_pct = 100 * dn / n_tested)
      }
    }
  }
}

pert <- bind_rows(pert_rows)
if (nrow(pert) == 0) stop("No data computed for the selected pathways.")
write_csv(pert, file.path(OUT_DIR, "temporal_KEGG_mainText_directional_table.csv"))

# ── 4. Build the figure: one panel per pathway, real 0–100% directional y-axis ─
cat("── Rendering main-text KEGG figure ───────────────────────────────────────\n")

# Pathway-strip labels carry the (mean) tested-gene count.
term_n   <- pert %>% group_by(term) %>% summarise(n = round(mean(n_tested)), .groups = "drop")
term_lab <- setNames(sprintf("%s\n(n≈%d genes)", term_n$term, term_n$n), term_n$term)

# Long ribbon table: induced positive (above 0), repressed negative (below 0).
rib <- bind_rows(
  pert %>% transmute(strain, phage, term, time, dir = "up",
                     ymin = 0,         ymax = up_pct,  edge = up_pct),
  pert %>% transmute(strain, phage, term, time, dir = "down",
                     ymin = -down_pct, ymax = 0,       edge = -down_pct)
) %>%
  mutate(
    strain   = factor(strain, levels = STRAINS),
    phage    = factor(PHAGE_LABEL[phage], levels = unname(PHAGE_LABEL[PHAGES])),
    term     = factor(term, levels = SEL_KEGG),
    fill_key = paste0(strain, ifelse(dir == "up", " induced", " repressed"))
  ) %>%
  arrange(time)

# Colour map + legend order (induced then repressed, grouped by strain).
fill_levels <- as.vector(rbind(paste0(STRAINS, " induced"), paste0(STRAINS, " repressed")))
fill_values <- setNames(
  unlist(lapply(STRAINS, function(s) c(PHAGE_COLS[[s]][["up"]], PHAGE_COLS[[s]][["down"]]))),
  fill_levels)
rib$fill_key <- factor(rib$fill_key, levels = fill_levels)

p <- ggplot(rib, aes(x = time)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax,
                  group = interaction(term, dir, strain, phage), fill = fill_key),
              alpha = 0.92) +
  geom_point(aes(y = edge, colour = fill_key,
                 group = interaction(term, dir, strain, phage)),
             size = 0.7, show.legend = FALSE) +
  facet_nested(term ~ strain + phage, switch = "y",
               labeller = labeller(term = term_lab),
               strip = strip_nested(
                 background_x = elem_list_rect(
                   fill = c(STRAIN_COLOURS[STRAINS],
                            rep("grey95", length(STRAINS) * length(PHAGES)))))) +
  scale_fill_manual(values = fill_values, name = NULL,
                    guide = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  scale_colour_manual(values = fill_values, guide = "none") +
  scale_y_continuous(
    limits = c(-100, 100), breaks = seq(-100, 100, 50),
    labels = function(b) paste0(abs(b), "%"),
    minor_breaks = seq(-100, 100, 25)) +
  scale_x_continuous(breaks = TIMES) +
  labs(
    title    = "Temporal, directional response of key KEGG pathways to phage infection (vs control)",
    subtitle = paste0("Band = % of a pathway's tested genes that are DE (padj<", SIG_PADJ,
                      ", |log2FC|>=", SIG_LFC, ") at each timepoint.\n",
                      "Above 0 = induced (up vs control), below 0 = repressed.  ",
                      "Colour = strain hue; darker = induced, lighter = repressed."),
    x = "Time post-infection (min)",
    y = "% of pathway genes  (▲ induced  /  ▼ repressed)") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.minor.y = element_line(colour = "grey94", linewidth = 0.2),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    strip.text   = element_text(face = "bold", size = 9, margin = margin(3, 3, 3, 3)),
    strip.text.y.left = element_text(face = "bold", size = 9, angle = 0),
    strip.placement = "outside",
    axis.text.y  = element_text(size = 7),
    axis.text.x  = element_text(size = 7),
    panel.spacing.x = unit(0.25, "lines"),
    panel.spacing.y = unit(0.4, "lines"),
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, colour = "grey40"),
    legend.position = "bottom")

fig_h <- length(SEL_KEGG) * 1.7 + 2.5
ggsave(FIG_OUT, p, width = 13, height = fig_h, dpi = 300, bg = "white")
cat("\nSaved figure:", FIG_OUT, "\n")
cat("Table in:", OUT_DIR, "\n")
