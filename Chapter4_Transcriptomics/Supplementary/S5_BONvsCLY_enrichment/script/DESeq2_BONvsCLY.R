# DESeq2_BONvsCLY.R ── Host gene DE: Bonnie vs Clyde at each timepoint
#
# Instead of comparing each phage to the uninfected control, this script
# contrasts the two phage infections directly: at every timepoint, the host
# transcriptome under Bonnie infection is compared to the host transcriptome
# under Clyde infection (same strain, same timepoint, 3 reps each).
#
# Clyde is the reference level, so:
#   positive log2FC = higher under Bonnie
#   negative log2FC = higher under Clyde
#
# DESeq2 (per strain × timepoint) and the figure are produced in this one script.
#
# ── RUN FROM ──────────────────────────────────────────────────────────────────
#   setwd("...2_Analysis/")
#   source("scripts/R/DESeq2_BONvsCLY.R")
#
# ── DATA IN ───────────────────────────────────────────────────────────────────
#   data/raw_counts/{STRAIN}_BON_full_raw_counts.tsv
#   data/raw_counts/{STRAIN}_CLY_full_raw_counts.tsv
#   Columns: Geneid, 0_R1..50_R3 (18 counts), Entity, Symbol
#
# ── OUTPUT ────────────────────────────────────────────────────────────────────
#   results/deseq2_BONvsCLY/DESeq2_{STRAIN}_BONvsCLY_T{tp}min.csv
#   results/figures/fig_host_BONvsCLY.pdf
# ══════════════════════════════════════════════════════════════════════════════

# ── 0. Packages ───────────────────────────────────────────────────────────────
if (!("DESeq2" %in% installed.packages())) BiocManager::install("DESeq2", update = FALSE)
if (!("ashr"   %in% installed.packages())) BiocManager::install("ashr",   update = FALSE)
if (!("tidyverse" %in% installed.packages())) install.packages("tidyverse")
if (!("ggrepel"   %in% installed.packages())) install.packages("ggrepel")
library(DESeq2)
library(ashr)
library(tidyverse)
library(ggrepel)

# ── 1. Configuration ──────────────────────────────────────────────────────────
STRAINS    <- c("C67", "D32", "D68")
TIMEPOINTS <- c(2, 10, 20, 30, 50)   # T=0 excluded (pre-infection baseline)
REF_PHAGE  <- "CLY"                  # reference level → +LFC = higher under Bonnie

# Self-contained: resolve relative to this script's location in paper_package
# (…/S5_BONvsCLY_enrichment/{script,data,figure}). No external paths.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT       <- normalizePath(file.path(.script_dir, ".."))
COUNTS_DIR <- file.path(ROOT, "data", "raw_counts")
CSV_DIR    <- file.path(ROOT, "data", "deseq2_BONvsCLY")
FIG_DIR    <- file.path(ROOT, "figure")
dir.create(CSV_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Samples to exclude (failed QC), as "{PHAGE}_{tp}_{rep}", e.g. "BON_0_R1".
# Set to character(0) if none.
EXCLUDE_SAMPLES <- character(0)

PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1      # below this: grey; at or above: coloured

# Optionally label individual genes on the figure (gene-number text via ggrepel).
# Toggle either independently; both off = clean dot plot.
LABEL_DEFENCE  <- TRUE
LABEL_PROPHAGE <- FALSE

# ── 2. Annotation lists ───────────────────────────────────────────────────────
defence_genes <- list(
  C67 = c(
    "C67_prot_00121","C67_prot_00122","C67_prot_00133","C67_prot_00134",
    "C67_prot_00284","C67_prot_00285","C67_prot_00286","C67_prot_00287",
    "C67_prot_00297","C67_prot_00298","C67_prot_01021","C67_prot_01022",
    "C67_prot_01023","C67_prot_01024","C67_prot_01211","C67_prot_01330",
    "C67_prot_01331","C67_prot_01418","C67_prot_01435","C67_prot_01476",
    "C67_prot_00060","C67_prot_00603","C67_prot_01356","C67_prot_00253",
    "C67_prot_00763","C67_prot_01000","C67_prot_01545","C67_prot_01598",
    "C67_prot_00210","C67_prot_00597","C67_prot_00619","C67_prot_01170",
    "C67_prot_01206","C67_prot_01226","C67_prot_01617","C67_prot_01819"
  ),
  D32 = c(
    "D32_prot_00742","D32_prot_00009","D32_prot_00085","D32_prot_00247",
    "D32_prot_00257","D32_prot_00379","D32_prot_00401","D32_prot_00920",
    "D32_prot_00940","D32_prot_01852","D32_prot_01911","D32_prot_01263",
    "D32_prot_01431","D32_prot_00042","D32_prot_01386","D32_prot_01328",
    "D32_prot_00172","D32_prot_00173","D32_prot_00292","D32_prot_00395",
    "D32_prot_01175","D32_prot_00925","D32_prot_00008","D32_prot_00010",
    "D32_prot_00011","D32_prot_00882","D32_prot_00883","D32_prot_00884",
    "D32_prot_01452","D32_prot_01453","D32_prot_01454","D32_prot_01455",
    "D32_prot_01970","D32_prot_01971","D32_prot_00160","D32_prot_00161",
    "D32_prot_01200","D32_prot_01201","D32_prot_01867","D32_prot_01369"
  ),
  D68 = c(
    "D68_prot_01256","D68_prot_01327","D68_prot_00729","D68_prot_00232",
    "D68_prot_00254","D68_prot_00291","D68_prot_00367","D68_prot_00528",
    "D68_prot_00548","D68_prot_00919","D68_prot_00929","D68_prot_01833",
    "D68_prot_01915","D68_prot_01201","D68_prot_01696","D68_prot_00324",
    "D68_prot_01102","D68_prot_01160","D68_prot_00237","D68_prot_00238",
    "D68_prot_01046","D68_prot_00533","D68_prot_00289","D68_prot_00292",
    "D68_prot_00293","D68_prot_00490","D68_prot_00491","D68_prot_00492",
    "D68_prot_01222","D68_prot_01223","D68_prot_01224","D68_prot_01225",
    "D68_prot_01973","D68_prot_01974","D68_prot_00442","D68_prot_00443",
    "D68_prot_01020","D68_prot_01021","D68_prot_01871","D68_prot_01119"
  )
)

prophage_genes <- list(
  C67 = c("C67_prot_02009","C67_prot_00947","C67_prot_00948","C67_prot_00949",
           "C67_prot_00950","C67_prot_00951","C67_prot_00952","C67_prot_00953",
           "C67_prot_00954","C67_prot_00955","C67_prot_00956","C67_prot_00957",
           "C67_prot_00958","C67_prot_00959","C67_prot_00960","C67_prot_00961",
           "C67_prot_00962","C67_prot_00963","C67_prot_00964","C67_prot_00965"),
  D32 = c("D32_prot_01483","D32_prot_01484","D32_prot_01485","D32_prot_01486",
           "D32_prot_01487","D32_prot_01488","D32_prot_01489","D32_prot_01490",
           "D32_prot_01491","D32_prot_01492","D32_prot_01493","D32_prot_01494",
           "D32_prot_01495","D32_prot_01496","D32_prot_01497","D32_prot_01498",
           "D32_prot_01499","D32_prot_01500","D32_prot_01501"),
  D68 = c("D68_prot_01332","D68_prot_01333","D68_prot_01334","D68_prot_01335",
           "D68_prot_01336","D68_prot_01337","D68_prot_01338","D68_prot_01339",
           "D68_prot_01340","D68_prot_01341","D68_prot_01342","D68_prot_01343",
           "D68_prot_01344","D68_prot_01345","D68_prot_01346","D68_prot_01347",
           "D68_prot_01348","D68_prot_01349","D68_prot_01350")
)

# ── 3. Helper: load host-gene counts for one strain+phage ─────────────────────
# Keeps only host rows (Geneid begins "{STRAIN}_") and drops the Entity/Symbol
# columns, leaving a gene × sample integer matrix with columns "{tp}_{rep}".
load_host_counts <- function(strain, phage) {
  f <- file.path(COUNTS_DIR, sprintf("%s_%s_full_raw_counts.tsv", strain, phage))
  if (!file.exists(f)) stop("Missing counts file: ", f)
  df  <- read.delim(f, row.names = 1, check.names = FALSE)
  df  <- df[, !colnames(df) %in% c("Entity", "Symbol"), drop = FALSE]
  df  <- df[grepl(paste0("^", strain, "_"), rownames(df)), , drop = FALSE]
  as.matrix(df)
}

# ── 4. DESeq2: Bonnie vs Clyde per strain × timepoint ─────────────────────────
cat("── Running DESeq2 (Bonnie vs Clyde) ──────────────────────────────────────\n")

all_results <- list()

for (strain in STRAINS) {
  bon <- load_host_counts(strain, "BON")
  cly <- load_host_counts(strain, "CLY")

  # Host gene sets are identical between the two files for a strain, but align
  # defensively on the shared gene set just in case.
  shared_genes <- intersect(rownames(bon), rownames(cly))
  bon <- bon[shared_genes, , drop = FALSE]
  cly <- cly[shared_genes, , drop = FALSE]

  # Tag columns by phage so the two matrices don't collide once combined.
  colnames(bon) <- paste0("BON_", colnames(bon))
  colnames(cly) <- paste0("CLY_", colnames(cly))
  mat <- cbind(bon, cly)

  for (tp in TIMEPOINTS) {
    # Pick the 6 columns (3 BON + 3 CLY) for this timepoint.
    cols <- colnames(mat)[grepl(paste0("_", tp, "_R[0-9]+$"), colnames(mat))]
    cols <- setdiff(cols, EXCLUDE_SAMPLES)

    sub_mat <- mat[, cols, drop = FALSE]
    phage   <- factor(sub("_.*", "", cols), levels = c(REF_PHAGE,
                       setdiff(c("BON", "CLY"), REF_PHAGE)))

    coldata <- data.frame(row.names = cols, phage = phage)

    dds <- DESeqDataSetFromMatrix(countData = sub_mat,
                                  colData   = coldata,
                                  design    = ~ phage)
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    dds <- DESeq(dds)

    coef_name <- grep("phage", resultsNames(dds), value = TRUE)
    res <- lfcShrink(dds, coef = coef_name, type = "ashr")

    out_csv <- file.path(CSV_DIR,
                         sprintf("DESeq2_%s_BONvsCLY_T%dmin.csv", strain, tp))
    write.csv(as.data.frame(res), out_csv)

    df <- as.data.frame(res)
    df$gene      <- rownames(df)
    df$strain    <- strain
    df$timepoint <- tp
    all_results[[paste(strain, tp)]] <- df

    cat(sprintf("  %s  T=%2dmin  (%s)  genes=%d  sig(padj<%.2g)=%d\n",
                strain, tp, coef_name, nrow(df),
                PADJ_THRESHOLD, sum(df$padj < PADJ_THRESHOLD, na.rm = TRUE)))
  }
}

full_df <- bind_rows(all_results)

# ── 5. Filter and categorise ──────────────────────────────────────────────────
plot_df <- full_df %>%
  filter(!is.na(padj), padj < PADJ_THRESHOLD) %>%
  mutate(
    is_defence  = mapply(function(g, s) g %in% defence_genes[[s]],  gene, strain),
    is_prophage = mapply(function(g, s) g %in% prophage_genes[[s]], gene, strain),
    category = case_when(
      is_defence  ~ "defence",
      is_prophage ~ "prophage",
      abs(log2FoldChange) >= LFC_THRESHOLD ~ "DE",
      TRUE ~ "small"
    ),
    category   = factor(category, levels = c("small", "DE", "prophage", "defence")),
    strain     = factor(strain, levels = STRAINS),
    timepoint  = factor(timepoint, levels = TIMEPOINTS),
    gene_label = sub(".*_", "", gene)   # compact id (e.g. "00121") for labelling
  )

cat(sprintf("\n  Genes passing padj < %g: %d rows\n", PADJ_THRESHOLD, nrow(plot_df)))
plot_df %>%
  dplyr::count(strain, category) %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0) %>%
  print()

tp_counts <- plot_df %>%
  group_by(strain, timepoint) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(strain = factor(strain, levels = STRAINS),
         label  = paste0("n = \n ", n))

# ── 6. Plot ───────────────────────────────────────────────────────────────────
pt_colours <- c(small = "grey72", DE = "#ccff00",
                prophage = "#f04a95", defence = "#4066ff")
pt_sizes <- c(small = 1, DE = 1.0, prophage = 1, defence = 1)
pt_alpha <- c(small = 0.20, DE = 0.85, prophage = 1, defence = 1)

# Visible y-window. Bonnie-vs-Clyde LFCs run roughly -5 to +17 (the long
# positive tail is rare multi-mapping artefacts), so the window is shallow
# below zero and tall above it — the mirror image of the phage-vs-control fig.
Y_BOT  <- -5     # data bottoms out near -4.8; nothing useful below this
Y_TOP  <- 20
NLAB_Y <- 19.5   # y-position of the "n =" counts (up in the top whitespace)

fig <- ggplot(plot_df,
              aes(x = timepoint, y = log2FoldChange,
                  colour = category, size = category, alpha = category)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
  geom_jitter(data = subset(plot_df, category == "small"),    width = 0.28, height = 0, shape = 19) +
  geom_jitter(data = subset(plot_df, category == "DE"),       width = 0.35, height = 0, shape = 20) +
  geom_jitter(data = subset(plot_df, category == "prophage"), width = 0.2,  height = 0, shape = 17) +
  geom_jitter(data = subset(plot_df, category == "defence"),  width = 0.2,  height = 0, shape = 15) +
  geom_text(data = tp_counts,
            aes(x = timepoint, y = NLAB_Y, label = label),
            inherit.aes = FALSE,
            vjust = 1, size = 3.4, colour = "grey30", fontface = "bold") +
  scale_colour_manual(
    values = pt_colours,
    labels = c(small    = paste0("Sig., |LFC| < ", LFC_THRESHOLD),
               DE       = paste0("DE,  |LFC| ≥ ", LFC_THRESHOLD),
               prophage = "Prophage",
               defence  = "Defence"),
    name = NULL
  ) +
  scale_size_manual(values  = pt_sizes, guide = "none") +
  scale_alpha_manual(values = pt_alpha, guide = "none") +
  scale_x_discrete(labels = paste0(TIMEPOINTS, " min")) +
  scale_y_continuous(
    trans  = scales::pseudo_log_trans(sigma = 1, base = 2),
    breaks = c(-5, -2, -1, 0, 1, 2, 5, 10, 20),
    labels = c(-5, -2, -1, 0, 1, 2, 5, 10, 20),
    expand = expansion(mult = c(0.01, 0.05))   # tight at bottom, small pad on top
  ) +
  coord_cartesian(ylim = c(Y_BOT, Y_TOP), clip = "off") +
  geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", colour = "black", linewidth = 0.45) +
  facet_grid(. ~ strain) +
  labs(
    title    = "Host gene differential expression: Bonnie vs Clyde",
    subtitle = paste0("padj < ", PADJ_THRESHOLD,
                      "  ·  log₂FC (Bonnie vs Clyde)  ·  pseudo-log₂ y-axis  ",
                      "·  +ve = higher under Bonnie"),
    x = "Time post-infection",
    y = expression(log[2]~"Fold Change  (pseudo-log scale)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.background   = element_rect(fill = "white"),
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5,
                                      margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 9, colour = "grey40", hjust = 0.5,
                                      margin = margin(b = 8)),
    plot.margin        = margin(t = 12, r = 12, b = 8, l = 8),
    strip.text         = element_text(face = "bold", size = 10),
    strip.background   = element_rect(fill = "grey93"),
    legend.position    = "bottom",
    legend.key.size    = unit(0.9, "lines"),
    legend.text        = element_text(size = 9),
    legend.margin      = margin(t = 4),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.spacing.x    = unit(0.9, "lines"),
    axis.text.x        = element_text(size = 8.5),
    axis.text.y        = element_text(size = 8),
    axis.title         = element_text(size = 10)
  )

# ── 7. Optional gene labels (defence / prophage) ──────────────────────────────
# Labels point to un-jittered positions, so they sit at the gene's true x/y
# rather than its jittered dot — close enough to read off which point is which.
if (LABEL_DEFENCE) {
  fig <- fig + geom_text_repel(
    data = subset(plot_df, category == "defence"),
    aes(x = timepoint, y = log2FoldChange, label = gene_label),
    inherit.aes = FALSE, size = 2.2, colour = "#1f3fb0",
    max.overlaps = Inf, box.padding = 0.25, segment.size = 0.2, segment.alpha = 0.5)
}
if (LABEL_PROPHAGE) {
  fig <- fig + geom_text_repel(
    data = subset(plot_df, category == "prophage"),
    aes(x = timepoint, y = log2FoldChange, label = gene_label),
    inherit.aes = FALSE, size = 2.2, colour = "#b01060",
    max.overlaps = Inf, box.padding = 0.25, segment.size = 0.2, segment.alpha = 0.5)
}

out_file <- file.path(FIG_DIR, "fig_host_BONvsCLY.pdf")
ggsave(out_file, fig, width = 11, height = 6, dpi = 300)
cat("\nSaved figure:", out_file, "\n")
cat("Saved CSVs to:", CSV_DIR, "\n")
browseURL(out_file)
