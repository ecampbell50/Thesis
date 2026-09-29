# Host gene log2FC vs control across the time course (Fig 4.3).

if (!"tidyverse" %in% installed.packages()) install.packages("tidyverse")
library(tidyverse)

# ── 1. Configuration ──────────────────────────────────────────────────────────
STRAINS    <- c("C67", "D32", "D68")
PHAGES     <- c("BON", "CLY")
TIMEPOINTS <- c(2, 10, 20, 30, 50)   # T=0 excluded (same samples across conditions)

# Self-contained: resolve relative to this script's location in paper_package
# (…/Figure2_host_expression/{script,data,figure}). No external paths.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT           <- normalizePath(file.path(.script_dir, ".."))
DESEQ_DIR      <- file.path(ROOT, "data", "DESeq2_vs_control")
OUT_DIR        <- file.path(ROOT, "figure")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1      # below this: grey; at or above: coloured

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

# ── 3. Load DESeq2 results ────────────────────────────────────────────────────
cat("── Reading DESeq2 CSVs ──────────────────────────────────────────────────\n")

full_df <- bind_rows(lapply(STRAINS, function(s) {
  bind_rows(lapply(PHAGES, function(p) {
    bind_rows(lapply(TIMEPOINTS, function(tp) {
      f <- file.path(DESEQ_DIR, sprintf("DESeq2_%s_%s_vs_CTLR_T%dmin.csv", s, p, tp))
      if (!file.exists(f)) { cat("  [MISSING]", f, "\n"); return(NULL) }
      res           <- read.csv(f, row.names = 1)
      res$gene      <- rownames(res)
      res$strain    <- s
      res$phage     <- p
      res$timepoint <- tp
      res
    }))
  }))
}))

cat(sprintf("  Total rows loaded: %d\n", nrow(full_df)))

# ── 4. Filter and categorise ──────────────────────────────────────────────────
plot_df <- full_df %>%
  filter(!is.na(padj), padj < PADJ_THRESHOLD) %>%
  mutate(
    is_defence  = mapply(function(g, s) g %in% defence_genes[[s]],  gene, strain),
    is_prophage = mapply(function(g, s) g %in% prophage_genes[[s]], gene, strain),
    category = case_when(
      is_defence  ~ "defence",
      is_prophage ~ "prophage",
      abs(log2FoldChange) >= LFC_THRESHOLD ~ "DE",
      TRUE ~ "small"   # sig but |LFC| < 1
    ),
    category    = factor(category, levels = c("small", "DE", "prophage", "defence")),
    phage_label = factor(if_else(phage == "BON", "Bonnie", "Clyde"), levels = c("Bonnie", "Clyde")),
    strain      = factor(strain, levels = STRAINS),
    timepoint   = factor(timepoint, levels = TIMEPOINTS)
  )

cat(sprintf("  Genes passing padj < %g: %d rows\n", PADJ_THRESHOLD, nrow(plot_df)))
plot_df %>%
  dplyr::count(strain, phage, category) %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0) %>%
  print()


# Per-timepoint count label: DE genes only (|LFC| >= LFC_THRESHOLD AND padj < PADJ_THRESHOLD)
tp_counts <- plot_df %>%
  group_by(strain, phage_label, timepoint) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(strain = factor(strain, levels = STRAINS),
         label  = paste0("n = \n ", n))

# ── 5. Plot ───────────────────────────────────────────────────────────────────
pt_colours <- c(
  small    = "grey72",
  DE       = "#ccff00",
  prophage = "#f04a95",
  defence  = "#4066ff"
)
pt_sizes <- c(small = 1, DE = 1.0, prophage = 1, defence = 1)
pt_alpha <- c(small = 0.20, DE = 0.85, prophage = 1, defence = 1)

# Cap y-axis: clips the handful of multi-mapping artefact genes (LFC ~ ±27)
# so they don't stretch the axis and hide everything else.
# Points beyond ±22 are simply not drawn (coord_cartesian, not ylim in scale).
Y_CAP <- 22

fig2 <- ggplot(plot_df,
               aes(x = timepoint, y = log2FoldChange,
                   colour = category, size = category, alpha = category)) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
  #geom_jitter(width = 0.28, height = 0, shape = 19) +
  geom_jitter(data = subset(plot_df, category == "small"), width = 0.28, height = 0, shape = 19) +
  geom_jitter(data = subset(plot_df, category == "DE"), width = 0.35, height = 0, shape = 20) +
  geom_jitter(data = subset(plot_df, category == "prophage"), width = 0.2, height = 0, shape = 17) +
  geom_jitter(data = subset(plot_df, category == "defence"), width = 0.2, height = 0, shape = 15) +

  # n= labels: placed at the top visible break (y = Y_CAP * 0.92) inside the panel
  # so they always show regardless of clip settings
  geom_text(data = tp_counts,
            aes(x = timepoint, y = 7, label = label),
            inherit.aes = FALSE,
            vjust = 0, size = 3.4, colour = "grey30", fontface = "bold") +
  scale_colour_manual(
    values = pt_colours,
    labels = c(small    = paste0("Significant (|log2FC| < ", LFC_THRESHOLD, ")"),
               DE       = paste0("Differentially expressed (|log2FC| ≥ ", LFC_THRESHOLD, ")"),
               prophage = "Prophage gene",
               defence  = "Defence gene"),
    name = NULL,
    # Make the legend keys match the shapes actually plotted (small=circle,
    # DE=circle, prophage=triangle, defence=square) at full opacity and a
    # readable size, rather than faint generic dots.
    guide = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(
        shape = c(19, 20, 17, 15),
        size  = 3.5,
        alpha = 1
      )
    )
  ) +
  scale_size_manual(values  = pt_sizes, guide = "none") +
  scale_alpha_manual(values = pt_alpha, guide = "none") +
  scale_x_discrete(labels = paste0(TIMEPOINTS, " min")) +
  scale_y_continuous(
    trans  = scales::pseudo_log_trans(sigma = 1, base = 2),
    breaks = c(-20, -10, -5, -2, -1, 0, 1, 2, 5),
    labels = c(-20, -10, -5, -2, -1, 0, 1, 2, 5),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  coord_cartesian(ylim = c(-Y_CAP, 7.5), clip = "off") +
  geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
           linetype = "dashed", colour = "black", linewidth = 0.45) +
  facet_grid(phage_label ~ strain) +
  labs(
    title    = "Host gene differential expression during phage infection",
    subtitle = paste0("padj < ", PADJ_THRESHOLD,
                      "  ·  log₂FC (phage vs control)  ·  pseudo-log₂ y-axis"),
    x = "Time post-infection",
    y = expression(log[2]~"Fold Change  (pseudo-log scale)")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5,
                                      margin = margin(b = 4)),
    plot.subtitle      = element_text(size = 9, colour = "grey40", hjust = 0.5,
                                      margin = margin(b = 8)),
    plot.margin        = margin(t = 12, r = 12, b = 8, l = 8),
    strip.text         = element_text(face = "bold", size = 10),
    strip.background   = element_rect(fill = "grey93"),
    legend.position    = "bottom",
    legend.key.size    = unit(1.3, "lines"),
    legend.text        = element_text(size = 11),
    legend.spacing.x   = unit(0.4, "cm"),
    legend.margin      = margin(t = 6),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.spacing.x    = unit(0.9, "lines"),
    panel.spacing.y    = unit(0.7, "lines"),
    axis.text.x        = element_text(size = 8.5),
    axis.text.y        = element_text(size = 8),
    axis.title         = element_text(size = 10)
  )

out_file <- file.path(OUT_DIR, "fig2_host_trajectories.pdf")
ggsave(out_file, fig2, width = 11, height = 12, dpi = 300)
cat("\nSaved:", out_file, "\n")
browseURL(out_file)
