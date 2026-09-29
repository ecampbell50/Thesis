# Read coverage over the resident prophage in each strain (Fig S4.4-S4.6).

library(tidyverse)

# Self-contained: resolve relative to this script's location in paper_package
# (…/S4_prophage_coverage/{script,data,figure}). No external paths.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT     <- normalizePath(file.path(.script_dir, ".."))
COV_DIR  <- file.path(ROOT, "data", "coverage_prophage")
FIG_DIR  <- file.path(ROOT, "figure")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Prophage boundaries (for the shaded annotation box)
PROPH <- tribble(
  ~strain, ~proph_start, ~proph_end,
  "C67",   55949,        69028,
  "D32",   5291,         18370,
  "D68",   5291,         18370
)

# ── Load normalization factors ────────────────────────────────────────────────
norm <- read_tsv(file.path(COV_DIR, "normalization_factors.tsv"),
                 show_col_types = FALSE) %>%
  mutate(timepoint = as.integer(timepoint),
         replicate = as.integer(replicate))

# ── Load all region TSVs ──────────────────────────────────────────────────────
files <- list.files(COV_DIR, pattern = "\\.tsv$", full.names = TRUE)
files <- files[!grepl("normalization", files)]

cov_raw <- map_dfr(files, function(f) {
  read_tsv(f, col_names = c("timepoint","strain","condition","replicate",
                            "start","end","depth"),
           col_types = "icciiii")
}) %>%
  mutate(timepoint = as.integer(timepoint),
         replicate = as.integer(replicate),
         # midpoint for plotting
         pos = (start + end) / 2)

# ── Normalise: depth per million bases mapped (like RPM) ─────────────────────
cov <- cov_raw %>%
  left_join(norm %>% select(timepoint, strain, condition, replicate,
                            total_bases_mapped),
            by = c("timepoint","strain","condition","replicate")) %>%
  mutate(depth_norm = depth / total_bases_mapped * 1e6)

# ── Average replicates → mean + ribbon ───────────────────────────────────────
cov_avg <- cov %>%
  group_by(strain, condition, timepoint, pos) %>%
  summarise(
    mean_depth = mean(depth_norm),
    se_depth   = sd(depth_norm) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    timepoint = factor(timepoint, levels = c(0,2,10,20,30,50),
                       labels = paste0(c(0,2,10,20,30,50)," min")),
    condition = factor(condition, levels = c("BON","CLY","CTLR"))
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
COND_COLS <- c(BON = "#1f77b4", CLY = "#d62728", CTLR = "grey50")

for (s in c("C67","D32","D68")) {

  proph_s <- PROPH %>% filter(strain == s)

  p <- ggplot(cov_avg %>% filter(strain == s),
              aes(x = pos, colour = condition, fill = condition)) +

    # Prophage region shading
    annotate("rect",
             xmin = proph_s$proph_start, xmax = proph_s$proph_end,
             ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "orange") +

    # Replicate ribbon (mean ± 1 SE)
    geom_ribbon(aes(ymin = mean_depth - se_depth,
                    ymax = mean_depth + se_depth),
                alpha = 0.18, colour = NA) +

    # Mean line
    geom_line(aes(y = mean_depth), linewidth = 0.5) +

    scale_colour_manual(values = COND_COLS, name = "Condition") +
    scale_fill_manual(values = COND_COLS,   name = "Condition") +

    facet_wrap(~ timepoint, ncol = 3) +
    # After the facet_wrap line, add:
    #coord_cartesian(ylim = c(0, 0.5)) +   # zoom in to see prophage-level expression
    labs(
      title    = paste0(s, " — prophage locus coverage (± 5 kb flanking)"),
      subtitle = "Normalised depth (per million bases mapped). Ribbon = ±1 SE across 3 replicates.\nOrange shading = prophage region.",
      x = "Genomic position",
      y = "Normalised depth"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey93"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(FIG_DIR, paste0("fig_prophage_coverage_", s, ".png")),
         p, width = 10, height = 7, dpi = 300)
  cat("Saved:", s, "\n")
}