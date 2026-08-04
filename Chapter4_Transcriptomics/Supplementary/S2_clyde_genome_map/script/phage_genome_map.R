# phage_genome_map.R
# Circular genome map for Bonnie and Clyde phages from GFF3 annotation.
# Forward-strand genes on outer ring, reverse-strand on inner ring, as arrows.
#
# Run from: 2_Analysis/
#   source("scripts/R/phage_genome_map.R")

library(tidyverse)

# Self-contained: resolve relative to this script's location in paper_package
# (…/{script,data,figure}). GFF3 inputs live in data/. No external paths.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT <- normalizePath(file.path(.script_dir, ".."))
DATA <- file.path(ROOT, "data")
FIG  <- file.path(ROOT, "figure")

# ── Configuration — edit paths here ──────────────────────────────────────────
PHAGE_CONFIGS <- list(
  list(
    name       = "Clyde",
    gff        = file.path(DATA, "Clyde_allannos.gff3"),
    genome_len = 34734,
    out        = file.path(FIG, "clyde_genome_map.pdf")
  ),
  list(
    name       = "Bonnie",
    gff        = file.path(DATA, "Bonnie_allannos.gff3"),
    genome_len = NA,   # leave NA to auto-read from ##sequence-region header
    out        = file.path(FIG, "bonnie_genome_map.pdf")
  )
)

dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# ── Colour scheme ─────────────────────────────────────────────────────────────
CAT_COLS <- c(
  "Structural"     = "#2196F3",
  "DNA metabolism" = "#FF9800",
  "Lysis"          = "#F44336",
  "Integration"    = "#9C27B0",
  "Regulation"     = "#4CAF50",
  "Defence"        = "#E91E63",
  "Hypothetical"   = "grey78",
  "Other"          = "#795548"
)

# ── Functional categorisation ─────────────────────────────────────────────────
categorise <- function(p) {
  case_when(
    str_detect(p, regex("hypoth|unknown|uncharacteri", ignore_case = TRUE))
      ~ "Hypothetical",
    str_detect(p, regex("tail|head|capsid|portal|terminase|baseplate|fiber|tape measure|Dit|neck|connector|decorator|virion|structural|coat|spike|receptor|adsorption|major capsid|prohead|scaffold", ignore_case = TRUE))
      ~ "Structural",
    str_detect(p, regex("DNA|polyme|helica|primase|replicat|exonuclease|recombinase|RecT|annealing|topoisomerase|ligase|nuclease|acetyltransferase|Holliday|methyltransferase", ignore_case = TRUE))
      ~ "DNA metabolism",
    str_detect(p, regex("lysin|holin|lysis|spanin|endolysin", ignore_case = TRUE))
      ~ "Lysis",
    str_detect(p, regex("integrase|excisionase|transposase|recombination directionality", ignore_case = TRUE))
      ~ "Integration",
    str_detect(p, regex("transcri|repressor|regulator|anti-repressor|binding|activator|sigma", ignore_case = TRUE))
      ~ "Regulation",
    str_detect(p, regex("defense|DefenseFinder|restriction|toxin|antitoxin|immunity|anti-phage|abortive", ignore_case = TRUE))
      ~ "Defence",
    TRUE ~ "Other"
  )
}

# ── Arrow polygon builder ─────────────────────────────────────────────────────
# Returns a data frame of polygon vertices for one gene arrow.
# + strand: arrow tip at deg_end (clockwise)
# - strand: arrow tip at deg_start (counter-clockwise)
gene_arrow <- function(deg_start, deg_end, r_lo, r_hi, strand,
                       category, locus, product) {
  span  <- deg_end - deg_start
  if (span <= 0) return(NULL)
  notch <- min(span * 0.3, 2.5)        # arrowhead depth in degrees
  r_mid <- (r_lo + r_hi) / 2
  npts  <- max(5, ceiling(span * 4))   # arc interpolation steps

  if (strand == "+") {
    body_end <- max(deg_start + span * 0.05, deg_end - notch)
    pts <- bind_rows(
      tibble(angle = seq(deg_start, body_end, length.out = npts), r = r_lo),
      tibble(angle = deg_end, r = r_mid),
      tibble(angle = seq(body_end, deg_start, length.out = npts), r = r_hi)
    )
  } else {
    body_start <- min(deg_end - span * 0.05, deg_start + notch)
    pts <- bind_rows(
      tibble(angle = deg_start, r = r_mid),
      tibble(angle = seq(body_start, deg_end, length.out = npts), r = r_lo),
      tibble(angle = seq(deg_end, body_start, length.out = npts), r = r_hi)
    )
  }
  pts %>% mutate(category = category, locus = locus, product = product,
                 strand = strand)
}

# ── Main plot function ────────────────────────────────────────────────────────
make_genome_map <- function(cfg) {

  # Auto-read genome length from GFF header if not supplied
  genome_len <- cfg$genome_len
  if (is.na(genome_len)) {
    hdr        <- readLines(path.expand(cfg$gff), n = 10)
    seq_line   <- hdr[grepl("##sequence-region", hdr)][1]
    genome_len <- as.integer(str_extract(seq_line, "\\d+$"))
    message(cfg$name, ": genome length read as ", genome_len, " bp")
  }

  # Parse GFF3
  gff <- read_tsv(path.expand(cfg$gff), comment = "#",
                  col_names = c("seqid","source","type","start","end",
                                "score","strand","phase","attributes"),
                  show_col_types = FALSE) %>%
    filter(type == "CDS") %>%
    mutate(
      locus     = str_extract(attributes, "(?<=locus_tag=)[^;]+"),
      product   = coalesce(str_extract(attributes, "(?<=product=)[^;]+"),
                           "hypothetical protein"),
      category  = factor(categorise(product), levels = names(CAT_COLS)),
      deg_start = (start - 1) / genome_len * 360,
      deg_end   =  end        / genome_len * 360,
      r_lo      = if_else(strand == "+", 0.82, 0.62),
      r_hi      = if_else(strand == "+", 0.94, 0.74)
    )

  # Build arrow polygons
  arcs <- pmap_dfr(
    gff %>% select(deg_start, deg_end, r_lo, r_hi, strand,
                   category, locus, product),
    gene_arrow
  ) %>% mutate(category = factor(category, levels = names(CAT_COLS)))

  # Callout labels for all non-hypothetical genes
  labels <- gff %>%
    filter(!str_detect(product, regex("hypoth", ignore_case = TRUE))) %>%
    mutate(
      angle_mid   = (deg_start + deg_end) / 2,
      short_label = str_trunc(product, 24, ellipsis = ".."),
      # Rotate text radially; flip on left side so it stays readable
      text_angle  = case_when(
        angle_mid <= 180 ~ 90  - angle_mid,
        TRUE             ~ 270 - angle_mid
      ),
      text_hjust  = if_else(angle_mid <= 180, 0, 1)
    )

  # Tick marks every 5 kb
  tick_at <- seq(0, genome_len - 1, by = 5000)
  ticks <- tibble(
    angle = tick_at / genome_len * 360,
    label = paste0(tick_at / 1000, " kb")
  )

  # ── Build plot ──────────────────────────────────────────────────────────────
  ggplot() +

    # Genome backbone circles (two thin rings)
    geom_rect(aes(xmin = 0, xmax = 360, ymin = 0.796, ymax = 0.806),
              fill = "grey60", colour = NA) +
    geom_rect(aes(xmin = 0, xmax = 360, ymin = 0.594, ymax = 0.606),
              fill = "grey75", colour = NA) +

    # Gene arrows with thin black border
    geom_polygon(data  = arcs,
                 aes(x = angle, y = r, group = locus, fill = category),
                 colour = "black", linewidth = 0.15, alpha = 0.95) +

    # Callout lines from just outside outer ring to label start
    geom_segment(data = labels,
                 aes(x = angle_mid, xend = angle_mid,
                     y = 0.955, yend = 1.07),
                 linewidth = 0.2, colour = "grey55") +

    # Gene product labels (radial, readable on both sides)
    geom_text(data  = labels,
              aes(x = angle_mid, y = 1.10, label = short_label,
                  angle = text_angle, hjust = text_hjust),
              size = 2.2, colour = "grey15") +

    # Tick marks: thin rects on the backbone circle — geom_rect is more
    # reliable than geom_segment in coord_polar (segments can drift off-angle).
    # 0.4 degrees either side gives a thin but visible mark.
    geom_rect(data = ticks,
              aes(xmin = angle - 0.4, xmax = angle + 0.4,
                  ymin = 0.790, ymax = 0.816),
              fill = "grey25", colour = NA) +

    # Tick labels inside the ring, in the gap between inner gene ring (0.74)
    # and backbone (0.790), to keep the outside clear for gene callouts.
    geom_text(data = ticks,
              aes(x = angle, y = 0.77, label = label),
              size = 2.3, colour = "grey40") +

    scale_fill_manual(values = CAT_COLS, name = "Function", drop = FALSE) +
    coord_polar(theta = "x", start = 0, direction = 1) +
    scale_x_continuous(limits = c(0, 360), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.55), expand = c(0, 0)) +
    labs(
      title    = paste0(cfg$name, " phage genome"),
      subtitle = paste0(format(genome_len, big.mark = ","), " bp  ·  ",
                        nrow(gff), " CDS")
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 14,
                                      margin = margin(b = 2)),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey50", size = 9,
                                      margin = margin(b = 8)),
      legend.position  = "right",
      legend.key.size  = unit(0.5, "cm"),
      legend.text      = element_text(size = 9),
      legend.title     = element_text(size = 10, face = "bold"),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

# ── Generate maps for all configured phages ───────────────────────────────────
for (cfg in PHAGE_CONFIGS) {
  gff_path <- path.expand(cfg$gff)
  if (!file.exists(gff_path)) {
    message("Skipping ", cfg$name, ": GFF not found at ", gff_path)
    next
  }
  p <- make_genome_map(cfg)
  ggsave(cfg$out, p, width = 10, height = 10, dpi = 300)
  cat("Saved:", cfg$out, "\n")
}
