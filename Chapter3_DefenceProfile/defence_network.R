#!/usr/bin/env Rscript
# ==============================================================================
# Chapter 3 - Defence-system protein-similarity network  (Figure 1 + Figure S2)
# ==============================================================================
# Merged from v15_defence_network.R + v15_community_panels.R.
# The plotting code of both is preserved VERBATIM; only the data loading and
# output paths were lifted out so the two share one pass over the inputs.
#
# Run from this directory:   Rscript defence_network.R
#
# Outputs -> figures/network/<system>/
#     <system>_network.pdf/.png    full network for that system
#     <system>_nodes.csv           node metadata
#     community_panels/*.png       one panel per Louvain community
#
# Figure S2 is the "all" run (no filter). Figure 1 is a composite assembled in a
# vector editor from the per-system networks and community panels below.
# ==============================================================================

.pkgs <- c("igraph", "ggraph", "tidyverse", "ggrepel")
.missing <- .pkgs[!.pkgs %in% rownames(installed.packages())]
if (length(.missing)) install.packages(.missing, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(igraph); library(ggraph); library(tidyverse); library(ggrepel)
})

# ---- Inputs ------------------------------------------------------------------
edge_file   <- "data/network_edges.tsv"
size_file   <- "data/node_sizes.tsv"
member_file <- "data/cluster_membership.tsv"
annot_file  <- "data/Ssuis_rawptolemaeaoutput.csv"
OUTROOT     <- "figures/network"

# Systems drawn in Figure 1: every defence type with >= 100 annotated genes.
# NULL is the unfiltered whole network = Figure S2.
# NOTE: filtering is a case-insensitive SUBSTRING match, so "RM" also captures
# composite annotations such as "(p::Other|d::RM)". This is intentional and is
# why the RM panel reports 29,535 proteins rather than 29,521.
SYSTEMS <- list(
  NULL,
  "RM", "Other", "DRT", "Gabija", "RosmerTA", "CRISPR-Cas", "Thoeris",
  "PD-T4-6", "Abi2DF", "Aditi", "CBASS", "SoFIC", "Paris", "Hachiman",
  "AbiH", "Kiwa", "Lamassu", "AbiE", "MazEF", "PD-T7-2"
)

# ---- Shared load: annotations -> per-representative summary ------------------
# (computed once and reused by every system below)
annotations <- read_csv(annot_file, show_col_types = FALSE) %>%
  mutate(across(where(is.character), str_trim)) %>%
  select(protein_id, final_system_type, final_system_subtype)

members <- read_tsv(member_file, show_col_types = FALSE,
                    col_names = c("representative", "member")) %>%
  left_join(annotations, by = c("member" = "protein_id"))

mode_or_unknown <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) "Unknown" else names(sort(table(x), decreasing = TRUE))[1]
}

rep_annot <- members %>%
  group_by(representative) %>%
  summarise(
    majority_type    = mode_or_unknown(final_system_type),
    majority_subtype = mode_or_unknown(final_system_subtype),
    all_types        = paste(unique(na.omit(final_system_type)),    collapse = "|"),
    all_subtypes     = paste(unique(na.omit(final_system_subtype)), collapse = "|"),
    .groups = "drop"
  )

# ==============================================================================
# FUNCTION 1 - full network per system   (body verbatim from v15_defence_network.R)
# ==============================================================================
plot_network <- function(FILTER_TYPES, FILTER_SUBTYPES = NULL) {
  SIZE_RANGE      <- c(1, 30)
  CHARGE          <- 0.9
  NITER           <- 1000
  NODE_SCALE      <- 0.8
  PACK_PADDING    <- 0.6
  LABEL_THRESHOLD <- 10

  filter_tag <- paste(c(FILTER_TYPES, FILTER_SUBTYPES), collapse = "_")
  if (filter_tag == "") filter_tag <- "all"
  outdir <- file.path(OUTROOT, filter_tag)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # ---- Load edges & sizes, apply filters ---------------------------------------
  edges  <- read_tsv(edge_file, show_col_types = FALSE)
  sizes  <- read_tsv(size_file, show_col_types = FALSE)
  nodes  <- sizes %>% left_join(rep_annot, by = "representative")

  match_any <- function(col, patterns) {
    if (is.null(patterns)) return(rep(TRUE, length(col)))
    str_detect(col, regex(paste(patterns, collapse = "|"), ignore_case = TRUE))
  }
  nodes <- nodes %>%
    filter(match_any(all_types, FILTER_TYPES),
           match_any(all_subtypes, FILTER_SUBTYPES))

  if (nrow(nodes) == 0) stop("No nodes remain after filtering.")

  edges <- edges %>% filter(query %in% nodes$representative,
                            target %in% nodes$representative)

  cat(sprintf("After filter: %d nodes, %d edges (%s proteins)\n",
              nrow(nodes), nrow(edges),
              format(sum(nodes$cluster_size), big.mark = ",")))

  # ---- Build graph -------------------------------------------------------------
  if (nrow(edges) == 0) {
    g <- make_empty_graph(directed = FALSE) %>% add_vertices(nrow(nodes))
    V(g)$name <- nodes$representative
    for (col in setdiff(names(nodes), "representative")) {
      vertex_attr(g, col) <- nodes[[col]]
    }
  } else {
    g <- graph_from_data_frame(edges %>% select(query, target),
                               directed = FALSE, vertices = nodes)
    E(g)$weight <- edges$pident
    g <- igraph::simplify(g)
  }

  V(g)$size_scaled <- sqrt(V(g)$cluster_size)

  # ---- Louvain communities -----------------------------------------------------
  if (ecount(g) > 0) {
    set.seed(42)
    V(g)$community <- as.character(igraph::membership(
      igraph::cluster_louvain(g, weights = E(g)$weight)))
  } else {
    V(g)$community <- as.character(seq_len(vcount(g)))
  }

  # ---- Community labels (one row per community) --------------------------------
  comm_labels <- tibble(
    community     = V(g)$community,
    majority_type = V(g)$majority_type,
    majority_sub  = V(g)$majority_subtype
  ) %>%
    group_by(community) %>%
    summarise(
      n_members   = n(),
      top_type    = mode_or_unknown(majority_type[majority_type != "Unknown"]),
      top_subtype = mode_or_unknown(majority_sub [majority_sub  != "Unknown"]),
      .groups = "drop"
    ) %>%
    mutate(label = case_when(
      str_detect(top_type, regex("^RM$", ignore_case = TRUE)) &
        top_subtype != "Unknown"                         ~ top_subtype,
      top_subtype != "Unknown" & top_subtype != top_type ~ paste0(top_type, " (", top_subtype, ")"),
      top_type != "Unknown"                              ~ top_type,
      TRUE                                               ~ "Unknown"
    ))

  # ---- RM-family colour palette ------------------------------------------------
  rm_hues <- c(RM_I = 10, RM_II = 130, RM_IIG = 200, RM_III = 260, RM_IV = 310, RM_V = 50, RM_other = 185)

  rm_family_of <- function(subtype) {
    case_when(
      str_detect(subtype, regex("RM_IV",  ignore_case = TRUE)) ~ "RM_IV",
      str_detect(subtype, regex("RM_III", ignore_case = TRUE)) ~ "RM_III",
      str_detect(subtype, regex("RM_IIG", ignore_case = TRUE)) ~ "RM_IIG",
      str_detect(subtype, regex("RM_II",  ignore_case = TRUE)) ~ "RM_II",
      str_detect(subtype, regex("RM_I",   ignore_case = TRUE)) ~ "RM_I",
      str_detect(subtype, regex("RM_V",   ignore_case = TRUE)) ~ "RM_V",
      str_detect(subtype, regex("RM",     ignore_case = TRUE)) ~ "RM_other",
      TRUE                                                      ~ "Unknown"
    )
  }

  pdc_hues <- c(PDC_M01 = 10,
    PDC_M09 = 20,
    PDC_M11 = 30, 
    PDC_M12 = 40,
    PDC_M14 = 50,
    PDC_M22 = 60,
    PDC_M24 = 70,
    PDC_M37 = 80,
    PDC_M39 = 90,
    PDC_M46 = 100,
    PDC_M47 = 110,
    PDC_M62 = 120,
    PDC_M67 = 130,
    PDC_S02 = 140,
    PDC_S04 = 150,
    PDC_S05 = 160,
    PDC_S07 = 170,
    PDC_S11 = 180,
    PDC_S13 = 190,
    PDC_S14 = 200,
    PDC_S15 = 210,
    PDC_S16 = 220,
    PDC_S18 = 230,
    PDC_S25 = 240,
    PDC_S30 = 250,
    PDC_S32 = 260,
    PDC_S35 = 270,
    PDC_S47 = 280,
    PDC_S49 = 290,
    PDC_S50 = 300,
    PDC_S51 = 310,
    PDC_S53 = 320,
    PDC_S56 = 330,
    PDC_S58 = 340,
    PDC_S60 = 350,
    PDC_S73 = 360
  )
  pdc_family_of <- function(subtype) {
    case_when(
      str_detect(subtype, regex("PDC-M01", ignore_case = TRUE)) ~ "PDC_M01",
      str_detect(subtype, regex("PDC-M09", ignore_case = TRUE)) ~ "PDC_M09",
      str_detect(subtype, regex("PDC-M11", ignore_case = TRUE)) ~ "PDC_M11",
      str_detect(subtype, regex("PDC-M12", ignore_case = TRUE)) ~ "PDC_M12",
      str_detect(subtype, regex("PDC-M14", ignore_case = TRUE)) ~ "PDC_M14",
      str_detect(subtype, regex("PDC-M22", ignore_case = TRUE)) ~ "PDC_M22",
      str_detect(subtype, regex("PDC-M24", ignore_case = TRUE)) ~ "PDC_M24",
      str_detect(subtype, regex("PDC-M37", ignore_case = TRUE)) ~ "PDC_M37",
      str_detect(subtype, regex("PDC-M39", ignore_case = TRUE)) ~ "PDC_M39",
      str_detect(subtype, regex("PDC-M46", ignore_case = TRUE)) ~ "PDC_M46",
      str_detect(subtype, regex("PDC-M47", ignore_case = TRUE)) ~ "PDC_M47",
      str_detect(subtype, regex("PDC-M62", ignore_case = TRUE)) ~ "PDC_M62",
      str_detect(subtype, regex("PDC-M67", ignore_case = TRUE)) ~ "PDC_M67",
      str_detect(subtype, regex("PDC-S02", ignore_case = TRUE)) ~ "PDC_S02",
      str_detect(subtype, regex("PDC-S04", ignore_case = TRUE)) ~ "PDC_S04",
      str_detect(subtype, regex("PDC-S05", ignore_case = TRUE)) ~ "PDC_S05",
      str_detect(subtype, regex("PDC-S07", ignore_case = TRUE)) ~ "PDC_S07",
      str_detect(subtype, regex("PDC-S11", ignore_case = TRUE)) ~ "PDC_S11",
      str_detect(subtype, regex("PDC-S13", ignore_case = TRUE)) ~ "PDC_S13",
      str_detect(subtype, regex("PDC-S14", ignore_case = TRUE)) ~ "PDC_S14",
      str_detect(subtype, regex("PDC-S15", ignore_case = TRUE)) ~ "PDC_S15",
      str_detect(subtype, regex("PDC-S16", ignore_case = TRUE)) ~ "PDC_S16",
      str_detect(subtype, regex("PDC-S18", ignore_case = TRUE)) ~ "PDC_S18",
      str_detect(subtype, regex("PDC-S25", ignore_case = TRUE)) ~ "PDC_S25",
      str_detect(subtype, regex("PDC-S30", ignore_case = TRUE)) ~ "PDC_S30",
      str_detect(subtype, regex("PDC-S32", ignore_case = TRUE)) ~ "PDC_S32",
      str_detect(subtype, regex("PDC-S35", ignore_case = TRUE)) ~ "PDC_S35",
      str_detect(subtype, regex("PDC-S47", ignore_case = TRUE)) ~ "PDC_S47",
      str_detect(subtype, regex("PDC-S49", ignore_case = TRUE)) ~ "PDC_S49",
      str_detect(subtype, regex("PDC-S50", ignore_case = TRUE)) ~ "PDC_S50",
      str_detect(subtype, regex("PDC-S51", ignore_case = TRUE)) ~ "PDC_S51",
      str_detect(subtype, regex("PDC-S53", ignore_case = TRUE)) ~ "PDC_S53",
      str_detect(subtype, regex("PDC-S56", ignore_case = TRUE)) ~ "PDC_S56",
      str_detect(subtype, regex("PDC-S58", ignore_case = TRUE)) ~ "PDC_S58",
      str_detect(subtype, regex("PDC-S60", ignore_case = TRUE)) ~ "PDC_S60",
      str_detect(subtype, regex("PDC-S73", ignore_case = TRUE)) ~ "PDC_S73"
    )
  }

  all_hues   <- c(rm_hues, pdc_hues)
  family_of  <- function(subtype) {
    fam <- pdc_family_of(subtype)
    if_else(is.na(fam), rm_family_of(subtype), fam)
  }

  comm_labels <- comm_labels %>%
    mutate(sys_family = replace_na(family_of(top_subtype), "Unknown")) %>%
    group_by(sys_family) %>%
    mutate(fam_rank = row_number(), fam_n = n()) %>%
    ungroup() %>%
    mutate(
      l_val = 42 + 30 * (fam_rank - 1) / pmax(fam_n - 1, 2),
      colour = if_else(
        sys_family %in% names(all_hues),
        hcl(h = all_hues[sys_family], l = l_val, c = 78),
        hcl(h = 360 * (fam_rank - 1) / pmax(fam_n, 1), l = 60, c = 70)
      )
    )

  community_palette <- setNames(comm_labels$colour, comm_labels$community)

  # ---- Layout: component packing, each sized by sqrt(n_nodes) ------------------
  set.seed(123)
  comp_info <- igraph::components(g)
  comp_xy   <- vector("list", comp_info$no)
  comp_r    <- numeric(comp_info$no)

  for (ci in seq_len(comp_info$no)) {
    idx <- which(comp_info$membership == ci)
    sub <- igraph::induced_subgraph(g, idx)
    n   <- length(idx)

    xy <- if (n == 1) {
      matrix(0, nrow = 1, ncol = 2)
    } else if (ecount(sub) == 0) {
      a <- head(seq(0, 2 * pi, length.out = n + 1), n)
      cbind(cos(a), sin(a)) * NODE_SCALE * sqrt(n)
    } else {
      igraph::layout_with_graphopt(sub, charge = CHARGE, niter = NITER)
    }

    # Centre + rescale to radius NODE_SCALE * sqrt(n)
    xy    <- scale(xy, center = TRUE, scale = FALSE)
    r_raw <- max(sqrt(rowSums(xy^2)))
    target_r <- NODE_SCALE * sqrt(n)
    if (r_raw > 0) xy <- xy * (target_r / r_raw)

    comp_r[ci]   <- max(target_r, NODE_SCALE * 0.5)
    comp_xy[[ci]] <- tibble(name = V(sub)$name, x = xy[, 1], y = xy[, 2], ci = ci)
  }

  # Greedy circle-packing of components (largest at origin, rest closest to centre)
  padded <- comp_r * (1 + PACK_PADDING)
  cx <- numeric(comp_info$no); cy <- numeric(comp_info$no)
  order_idx <- order(comp_r, decreasing = TRUE)

  for (k in seq_along(order_idx)[-1]) {
    ci <- order_idx[k]; ri <- padded[ci]
    best <- c(Inf, 0, 0)
    for (j in seq_len(k - 1)) {
      cj <- order_idx[j]
      for (ang in seq(0, 2 * pi, length.out = 36)) {
        tx <- cx[cj] + (ri + padded[cj]) * cos(ang)
        ty <- cy[cj] + (ri + padded[cj]) * sin(ang)
        placed <- order_idx[seq_len(k - 1)]
        if (all(sqrt((tx - cx[placed])^2 + (ty - cy[placed])^2) >= ri + padded[placed] - 0.01)) {
          d <- sqrt(tx^2 + ty^2)
          if (d < best[1]) best <- c(d, tx, ty)
        }
      }
    }
    cx[ci] <- best[2]; cy[ci] <- best[3]
  }

  coords <- bind_rows(comp_xy) %>% mutate(x = x + cx[ci], y = y + cy[ci])

  layout_df <- create_layout(g, layout = "nicely")
  ord <- match(V(g)$name, coords$name)
  layout_df$x <- coords$x[ord]
  layout_df$y <- coords$y[ord]

  plot_r <- max(sqrt(layout_df$x^2 + layout_df$y^2), na.rm = TRUE) * 1.15

  # ---- Centroids for labels ----------------------------------------------------
  centroids <- layout_df %>%
    as_tibble() %>%
    mutate(community = V(g)$community) %>%
    group_by(community) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
    left_join(comm_labels, by = "community") %>%
    filter(n_members >= LABEL_THRESHOLD)

  # ---- Size legend breaks (include largest node) -------------------------------
  max_cs       <- max(V(g)$cluster_size)
  size_breaks  <- sort(unique(c(1, 5, 25, 100, 500, max_cs)))
  size_labels  <- as.character(size_breaks)

  # ---- Plot --------------------------------------------------------------------
  has_edges <- ecount(g) > 0
  filter_label <- paste(c(FILTER_TYPES, FILTER_SUBTYPES), collapse = ", ")

  p <- ggraph(layout_df) +
    { if (has_edges) geom_edge_link(aes(alpha = weight), colour = "black", width = 0.15, show.legend = FALSE) } +
    geom_node_point(aes(size = size_scaled, colour = community), alpha = 0.9) +
    geom_label_repel(data = centroids,
                     aes(x = x, y = y, label = paste0("C", community, ": ", label, "\n(n=", n_members, ")")),
                     size = 5.2, fontface = "bold", fill = alpha("white", 0.85),
                     label.size = 0.15, label.padding = unit(0.12, "lines"),
                     max.overlaps = 50, seed = 42, min.segment.length = 0,
                     force = 2, force_pull = 0.5) +
    scale_size_continuous(
      name   = "Proteins in cluster",
      range  = SIZE_RANGE,
      breaks = sqrt(size_breaks),
      labels = size_labels,
      guide  = guide_legend(
        title.position = "top",
        label.position = "bottom",
        direction      = "horizontal",
        override.aes   = list(alpha = 0.9, colour = "grey40"),
        nrow           = 1
      )
    ) +
    scale_colour_manual(values = community_palette, guide = "none") +
    { if (has_edges) scale_edge_alpha_continuous(range = c(0.5, 1)) } +
    coord_fixed(ratio = 1, xlim = c(-plot_r, plot_r), ylim = c(-plot_r, plot_r)) +
    labs(
      title = if (filter_label == "") "Defence Network" else paste("Defence Network:", filter_label),
      subtitle = sprintf("%s representatives | %s proteins | %d communities",
                         format(vcount(g), big.mark = ","),
                         format(sum(V(g)$cluster_size), big.mark = ","),
                         n_distinct(V(g)$community)),
      caption = "red=RM_I  green=RM_II  cyan=RM_IIG  blue=RM_III  purple=RM_IV  orange=RM_V  teal=RM_other  grey=unknown"
    ) +
    theme_void(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 18, hjust = 0.5),
          plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
          plot.caption  = element_text(size = 8,  hjust = 0.5, colour = "grey50"),
          plot.margin   = margin(15, 15, 15, 15),
          legend.position = "bottom",
          legend.direction = "horizontal")


  ggsave(file.path(outdir, paste0(filter_tag, "_network.pdf")), p, width = 22, height = 22, dpi = 300)
  ggsave(file.path(outdir, paste0(filter_tag, "_network.png")), p, width = 22, height = 22, dpi = 300)

  # ---- Export node metadata ----------------------------------------------------
  tibble(
    representative   = V(g)$name,
    community        = V(g)$community,
    cluster_size     = V(g)$cluster_size,
    majority_type    = V(g)$majority_type,
    majority_subtype = V(g)$majority_subtype
  ) %>%
    left_join(comm_labels %>% select(community, label, sys_family, n_members),
              by = "community") %>%
    write_csv(file.path(outdir, paste0(filter_tag, "_nodes.csv")))

  cat(sprintf("Done. Output in %s/\n", outdir))
}

# ==============================================================================
# FUNCTION 2 - per-community panels   (body verbatim from v15_community_panels.R)
# ==============================================================================
plot_community_panels <- function(FILTER_TYPES, FILTER_SUBTYPES = NULL) {
  SIZE_RANGE      <- c(3, 60)
  CHARGE          <- 0.9
  NITER           <- 1000
  NODE_SCALE      <- 0.8
  LABEL_THRESHOLD <- 10
  MIN_PANEL_IN <- 6
  MAX_PANEL_IN <- 16
  DPI          <- 300

  filter_tag <- paste(c(FILTER_TYPES, FILTER_SUBTYPES), collapse = "_")
  if (filter_tag == "") filter_tag <- "all"
  outdir_full   <- file.path(OUTROOT, filter_tag)
  outdir_panels <- file.path(outdir_full, "community_panels")
  dir.create(outdir_panels, showWarnings = FALSE, recursive = TRUE)

  match_any <- function(col, patterns) {
    if (is.null(patterns)) return(rep(TRUE, length(col)))
    str_detect(col, regex(paste(patterns, collapse = "|"), ignore_case = TRUE))
  }

  edges <- read_tsv(edge_file, show_col_types = FALSE)
  sizes <- read_tsv(size_file, show_col_types = FALSE)
  nodes <- sizes %>%
    left_join(rep_annot, by = "representative") %>%
    filter(match_any(all_types, FILTER_TYPES),
           match_any(all_subtypes, FILTER_SUBTYPES))

  if (nrow(nodes) == 0) stop("No nodes remain after filtering.")
  edges <- edges %>% filter(query %in% nodes$representative,
                            target %in% nodes$representative)

  cat(sprintf("After filter: %d nodes, %d edges (%s proteins)\n",
              nrow(nodes), nrow(edges),
              format(sum(nodes$cluster_size), big.mark = ",")))

  # ---- Build graph (identical to v15) ------------------------------------------
  if (nrow(edges) == 0) {
    g <- make_empty_graph(directed = FALSE) %>% add_vertices(nrow(nodes))
    V(g)$name <- nodes$representative
    for (col in setdiff(names(nodes), "representative")) {
      vertex_attr(g, col) <- nodes[[col]]
    }
  } else {
    g <- graph_from_data_frame(edges %>% select(query, target),
                               directed = FALSE, vertices = nodes)
    E(g)$weight <- edges$pident
    g <- igraph::simplify(g)
  }
  V(g)$size_scaled <- sqrt(V(g)$cluster_size)

  # ---- Louvain communities (identical to v15) ----------------------------------
  if (ecount(g) > 0) {
    set.seed(42)
    V(g)$community <- as.character(igraph::membership(
      igraph::cluster_louvain(g, weights = E(g)$weight)))
  } else {
    V(g)$community <- as.character(seq_len(vcount(g)))
  }

  # ---- Community labels & palette (identical to v15) ---------------------------
  comm_labels <- tibble(
    community     = V(g)$community,
    majority_type = V(g)$majority_type,
    majority_sub  = V(g)$majority_subtype
  ) %>%
    group_by(community) %>%
    summarise(
      n_members   = n(),
      top_type    = mode_or_unknown(majority_type[majority_type != "Unknown"]),
      top_subtype = mode_or_unknown(majority_sub [majority_sub  != "Unknown"]),
      .groups = "drop"
    ) %>%
    mutate(label = case_when(
      str_detect(top_type, regex("^RM$", ignore_case = TRUE)) &
        top_subtype != "Unknown"                         ~ top_subtype,
      top_subtype != "Unknown" & top_subtype != top_type ~ paste0(top_type, " (", top_subtype, ")"),
      top_type != "Unknown"                              ~ top_type,
      TRUE                                               ~ "Unknown"
    ))

  rm_hues <- c(RM_I = 10, RM_II = 130, RM_IIG = 200, RM_III = 260, RM_IV = 310, RM_V = 50, RM_other = 185)

  rm_family_of <- function(subtype) {
    case_when(
      str_detect(subtype, regex("RM_IV",  ignore_case = TRUE)) ~ "RM_IV",
      str_detect(subtype, regex("RM_III", ignore_case = TRUE)) ~ "RM_III",
      str_detect(subtype, regex("RM_IIG", ignore_case = TRUE)) ~ "RM_IIG",
      str_detect(subtype, regex("RM_II",  ignore_case = TRUE)) ~ "RM_II",
      str_detect(subtype, regex("RM_I",   ignore_case = TRUE)) ~ "RM_I",
      str_detect(subtype, regex("RM_V",   ignore_case = TRUE)) ~ "RM_V",
      str_detect(subtype, regex("RM",     ignore_case = TRUE)) ~ "RM_other",
      TRUE                                                      ~ "Unknown"
    )
  }

  pdc_hues <- c(
    PDC_M01 = 10, PDC_M09 = 20, PDC_M11 = 30, PDC_M12 = 40, PDC_M14 = 50,
    PDC_M22 = 60, PDC_M24 = 70, PDC_M37 = 80, PDC_M39 = 90, PDC_M46 = 100,
    PDC_M47 = 110, PDC_M62 = 120, PDC_M67 = 130, PDC_S02 = 140, PDC_S04 = 150,
    PDC_S05 = 160, PDC_S07 = 170, PDC_S11 = 180, PDC_S13 = 190, PDC_S14 = 200,
    PDC_S15 = 210, PDC_S16 = 220, PDC_S18 = 230, PDC_S25 = 240, PDC_S30 = 250,
    PDC_S32 = 260, PDC_S35 = 270, PDC_S47 = 280, PDC_S49 = 290, PDC_S50 = 300,
    PDC_S51 = 310, PDC_S53 = 320, PDC_S56 = 330, PDC_S58 = 340, PDC_S60 = 350,
    PDC_S73 = 360
  )

  pdc_family_of <- function(subtype) {
    case_when(
      str_detect(subtype, regex("PDC-M01", ignore_case = TRUE)) ~ "PDC_M01",
      str_detect(subtype, regex("PDC-M09", ignore_case = TRUE)) ~ "PDC_M09",
      str_detect(subtype, regex("PDC-M11", ignore_case = TRUE)) ~ "PDC_M11",
      str_detect(subtype, regex("PDC-M12", ignore_case = TRUE)) ~ "PDC_M12",
      str_detect(subtype, regex("PDC-M14", ignore_case = TRUE)) ~ "PDC_M14",
      str_detect(subtype, regex("PDC-M22", ignore_case = TRUE)) ~ "PDC_M22",
      str_detect(subtype, regex("PDC-M24", ignore_case = TRUE)) ~ "PDC_M24",
      str_detect(subtype, regex("PDC-M37", ignore_case = TRUE)) ~ "PDC_M37",
      str_detect(subtype, regex("PDC-M39", ignore_case = TRUE)) ~ "PDC_M39",
      str_detect(subtype, regex("PDC-M46", ignore_case = TRUE)) ~ "PDC_M46",
      str_detect(subtype, regex("PDC-M47", ignore_case = TRUE)) ~ "PDC_M47",
      str_detect(subtype, regex("PDC-M62", ignore_case = TRUE)) ~ "PDC_M62",
      str_detect(subtype, regex("PDC-M67", ignore_case = TRUE)) ~ "PDC_M67",
      str_detect(subtype, regex("PDC-S02", ignore_case = TRUE)) ~ "PDC_S02",
      str_detect(subtype, regex("PDC-S04", ignore_case = TRUE)) ~ "PDC_S04",
      str_detect(subtype, regex("PDC-S05", ignore_case = TRUE)) ~ "PDC_S05",
      str_detect(subtype, regex("PDC-S07", ignore_case = TRUE)) ~ "PDC_S07",
      str_detect(subtype, regex("PDC-S11", ignore_case = TRUE)) ~ "PDC_S11",
      str_detect(subtype, regex("PDC-S13", ignore_case = TRUE)) ~ "PDC_S13",
      str_detect(subtype, regex("PDC-S14", ignore_case = TRUE)) ~ "PDC_S14",
      str_detect(subtype, regex("PDC-S15", ignore_case = TRUE)) ~ "PDC_S15",
      str_detect(subtype, regex("PDC-S16", ignore_case = TRUE)) ~ "PDC_S16",
      str_detect(subtype, regex("PDC-S18", ignore_case = TRUE)) ~ "PDC_S18",
      str_detect(subtype, regex("PDC-S25", ignore_case = TRUE)) ~ "PDC_S25",
      str_detect(subtype, regex("PDC-S30", ignore_case = TRUE)) ~ "PDC_S30",
      str_detect(subtype, regex("PDC-S32", ignore_case = TRUE)) ~ "PDC_S32",
      str_detect(subtype, regex("PDC-S35", ignore_case = TRUE)) ~ "PDC_S35",
      str_detect(subtype, regex("PDC-S47", ignore_case = TRUE)) ~ "PDC_S47",
      str_detect(subtype, regex("PDC-S49", ignore_case = TRUE)) ~ "PDC_S49",
      str_detect(subtype, regex("PDC-S50", ignore_case = TRUE)) ~ "PDC_S50",
      str_detect(subtype, regex("PDC-S51", ignore_case = TRUE)) ~ "PDC_S51",
      str_detect(subtype, regex("PDC-S53", ignore_case = TRUE)) ~ "PDC_S53",
      str_detect(subtype, regex("PDC-S56", ignore_case = TRUE)) ~ "PDC_S56",
      str_detect(subtype, regex("PDC-S58", ignore_case = TRUE)) ~ "PDC_S58",
      str_detect(subtype, regex("PDC-S60", ignore_case = TRUE)) ~ "PDC_S60",
      str_detect(subtype, regex("PDC-S73", ignore_case = TRUE)) ~ "PDC_S73"
    )
  }

  all_hues  <- c(rm_hues, pdc_hues)
  family_of <- function(subtype) {
    fam <- pdc_family_of(subtype)
    if_else(is.na(fam), rm_family_of(subtype), fam)
  }

  comm_labels <- comm_labels %>%
    mutate(sys_family = replace_na(family_of(top_subtype), "Unknown")) %>%
    group_by(sys_family) %>%
    mutate(fam_rank = row_number(), fam_n = n()) %>%
    ungroup() %>%
    mutate(
      l_val  = 42 + 30 * (fam_rank - 1) / pmax(fam_n - 1, 2),
      colour = if_else(
        sys_family %in% names(all_hues),
        hcl(h = all_hues[sys_family], l = l_val, c = 78),
        hcl(h = 360 * (fam_rank - 1) / pmax(fam_n, 1), l = 60, c = 70)
      )
    )

  community_palette <- setNames(comm_labels$colour, comm_labels$community)

  # ---- Shared size-scale breaks (based on the FULL network) --------------------
  # Keeping these fixed across all panels means node sizes are directly comparable.
  max_cs      <- max(V(g)$cluster_size)
  size_breaks <- sort(unique(c(1, 5, 25, 100, 500, max_cs)))
  size_labels <- as.character(size_breaks)

  # ---- Helper: layout for a single (sub)graph ----------------------------------
  layout_subgraph <- function(sg) {
    n <- vcount(sg)
    if (n == 1) {
      xy <- matrix(0, nrow = 1, ncol = 2)
    } else if (ecount(sg) == 0) {
      a  <- head(seq(0, 2 * pi, length.out = n + 1), n)
      xy <- cbind(cos(a), sin(a)) * NODE_SCALE * sqrt(n)
    } else {
      xy <- igraph::layout_with_graphopt(sg, charge = CHARGE, niter = NITER)
    }
    xy <- scale(xy, center = TRUE, scale = FALSE)
    r_raw    <- max(sqrt(rowSums(xy^2)))
    target_r <- NODE_SCALE * sqrt(n)
    if (r_raw > 0) xy <- xy * (target_r / r_raw)
    xy
  }

  # ---- Per-community plots -----------------------------------------------------
  comm_ids <- sort(unique(V(g)$community), method = "radix",
                   decreasing = TRUE)   # largest communities first when sorted by id

  for (cid in comm_ids) {

    # Induce subgraph for this community
    idx <- which(V(g)$community == cid)
    sg  <- igraph::induced_subgraph(g, idx)
    n   <- vcount(sg)

    # Skip singletons
    if (n < 2) {
      cat(sprintf("  Skipping singleton community %s\n", cid))
      next
    }

    meta <- comm_labels %>% filter(community == cid)
    lbl  <- meta$label
    n_members <- meta$n_members

    cat(sprintf("  Community %s: %d nodes (%s)\n", cid, n, lbl))

    # Layout
    set.seed(as.integer(cid) + 100)
    xy <- layout_subgraph(sg)

    # Build ggraph layout manually
    layout_sg <- create_layout(sg, layout = "nicely")
    layout_sg$x <- xy[, 1]
    layout_sg$y <- xy[, 2]

    has_edges <- ecount(sg) > 0
    plot_r    <- max(sqrt(layout_sg$x^2 + layout_sg$y^2), na.rm = TRUE) * 1.2
    if (!is.finite(plot_r) || plot_r == 0) plot_r <- 1

    p_comm <- ggraph(layout_sg) +
      { if (has_edges) geom_edge_link(colour = "black", alpha = 0.15,
                                      width = 0.6, show.legend = FALSE) } +
      geom_node_point(aes(size = size_scaled, colour = community), alpha = 0.9) +
      scale_size_continuous(
        name   = "Proteins in cluster",
        range  = SIZE_RANGE,           # same range as full-network plot
        limits = c(1, sqrt(max_cs)),   # fix to full-network range so breaks aren't dropped
        breaks = sqrt(size_breaks),    # same breaks, so legend is identical
        labels = size_labels,
        guide  = guide_legend(
          title.position = "top",
          label.position = "bottom",
          direction      = "horizontal",
          override.aes   = list(alpha = 0.9, colour = "grey40"),
          nrow           = 1
        )
      ) +
      scale_colour_manual(values = community_palette, guide = "none") +
      coord_fixed(ratio = 1, xlim = c(-plot_r, plot_r), ylim = c(-plot_r, plot_r)) +
      theme_void(base_size = 14) +
      theme(
        plot.margin      = margin(10, 10, 10, 10),
        legend.position  = "none"
      )

    # Panel size proportional to community, bounded
    panel_in <- clamp <- function(x, lo, hi) max(lo, min(hi, x))
    panel_in <- clamp(MIN_PANEL_IN * sqrt(n / 5), MIN_PANEL_IN, MAX_PANEL_IN)

    # Save into a subdirectory named after the community's top subtype
    safe_lbl  <- gsub("[^A-Za-z0-9_-]", "_", lbl)
    type_dir  <- file.path(outdir_panels, safe_lbl)
    dir.create(type_dir, showWarnings = FALSE, recursive = TRUE)
    fname     <- sprintf("community_%s_%s.png", cid, safe_lbl)
    out_path  <- file.path(type_dir, fname)

    ggsave(out_path, p_comm, width = panel_in, height = panel_in, dpi = DPI)
  }

  cat(sprintf("\nDone. %d community PNGs saved to %s/\n",
              length(comm_ids), outdir_panels))
}

# ==============================================================================
# MAIN
# ==============================================================================
for (sys in SYSTEMS) {
  tag <- if (is.null(sys)) "all (Figure S2)" else sys
  cat("\n========================================\n")
  cat("  ", tag, "\n", sep = "")
  cat("========================================\n")
  ok <- tryCatch({ plot_network(sys); TRUE },
                 error = function(e) { cat("  ! network failed: ", conditionMessage(e), "\n"); FALSE })
  if (ok) tryCatch(plot_community_panels(sys),
                   error = function(e) cat("  ! panels failed: ", conditionMessage(e), "\n"))
}
cat("\nDone. Outputs in ", OUTROOT, "/\n", sep = "")
