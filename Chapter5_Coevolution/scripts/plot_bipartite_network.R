#!/usr/bin/env Rscript
# Bacteria-prophage network (Fig 5.8).

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(ggnewscale)
})
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
in_dir   <- "."
out_dir  <- "figures"

edge_file <- file.path(in_dir, "data/Bipartite_BAC_PRO_genome_edgetable.csv")
node_file <- file.path(in_dir, "data/GenomeType_nodetable.csv")
bac_map   <- file.path(in_dir, "data/genome_community_mapping_49_BAC.csv")
pro_map   <- file.path(in_dir, "data/genome_community_mapping_44_PRO.csv")

seed         <- 42
top_k_each   <- 15     # full figure: colour the N largest BAC communities and N largest PRO communities
fr_niter     <- 600    # iterations for Fruchterman-Reingold layouts
region_gap   <- 0.5    # horizontal gap between the BAC and PRO regions (full figure)
clump_spread <- 0.004  # full figure: how far to inflate each community disc. Kept small ON PURPOSE:
                       # compact communities let the host links stack into legible bundles instead
                       # of fanning out into a red haze. Bigger = more spread but fuzzier connections.
node_bg_size <- 1.0    # point size for unclustered / background genomes
node_top_size<- 2.4    # point size for the coloured top-15 community genomes
rim_frac     <- 0.80   # full figure: radius (as fraction of region) where the coloured top-15 communities sit
core_frac    <- 0.52   # full figure: radius of the central grey core (unclustered / non-top-15 genomes)

# Community summary thresholds / options
min_summary_sz   <- 5     # only treat communities with >= this many genomes as "real"
min_summary_edge <- 3     # only draw between-community links backed by >= this many edges
drop_isolated    <- TRUE  # drop summary nodes with no qualifying between-community link
show_sim_edges   <- TRUE  # show grey BAC-BAC / PRO-PRO between-community links (FALSE = host only)
main_comp_only   <- TRUE  # summary: show only the main connected component (drops isolated pairs -> fills frame)
host_bundle_min  <- 3     # full figure: draw one red line per BAC-community<->PRO-community host coupling
                          #   backed by >= this many host links (centroid-to-centroid, width ~ count).
                          #   Replaces 3,216 individual host lines (which fan into a red haze).

full_png    <- file.path(out_dir, "bipartite_full_network.png")
summary_png <- file.path(out_dir, "bipartite_community_summary.png")
layout_csv  <- file.path(out_dir, "data/full_network_layout.csv")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

# Two 15-colour qualitative palettes (BAC vs PRO get their own legend)
bac_pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF",
             "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#1B9E77","#E7298A")
pro_pal <- c("#1F78B4","#33A02C","#6A3D9A","#B15928","#D7263D","#FF7F00","#2166AC",
             "#762A83","#01665E","#C51B7D","#4D9221","#8C510A","#5E3C99","#E08214","#7B3294")

# ----------------------------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------------------------
message("Loading data ...")
edges <- fread(edge_file); setnames(edges, c("Source", "Target", "Value"))
nodes <- fread(node_file); setnames(nodes, c("Genome", "Type"))
bac <- fread(bac_map); setnames(bac, c("Genome", "Community")); bac[, comm := paste0("BAC_com_", Community)]
pro <- fread(pro_map); setnames(pro, c("Genome", "Community")); pro[, comm := paste0("PRO_com_", Community)]
comm_map <- rbindlist(list(bac[, .(Genome, comm)], pro[, .(Genome, comm)]))

node_attr <- merge(nodes, comm_map, by = "Genome", all.x = TRUE)
node_attr[is.na(comm), comm := "unassigned"]

# ----------------------------------------------------------------------------
# 2. Classify edges and build the master graph
# ----------------------------------------------------------------------------
type_lookup <- setNames(node_attr$Type, node_attr$Genome)
edges[, type_s := type_lookup[Source]]
edges[, type_t := type_lookup[Target]]
edges[, edge_type := fifelse(type_s == "BAC" & type_t == "BAC", "BAC-BAC",
                      fifelse(type_s == "PRO" & type_t == "PRO", "PRO-PRO", "host"))]

g <- graph_from_data_frame(d = edges[, .(Source, Target, weight = Value, edge_type)],
                           vertices = node_attr[, .(Genome, Type, comm)],
                           directed = FALSE)

# ---- Community census (answers "where are the others?") --------------------
census <- function(tag) {
  s <- node_attr[Type == tag & comm != "unassigned", .N, by = comm]$N
  message(sprintf("  %s: %d community labels | size 1: %d, size 2: %d, size 3-4: %d, size 5-10: %d, size >10: %d  (real clusters size>=%d: %d)",
                  tag, length(s), sum(s==1), sum(s==2), sum(s %in% 3:4),
                  sum(s %in% 5:10), sum(s>10), min_summary_sz, sum(s>=min_summary_sz)))
}
message(sprintf("Nodes: %d (BAC %d, PRO %d) | Edges: %d (BAC-BAC %d, PRO-PRO %d, host %d, all host Value==1)",
                nrow(node_attr), sum(node_attr$Type=="BAC"), sum(node_attr$Type=="PRO"),
                nrow(edges), sum(edges$edge_type=="BAC-BAC"),
                sum(edges$edge_type=="PRO-PRO"), sum(edges$edge_type=="host")))
message("Community census:")
census("BAC"); census("PRO")

# ----------------------------------------------------------------------------
# 3. Per-set top-N community colour schemes (full figure)
#    BAC and PRO get independent palettes + independent legends.
# ----------------------------------------------------------------------------
comm_sizes <- node_attr[comm != "unassigned", .N, by = comm][order(-N)]
bac_top <- head(comm_sizes[grepl("^BAC", comm), comm], top_k_each)
pro_top <- head(comm_sizes[grepl("^PRO", comm), comm], top_k_each)

node_attr[, bac_fill := NA_character_]
node_attr[Type == "BAC", bac_fill := fifelse(comm %in% bac_top, comm, "Other bacteria")]
node_attr[, pro_fill := NA_character_]
node_attr[Type == "PRO", pro_fill := fifelse(comm %in% pro_top, comm, "Other prophage")]

bac_levels <- c(bac_top, "Other bacteria")
pro_levels <- c(pro_top, "Other prophage")
node_attr[, bac_fill := factor(bac_fill, levels = bac_levels)]
node_attr[, pro_fill := factor(pro_fill, levels = pro_levels)]
bac_fill_vals <- setNames(c(bac_pal[seq_along(bac_top)], "grey80"), bac_levels)
pro_fill_vals <- setNames(c(pro_pal[seq_along(pro_top)], "grey80"), pro_levels)

# ----------------------------------------------------------------------------
# 4. Bipartite full-network layout
#    Host edges (Value==1) are NOT used for positioning - otherwise every
#    prophage collapses onto its host. Instead: lay out the BAC subgraph and
#    the PRO subgraph independently (each component packed to use space), then
#    place BAC on the left and PRO on the right; host links bridge the two.
# ----------------------------------------------------------------------------
# Scale a coord matrix into a disc of radius h (preserving aspect). Scaling uses
# the q-th percentile radius so the dense core fills the disc; the few outlier
# singletons that FR flings out are clamped onto the rim (a natural sparse halo)
# instead of stretching the canvas and emptying the centre.
fit_box <- function(m, h, q = 0.95) {
  m <- sweep(m, 2, colMeans(m))
  r <- sqrt(rowSums(m^2))
  s <- as.numeric(quantile(r, q)); if (!is.finite(s) || s == 0) s <- max(r, 1)
  m <- m / s * h
  r2 <- sqrt(rowSums(m^2)); over <- r2 > h & r2 > 0
  m[over, ] <- m[over, ] * (h / r2[over])
  m
}

if (file.exists(layout_csv)) {
  message("Loading cached layout from ", layout_csv)
  lay_dt <- fread(layout_csv)
} else {
  message("Computing bipartite layout (Fruchterman-Reingold per set) ...")
  g_bac <- induced_subgraph(g, V(g)[V(g)$Type == "BAC"])  # keeps only BAC-BAC edges
  g_pro <- induced_subgraph(g, V(g)[V(g)$Type == "PRO"])  # keeps only PRO-PRO edges

  # Plain FR (not component-packing): the connected core fills the centre and
  # disconnected singletons are pushed outward as a halo -> no hollow "wheel".
  set.seed(seed)
  lay_bac <- layout_with_fr(g_bac, niter = fr_niter)
  set.seed(seed)
  lay_pro <- layout_with_fr(g_pro, niter = fr_niter)

  # keep per-node density comparable: box half-width ~ sqrt(n), normalised so the larger set = 1
  sb <- sqrt(vcount(g_bac)); sp <- sqrt(vcount(g_pro))
  h_bac <- sb / max(sb, sp); h_pro <- sp / max(sb, sp)

  lay_bac <- fit_box(lay_bac, h_bac)
  lay_pro <- fit_box(lay_pro, h_pro)
  lay_bac[, 1] <- lay_bac[, 1] - (h_bac + region_gap/2)   # shift BAC left
  lay_pro[, 1] <- lay_pro[, 1] + (h_pro + region_gap/2)   # shift PRO right

  lay_dt <- rbind(
    data.table(Genome = V(g_bac)$name, x = lay_bac[, 1], y = lay_bac[, 2]),
    data.table(Genome = V(g_pro)$name, x = lay_pro[, 1], y = lay_pro[, 2]))
  fwrite(lay_dt, layout_csv)
  message("  Layout cached to ", layout_csv, " (delete to recompute)")
}

# --- Core/rim re-arrangement ----------------------------------------------
# Rebuild each region so the unclustered / non-top-15 genomes form a compact
# grey CORE in the centre and the coloured top-15 communities sit around the
# RIM. Because the host-coupling lines attach to community centroids, pushing
# the big communities outward moves their endpoints to the periphery, so the
# thick links run rim-to-rim instead of all piling through the centre - much
# easier to trace which prophage cluster couples to which bacterial cluster.
# This is a post-layout transform (works off the cached FR coords).
lay_dt <- merge(lay_dt, node_attr[, .(Genome, comm, Type)], by = "Genome", sort = FALSE)

# same per-region scale rule as the layout step, recomputed so this works off cache
sb <- sqrt(sum(node_attr$Type == "BAC")); sp <- sqrt(sum(node_attr$Type == "PRO"))
h_bac <- sb / max(sb, sp); h_pro <- sp / max(sb, sp)

arrange_region <- function(d, top_comms, h, cx0) {
  d <- copy(d)
  d[, `:=`(rx = x - mean(x), ry = y - mean(y))]    # local coords, centred on FR centroid
  d[, is_top := comm %in% top_comms]

  # CORE: grey (non-top-15) genomes, FR cloud rescaled into the central disc
  core <- d[is_top == FALSE]
  if (nrow(core)) {
    rr <- sqrt(core$rx^2 + core$ry^2)
    s  <- as.numeric(quantile(rr, 0.95)); if (!is.finite(s) || s == 0) s <- max(rr, 1)
    core[, `:=`(nx = rx / s * (h * core_frac), ny = ry / s * (h * core_frac))]
    r2 <- sqrt(core$nx^2 + core$ny^2); ov <- r2 > h * core_frac & r2 > 0
    core[ov, `:=`(nx = nx * (h * core_frac / r2[ov]), ny = ny * (h * core_frac / r2[ov]))]
  }

  # RIM: each top-15 community gets an angular slot (width ~ sqrt(size)) on the
  # rim ring, members jittered into a compact blob around that point.
  rim <- d[is_top == TRUE]
  if (nrow(rim)) {
    cc <- rim[, .(a = atan2(mean(ry), mean(rx)), sz = .N), by = comm]
    setorder(cc, a)                                # keep FR angular neighbours adjacent
    cc[, w := sqrt(sz)]
    cc[, theta := 2 * pi * (cumsum(w) - w / 2) / sum(w)]
    rim <- merge(rim, cc[, .(comm, theta, sz)], by = "comm")
    set.seed(seed)
    rim[, jr  := clump_spread * sqrt(sz) * sqrt(runif(.N))]
    rim[, jth := runif(.N, 0, 2 * pi)]
    rim[, `:=`(nx = h * rim_frac * cos(theta) + jr * cos(jth),
               ny = h * rim_frac * sin(theta) + jr * sin(jth))]
  }
  out <- rbindlist(list(core[, .(Genome, nx, ny)], rim[, .(Genome, nx, ny)]))
  out[, .(Genome, x = nx + cx0, y = ny)]
}

lay_dt <- rbindlist(list(
  arrange_region(lay_dt[Type == "BAC"], bac_top, h_bac, -(h_bac + region_gap / 2)),
  arrange_region(lay_dt[Type == "PRO"], pro_top, h_pro,  (h_pro + region_gap / 2))))

ndf  <- merge(lay_dt, node_attr, by = "Genome", sort = FALSE)
xs <- setNames(lay_dt$x, lay_dt$Genome); ys <- setNames(lay_dt$y, lay_dt$Genome)

# within-region similarity edges: still drawn individually (they define the blobs)
edf_sim <- edges[edge_type != "host",
                 .(x = xs[Source], y = ys[Source], xend = xs[Target], yend = ys[Target])]

# host edges: aggregate to ONE line per (BAC community <-> PRO community) pair,
# drawn centroid-to-centroid with width ~ number of host links. This is what
# makes the biology legible: distinct prophage clusters coupling to distinct
# bacterial clusters, instead of 3,216 individual lines fanning into a haze.
cent  <- ndf[comm != "unassigned", .(cx = mean(x), cy = mean(y)), by = comm]
cmx   <- setNames(cent$cx, cent$comm); cmy <- setNames(cent$cy, cent$comm)
hoste <- edges[edge_type == "host"]
hoste[, cs := comm_map[match(Source, Genome), comm]][, ct := comm_map[match(Target, Genome), comm]]
hoste[, `:=`(bac_c = fifelse(type_lookup[Source] == "BAC", cs, ct),
             pro_c = fifelse(type_lookup[Source] == "PRO", cs, ct))]
hoste <- hoste[!is.na(bac_c) & !is.na(pro_c) & bac_c != "unassigned" & pro_c != "unassigned"]
edf_host <- hoste[, .(n_links = .N), by = .(bac_c, pro_c)][n_links >= host_bundle_min]
edf_host[, `:=`(x = cmx[bac_c], y = cmy[bac_c], xend = cmx[pro_c], yend = cmy[pro_c])]
edf_host <- edf_host[is.finite(x) & is.finite(xend)]
message(sprintf("  host couplings: %d community-pairs with >=%d links (from %d individual host edges)",
                nrow(edf_host), host_bundle_min, nrow(hoste)))

# ----------------------------------------------------------------------------
# 5. Plot full network (two fill scales -> two legends, placed beside each set)
# ----------------------------------------------------------------------------
message("Rendering full network figure ...")
xr <- range(ndf$x)
ndf_bac <- ndf[Type == "BAC"]; ndf_pro <- ndf[Type == "PRO"]

# main panel, no legends (legends are extracted separately and placed L/R)
p_main <- ggplot() +
  geom_segment(data = edf_host, aes(x, y, xend = xend, yend = yend,
                                    linewidth = n_links, alpha = n_links),
               colour = "#D7263D", lineend = "round") +
  scale_linewidth_continuous(range = c(0.2, 2.8), guide = "none") +
  scale_alpha_continuous(range = c(0.18, 0.85), guide = "none") +
  # faint grey background = unclustered / non-top-15 genomes (recede)
  geom_point(data = ndf_bac[bac_fill == "Other bacteria"], aes(x, y),
             shape = 21, size = node_bg_size, stroke = 0, fill = "grey70", alpha = 0.45) +
  geom_point(data = ndf_pro[pro_fill == "Other prophage"], aes(x, y),
             shape = 24, size = node_bg_size, stroke = 0, fill = "grey70", alpha = 0.45) +
  # top-15 communities pop on top
  geom_point(data = ndf_bac[bac_fill %in% bac_top], aes(x, y, fill = bac_fill),
             shape = 21, size = node_top_size, stroke = 0.15, colour = "grey25") +
  scale_fill_manual(values = bac_fill_vals, limits = bac_top, na.translate = FALSE) +
  new_scale_fill() +
  geom_point(data = ndf_pro[pro_fill %in% pro_top], aes(x, y, fill = pro_fill),
             shape = 24, size = node_top_size, stroke = 0.15, colour = "grey25") +
  scale_fill_manual(values = pro_fill_vals, limits = pro_top, na.translate = FALSE) +
  annotate("text", x = xr[1], y = max(ndf$y), label = "BACTERIA",
           hjust = 0, fontface = "bold", colour = "grey30", size = 5) +
  annotate("text", x = xr[2], y = max(ndf$y), label = "PROPHAGES",
           hjust = 1, fontface = "bold", colour = "grey30", size = 5) +
  coord_equal() +
  theme_void(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "S. suis bacteria-prophage bipartite similarity network",
       subtitle = "Bacteria (left) and prophages (right) laid out separately; red links join coupled BAC<->PRO communities (width ~ no. of host links)")

# helper: pull the legend grob out of a plot
extract_legend <- function(p) {
  g <- ggplotGrob(p)
  i <- which(grepl("guide-box", g$layout$name))
  for (k in i) if (!inherits(g$grobs[[k]], "zeroGrob")) return(g$grobs[[k]])
  g$grobs[[i[1]]]
}
leg_theme <- theme_void(base_size = 12) +
  theme(legend.position = "right", legend.key.size = unit(4, "mm"),
        legend.title = element_text(face = "bold"))
bac_leg <- extract_legend(
  ggplot(ndf_bac[bac_fill %in% bac_top], aes(x, y, fill = bac_fill)) + geom_point(shape = 21, size = 3) +
    scale_fill_manual(values = bac_fill_vals, limits = bac_top, na.translate = FALSE,
                      name = paste0("Bacterial community\n(top ", top_k_each, " by size)")) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 3.5), ncol = 1)) + leg_theme)
pro_leg <- extract_legend(
  ggplot(ndf_pro[pro_fill %in% pro_top], aes(x, y, fill = pro_fill)) + geom_point(shape = 24, size = 3) +
    scale_fill_manual(values = pro_fill_vals, limits = pro_top, na.translate = FALSE,
                      name = paste0("Prophage community\n(top ", top_k_each, " by size)")) +
    guides(fill = guide_legend(override.aes = list(shape = 24, size = 3.5), ncol = 1)) + leg_theme)

p_full <- cowplot::plot_grid(bac_leg, p_main, pro_leg, nrow = 1,
                             rel_widths = c(0.15, 0.70, 0.15))
ggsave(full_png, p_full, width = 18, height = 9, dpi = 320, bg = "white", limitsize = FALSE)
message("  Saved ", full_png)

# ----------------------------------------------------------------------------
# 6. Community summary network
# ----------------------------------------------------------------------------
message("Building community summary ...")
comm_lookup <- setNames(node_attr$comm, node_attr$Genome)
ctype <- function(cm) fifelse(grepl("^BAC", cm), "BAC", "PRO")

agg <- copy(edges)
agg[, c1 := comm_lookup[Source]][, c2 := comm_lookup[Target]]
agg <- agg[c1 != "unassigned" & c2 != "unassigned"]
agg[, `:=`(ca = pmin(c1, c2), cb = pmax(c1, c2))]
agg_inter <- agg[ca != cb, .(n_edges = .N, w = sum(Value)), by = .(ca, cb)]
agg_inter[, edge_type := fifelse(ctype(ca) == "BAC" & ctype(cb) == "BAC", "BAC-BAC",
                          fifelse(ctype(ca) == "PRO" & ctype(cb) == "PRO", "PRO-PRO", "host"))]

comm_nodes <- node_attr[comm != "unassigned", .(size = .N), by = comm]
comm_nodes[, Type := ctype(comm)]
fwrite(comm_nodes, file.path(out_dir, "data/community_summary_nodes.csv"))
fwrite(agg_inter,  file.path(out_dir, "data/community_summary_edges.csv"))

# Apply thresholds
keep_nodes <- comm_nodes[size >= min_summary_sz, comm]
agg_show   <- agg_inter[ca %in% keep_nodes & cb %in% keep_nodes & n_edges >= min_summary_edge]
if (!show_sim_edges) agg_show <- agg_show[edge_type == "host"]
if (drop_isolated) {
  linked     <- unique(c(agg_show$ca, agg_show$cb))
  nodes_show <- comm_nodes[comm %in% keep_nodes & comm %in% linked]
} else {
  nodes_show <- comm_nodes[comm %in% keep_nodes]
}
message(sprintf("  Summary: %d/%d communities shown (size>=%d%s), %d between-community links (>=%d edges)",
                nrow(nodes_show), nrow(comm_nodes), min_summary_sz,
                if (drop_isolated) ", linked only" else "", nrow(agg_show), min_summary_edge))

gs <- graph_from_data_frame(d = agg_show[, .(ca, cb, n_edges, edge_type)],
                            vertices = nodes_show[, .(comm, Type, size)], directed = FALSE)
if (main_comp_only) {
  comp <- components(gs)
  gs <- induced_subgraph(gs, which(comp$membership == which.max(comp$csize)))
  message(sprintf("  main component: %d communities (%d smaller isolated component(s) of linked pairs omitted)",
                  vcount(gs), comp$no - 1))
}
set.seed(seed)
ls <- layout_with_fr(gs, niter = fr_niter)   # single connected graph -> fills the frame
sdf_n <- data.table(comm = V(gs)$name, Type = V(gs)$Type, size = V(gs)$size, x = ls[, 1], y = ls[, 2])
edt <- as.data.table(igraph::as_data_frame(gs, what = "edges"))   # from,to,n_edges,edge_type
xs2 <- setNames(sdf_n$x, sdf_n$comm); ys2 <- setNames(sdf_n$y, sdf_n$comm)
sdf_e <- edt[, .(x = xs2[from], y = ys2[from], xend = xs2[to], yend = ys2[to], n_edges, edge_type)]
sdf_e_sim  <- sdf_e[edge_type != "host"]
sdf_e_host <- sdf_e[edge_type == "host"]
lab_df <- sdf_n[size >= quantile(size, 0.75)]

message("Rendering community summary figure ...")
p_sum <- ggplot()
if (nrow(sdf_e_sim) > 0)
  p_sum <- p_sum + geom_segment(data = sdf_e_sim,
             aes(x, y, xend = xend, yend = yend, linewidth = n_edges),
             colour = "grey75", alpha = 0.55)
p_sum <- p_sum +
  geom_segment(data = sdf_e_host,
               aes(x, y, xend = xend, yend = yend, linewidth = n_edges),
               colour = "#D7263D", alpha = 0.75) +
  geom_point(data = sdf_n, aes(x, y, size = size, shape = Type, fill = Type),
             stroke = 0.3, colour = "grey20") +
  scale_shape_manual(values = c(BAC = 21, PRO = 24), name = "Community type") +
  scale_fill_manual(values = c(BAC = "#377EB8", PRO = "#FF7F00"), name = "Community type") +
  scale_size_continuous(range = c(2.5, 13), name = "Genomes in community") +
  scale_linewidth_continuous(range = c(0.2, 3), name = "Links between communities") +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  coord_equal() +
  theme_void(base_size = 13) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))

if (nrow(lab_df) > 0) {
  if (have_ggrepel) {
    p_sum <- p_sum + ggrepel::geom_text_repel(data = lab_df, aes(x, y, label = comm),
               size = 2.8, max.overlaps = 40, segment.size = 0.2, min.segment.length = 0)
  } else {
    p_sum <- p_sum + geom_text(data = lab_df, aes(x, y, label = comm), size = 2.6, vjust = -1.2)
  }
}
p_sum <- p_sum + labs(
  title = "S. suis bacteria-prophage community network (summary)",
  subtitle = sprintf("One node per community (>=%d genomes); red = prophage-host coupling between BAC & PRO communities",
                     min_summary_sz))

ggsave(summary_png, p_sum, width = 14, height = 11, dpi = 320, bg = "white")
message("  Saved ", summary_png)
message("\nDone. Outputs in ", normalizePath(out_dir))
