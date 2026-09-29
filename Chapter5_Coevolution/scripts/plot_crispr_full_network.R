#!/usr/bin/env Rscript
# Plot the full clustered CRISPR network from build_crispr_full_network.py.
# Circle = genome (ocean blue) | triangle = prophage community (sunset orange) |
# square = spacer cluster (forest green). Crimson edge = coexistence (a genome
# carrying a prophage community its own spacers target). Also exports GraphML
# for Cytoscape, where 487 genomes is comfortable to explore interactively.
suppressPackageStartupMessages({
  library(data.table); library(igraph); library(ggraph); library(tidygraph); library(ggplot2)
})
set.seed(42)
IN <- "data"
LAYOUT <- "stress"   # "stress" = even spacing (graphlayouts); try "drl" or "lgl" for more spread
nodes <- fread(file.path(IN, "crispr_full_nodes.csv"), colClasses = list(character = "name"))
edges <- fread(file.path(IN, "crispr_full_edges.csv"), colClasses = list(character = c("from","to")))

g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
COL  <- c(genome = "#0077B6", prophage_community = "#FB8500", spacer_cluster = "#228B22")
ECOL <- c("prophage in genome" = "#CED4DA", "spacer in genome" = "#74C69D",
          "protospacer hit" = "#F48C06",
          "coexistence (carries a phage it targets)" = "#9D0208")

p <- ggraph(g, layout = LAYOUT) +
  geom_edge_link(aes(edge_colour = etype, edge_width = etype, edge_alpha = etype)) +
  geom_node_point(aes(shape = type, fill = type, size = type), colour = "grey25", stroke = 0.3) +
  scale_shape_manual(values = c(genome = 21, prophage_community = 24, spacer_cluster = 22), name = "node") +
  scale_fill_manual(values = COL, name = "node") +
  scale_size_manual(values = c(genome = 2.6, prophage_community = 4.8, spacer_cluster = 1.8), name = "node") +
  scale_edge_colour_manual(values = ECOL, name = "link", breaks = names(ECOL)) +
  scale_edge_width_manual(values = c("prophage in genome"=0.3,"spacer in genome"=0.3,
                                     "protospacer hit"=0.4,"coexistence (carries a phage it targets)"=1.2),
                          guide = "none") +
  scale_edge_alpha_manual(values = c("prophage in genome"=0.35,"spacer in genome"=0.4,
                                     "protospacer hit"=0.5,"coexistence (carries a phage it targets)"=0.95),
                          guide = "none") +
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  labs(title = "Full CRISPR spacer-prophage network (spacers sequence-clustered, prophages by community)",
       subtitle = "crimson = a genome carrying a prophage community its own spacers target (coexistence/failure)  ·  orange = spacer -> prophage-community hits") +
  theme_void(base_size = 12) +
  theme(plot.subtitle = element_text(size = 8.5, colour = "grey30"), legend.position = "right")

ggsave(file.path("figures", "fig_crispr_full_network.png"), p, width = 19, height = 15, dpi = 200, bg = "white", limitsize = FALSE)
write_graph(as.igraph(g), file.path(IN, "crispr_full_network.graphml"), format = "graphml")
cat("wrote fig_crispr_full_network.png + crispr_full_network.graphml\n")
