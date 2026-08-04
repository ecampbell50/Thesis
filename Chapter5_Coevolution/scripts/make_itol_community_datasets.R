#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# make_itol_community_datasets.R
#
# Builds two iToL datasets for the bacterial phylogeny that use the SAME
# top-15-community selection and the SAME colours as the bipartite infection
# network figure (plot_bipartite_network_ANNOTATED.R):
#
#   1. Dataset_Colourstrip_BACcom_top15.txt  (DATASET_COLORSTRIP)
#        A coloured strip on each bacterial tip, coloured by its BAC community
#        (only the 15 largest BAC communities; others left blank).
#
#   2. Dataset_Symbol_PROcom_top15.txt       (DATASET_SYMBOL)
#        A filled triangle placed at each host tip for every top-15 PRO
#        community that host carries (echoing the prophage=triangle glyph in
#        the bipartite figure). A host with prophages from several top-15
#        communities gets several triangles.
#
# The top-15 selection and colour->community mapping are computed with the
# IDENTICAL code path as the bipartite script, so the figures stay consistent.
# ---------------------------------------------------------------------------

suppressMessages(library(data.table))

# ---- paths (mirror plot_bipartite_network_ANNOTATED.R) --------------------
in_dir    <- "Results"
node_file <- file.path(in_dir, "data/GenomeType_nodetable.csv")
bac_map   <- file.path(in_dir, "data/genome_community_mapping_49_BAC.csv")
pro_map   <- file.path(in_dir, "data/genome_community_mapping_44_PRO.csv")
tree_file <- file.path("iToL", "data/SsuisPhylo_FastTree_11Sep24.treefile")
out_dir   <- "iToL"

top_k_each <- 15

# ---- colour palettes (copied verbatim from the bipartite script) ----------
bac_pal <- c("#DBFFFF","#E9DCFF","#FFD2FB","#BEFFDA","#F1FFDB","#EFFFC5","#FFF4D1",
             "#D2ECFF","#FFDBDB","#FCFDAF","#FFE2D1","#FFBBDA","#E8CFF8","#90F1EF","#FFEF9F")
pro_pal <- c("#FF0A54","#D7FF00","#007BFF","#FFD500","#00BFFF","#FF008A","#A600FF",
             "#FF9100","#00FF4B","#E600FF","#A7FF00","#FF5768","#00F8FF","#5A00FF","#FF6F3C")

# ---- read data exactly as the bipartite script does -----------------------
nodes <- fread(node_file); setnames(nodes, c("Genome", "Type"))
bac <- fread(bac_map); setnames(bac, c("Genome", "Community")); bac[, comm := paste0("BAC_com_", Community)]
pro <- fread(pro_map); setnames(pro, c("Genome", "Community")); pro[, comm := paste0("PRO_com_", Community)]

comm_map  <- rbindlist(list(bac[, .(Genome, comm)], pro[, .(Genome, comm)]))
node_attr <- merge(nodes, comm_map, by = "Genome", all.x = TRUE)
node_attr[is.na(comm), comm := "unassigned"]

# ---- top-15 selection + colour maps (identical to bipartite script) -------
comm_sizes <- node_attr[comm != "unassigned", .N, by = comm][order(-N)]
bac_top <- head(comm_sizes[grepl("^BAC", comm), comm], top_k_each)
pro_top <- head(comm_sizes[grepl("^PRO", comm), comm], top_k_each)

bac_col <- setNames(bac_pal[seq_along(bac_top)], bac_top)   # comm -> hex
pro_col <- setNames(pro_pal[seq_along(pro_top)], pro_top)

# nice legend labels: "BAC_com_3" -> "BAC community 3"
pretty_lab <- function(x) sub("_com_", " community ", x)

# ---- tree tips (bacterial genome IDs, no .fna) ----------------------------
tree  <- ape::read.tree(tree_file)
tips  <- tree$tip.label
strip_fna <- function(x) sub("\\.fna$", "", x)

# =====================================================================
# 1. BAC community colour strip
# =====================================================================
bac[, id := strip_fna(Genome)]
bac_dat <- bac[comm %in% bac_top & id %in% tips]
bac_dat[, col := bac_col[comm]]
bac_dat[, lab := pretty_lab(comm)]
# keep top-15 legend order (largest first)
bac_dat[, comm := factor(comm, levels = bac_top)]
setorder(bac_dat, comm)

f1 <- file.path(out_dir, "Dataset_Colourstrip_BACcom_top15.txt")
con <- file(f1, "w")
writeLines(c(
  "DATASET_COLORSTRIP",
  "SEPARATOR SPACE",
  "DATASET_LABEL BAC_community_top15",
  "COLOR #1f78b4",
  "COLOR_BRANCHES 0",
  "STRIP_WIDTH 40",
  "MARGIN 2",
  paste("LEGEND_TITLE", "Bacterial_community_(top15)"),
  paste("LEGEND_SHAPES", paste(rep(1, length(bac_top)), collapse = " ")),
  paste("LEGEND_COLORS", paste(bac_col[bac_top], collapse = " ")),
  paste("LEGEND_LABELS", paste(gsub(" ", "_", pretty_lab(bac_top)), collapse = " ")),
  "DATA"
), con)
writeLines(sprintf("%s %s %s", bac_dat$id, bac_dat$col, gsub(" ", "_", bac_dat$lab)), con)
close(con)

# =====================================================================
# 2. PRO community external-shape grid (aligned columns beside the tips)
# =====================================================================
# host genome = substring of the prophage ID before the first "_".
# NOTE: we keep every prophage row (no dedup) and COUNT how many prophages of
# each top-15 community each host carries, so all instances are represented.
pro[, host := sub("_.*$", "", Genome)]
pro_dat <- pro[comm %in% pro_top & host %in% tips]

# host x community count matrix (rows = host tips, one column per top-15 comm)
counts <- dcast(pro_dat, host ~ comm, fun.aggregate = length, value.var = "Genome")
# make sure every top-15 community is present as a column, in largest-first order
for (cc in pro_top) if (!cc %in% names(counts)) counts[, (cc) := 0L]
setcolorder(counts, c("host", pro_top))

# External shape field definitions (one field per top-15 prophage community).
#   FIELD_SHAPES codes: 1 rectangle, 2 circle, 3 star, 4 right triangle,
#   5 left triangle, 6 checkmark. Use 4 (triangle) to echo the prophage glyph.
f2 <- file.path(out_dir, "Dataset_ExternalShape_PROcom_top15.txt")
con <- file(f2, "w")
writeLines(c(
  "DATASET_EXTERNALSHAPE",
  "SEPARATOR SPACE",
  "DATASET_LABEL PRO_community_top15",
  "COLOR #ff7f00",
  paste("FIELD_COLORS", paste(pro_col[pro_top], collapse = " ")),
  paste("FIELD_LABELS", paste(gsub(" ", "_", pretty_lab(pro_top)), collapse = " ")),
  paste("FIELD_SHAPES", paste(rep(4, length(pro_top)), collapse = " ")),
  "SHAPE_SPACING 2",
  "SHAPE_TYPE 2",      # 2 = all shapes uniform size (value>0 = shown); 1 = size scaled by count
  "COLOR_FILL 1",
  "HORIZONTAL_GRID 0",
  "VERTICAL_GRID 0",
  "SHOW_VALUES 0",
  paste("LEGEND_TITLE", "Prophage_community_(top15)"),
  paste("LEGEND_SHAPES", paste(rep(4, length(pro_top)), collapse = " ")),
  paste("LEGEND_COLORS", paste(pro_col[pro_top], collapse = " ")),
  paste("LEGEND_LABELS", paste(gsub(" ", "_", pretty_lab(pro_top)), collapse = " ")),
  "DATA"
), con)
# Each data row: host  <count for comm1>  <count for comm2> ... <count comm15>.
# A cell value of 0 draws nothing; larger counts draw a larger triangle, so a
# host carrying several prophages of one community shows a bigger symbol.
mat <- as.matrix(counts[, ..pro_top])
writeLines(paste(counts$host, apply(mat, 1, paste, collapse = " ")), con)
close(con)

# ---- report ---------------------------------------------------------------
message("Top-15 BAC communities (largest first) + colour:")
for (i in seq_along(bac_top)) message(sprintf("  %-12s %s", pretty_lab(bac_top[i]), bac_col[bac_top[i]]))
message("\nTop-15 PRO communities (largest first) + colour:")
for (i in seq_along(pro_top)) message(sprintf("  %-12s %s", pretty_lab(pro_top[i]), pro_col[pro_top[i]]))

message(sprintf("\nBAC strip: %d tips coloured (of %d tree tips)", nrow(bac_dat), length(tips)))
message(sprintf("PRO grid: %d prophage instances (top-15 comms) across %d host tips",
                nrow(pro_dat), nrow(counts)))
message(sprintf("  max prophages of one community in a single host: %d", max(mat)))
message(sprintf("\nWrote:\n  %s\n  %s", f1, f2))
