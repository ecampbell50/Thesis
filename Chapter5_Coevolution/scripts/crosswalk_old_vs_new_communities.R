#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# crosswalk_old_vs_new_communities.R
#
# For each of the new top-15 BAC/PRO communities (from the noBAC2/noPRO4
# Netpass re-runs), shows exactly which OLD community(ies) its member
# genomes/prophages came from, and in what proportions. This replaces
# eyeballing colour-matched rings on the tree (which is unreliable, since the
# old and new legends assign the SAME colour to DIFFERENT community numbers)
# with a direct genome-level crosswalk.
# ---------------------------------------------------------------------------

suppressMessages(library(data.table))

in_dir      <- "."
test_dir    <- "."
bac_map_old <- file.path(in_dir, "data/genome_community_mapping_49_BAC.csv")
pro_map_old <- file.path(in_dir, "data/genome_community_mapping_44_PRO.csv")
bac_map_new <- file.path(test_dir, "data/genome_community_mapping_48_BAC.csv")
pro_map_new <- file.path(test_dir, "data/genome_community_mapping_44_PRO.csv")

top_k_each <- 15

crosswalk <- function(old_file, new_file, prefix) {
  old <- fread(old_file); setnames(old, c("Genome", "Community"))
  old[, old_comm := paste0(prefix, "_com_", Community)]
  new <- fread(new_file); setnames(new, c("Genome", "Community"))
  new[, new_comm := paste0(prefix, "_com_", Community)]

  new_sizes <- new[, .N, by = new_comm][order(-N)]
  new_top <- head(new_sizes$new_comm, top_k_each)

  merged <- merge(new[, .(Genome, new_comm)], old[, .(Genome, old_comm)], by = "Genome", all.x = TRUE)
  n_unmatched <- merged[is.na(old_comm), .N]
  if (n_unmatched > 0) message(sprintf("[%s] WARNING: %d retained genomes had no old-community match (unexpected)", prefix, n_unmatched))

  message(sprintf("\n======================================================================"))
  message(sprintf("%s: composition of each NEW top-%d community, by OLD community", prefix, top_k_each))
  message(sprintf("======================================================================"))

  out <- list()
  for (nc in new_top) {
    sub <- merged[new_comm == nc]
    tab <- sub[, .N, by = old_comm][order(-N)]
    tab[, pct := round(100 * N / sum(N), 1)]
    total <- sum(tab$N)
    top_contributors <- head(tab, 5)
    contrib_str <- paste(sprintf("%s (%d, %.1f%%)", top_contributors$old_comm, top_contributors$N, top_contributors$pct),
                          collapse = "; ")
    n_sources <- nrow(tab)
    message(sprintf("  NEW %-14s n=%-4d  <- from %d old communities: %s%s",
                     nc, total, n_sources, contrib_str,
                     if (n_sources > 5) sprintf(" ... (+%d more)", n_sources - 5) else ""))
    tab[, new_comm := nc]
    out[[nc]] <- tab
  }

  result <- rbindlist(out)
  setcolorder(result, c("new_comm", "old_comm", "N", "pct"))
  result
}

bac_result <- crosswalk(bac_map_old, bac_map_new, "BAC")
pro_result <- crosswalk(pro_map_old, pro_map_new, "PRO")

fwrite(bac_result, file.path(test_dir, "data/crosswalk_BAC_new_vs_old.csv"))
fwrite(pro_result, file.path(test_dir, "data/crosswalk_PRO_new_vs_old.csv"))

message(sprintf("\nWrote:\n  %s\n  %s",
                 file.path(test_dir, "data/crosswalk_BAC_new_vs_old.csv"),
                 file.path(test_dir, "data/crosswalk_PRO_new_vs_old.csv")))
