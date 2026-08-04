#!/usr/bin/env Rscript
# ============================================================================
# community_phylo_distance.R  —  Q5: phylogenetic homogeneity per BAC community.
# Mean pairwise patristic (cophenetic) distance between member genomes on the
# FastTree tree. Low = tight/clonal community; high = phylogenetically diverse.
# Out: Results/community_questions/q5_community_phylo_distance.csv
# ============================================================================
suppressPackageStartupMessages({ library(ape); library(data.table) })

TREE <- "data/SsuisPhylo_FastTree_11Sep24.treefile"
BAC  <- "data/genome_community_mapping_49_BAC.csv"
OUT  <- "data/q5_community_phylo_distance.csv"

tr  <- read.tree(TREE)
# normalise tip + mapping IDs (strip common extensions) so they join
norm <- function(x) sub("\\.(fna|fasta|ref)$", "", x)
tip_norm <- norm(tr$tip.label)

bac <- fread(BAC, colClasses = list(character = 1))   # colClasses avoids the float-ID trap
setnames(bac, c("g", "c")); bac[, g := norm(g)]; bac[, bac_comm := paste0("BAC_com_", c)]
bac <- bac[g %in% tip_norm]
cat(sprintf("tips: %d | BAC genomes on tree: %d / %d\n",
            length(tr$tip.label), nrow(bac), length(tip_norm)))

# full patristic distance matrix (2124 x 2124 is fine), index by normalised tip
D <- cophenetic(tr)
rownames(D) <- tip_norm; colnames(D) <- tip_norm

res <- bac[, {
  tips <- g
  if (length(tips) < 2) {
    .(n_genomes = length(tips), mean_pair_dist = NA_real_, sd_pair_dist = NA_real_)
  } else {
    sub <- D[tips, tips]
    vals <- sub[upper.tri(sub)]            # unique pairwise distances
    .(n_genomes = length(tips), mean_pair_dist = mean(vals), sd_pair_dist = sd(vals))
  }
}, by = bac_comm]

setorder(res, mean_pair_dist, na.last = TRUE)
fwrite(res, OUT)
cat("wrote", OUT, "\n")
cat(sprintf("communities with >=2 genomes on tree: %d\n", sum(!is.na(res$mean_pair_dist))))
print(head(res[!is.na(mean_pair_dist)], 5))
