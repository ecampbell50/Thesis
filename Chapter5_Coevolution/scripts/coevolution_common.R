# Builds the bacterial x prophage community matrix. Used by analyse_modularity.R and analyse_nestedness.R.

suppressPackageStartupMessages({ library(data.table) })

# Build the BAC-community x PRO-community binary incidence matrix.
#   res_dir   = folder holding the three Netpass CSVs
#   min_hosts = keep only prophage communities linked to >= this many DISTINCT
#               host (bacterial) genomes; smaller ones are unclustered singletons
# Returns a numeric 0/1 matrix (rows = BAC communities, cols = PRO communities).
build_baccomm_procomm <- function(res_dir = "Results", min_hosts = 5L) {
  # --- community maps: genome -> namespaced community label -----------------
  bac <- fread(file.path(res_dir, "data/genome_community_mapping_49_BAC.csv"))
  pro <- fread(file.path(res_dir, "data/genome_community_mapping_44_PRO.csv"))
  setnames(bac, c("g", "c")); setnames(pro, c("g", "c"))
  bm <- setNames(paste0("BAC_", bac$c), bac$g)   # named vector: BAC genome -> "BAC_<community>"
  pm <- setNames(paste0("PRO_", pro$c), pro$g)   # named vector: PRO genome -> "PRO_<community>"

  # --- edgetable: keep only BAC<->PRO host edges ----------------------------
  e <- fread(file.path(res_dir, "data/Bipartite_BAC_PRO_genome_edgetable.csv"))
  # logical masks for each endpoint's type (membership in a mapping = that type)
  s_bac <- e$Source %in% names(bm); t_pro <- e$Target %in% names(pm)
  s_pro <- e$Source %in% names(pm); t_bac <- e$Target %in% names(bm)
  # two orientations of a host edge: (Source=BAC,Target=PRO) and (Source=PRO,Target=BAC)
  h1 <- data.table(bac_comm = bm[e$Source[s_bac & t_pro]],
                   pro_comm = pm[e$Target[s_bac & t_pro]],
                   host     = e$Source[s_bac & t_pro])   # host = the bacterial genome
  h2 <- data.table(bac_comm = bm[e$Target[s_pro & t_bac]],
                   pro_comm = pm[e$Source[s_pro & t_bac]],
                   host     = e$Target[s_pro & t_bac])
  hl <- rbind(h1, h2)

  # --- size filter on prophage communities ----------------------------------
  # size = number of DISTINCT host genomes a prophage community is linked to
  hc   <- hl[, uniqueN(host), by = pro_comm]
  keep <- hc[V1 >= min_hosts, pro_comm]
  sub  <- hl[pro_comm %in% keep]

  # --- binary incidence matrix ----------------------------------------------
  M <- table(sub$bac_comm, sub$pro_comm)   # counts of host links per (bac_comm, pro_comm)
  M <- (M > 0) * 1                          # collapse to presence/absence (0/1)
  storage.mode(M) <- "double"
  # drop any all-zero rows/cols (shouldn't occur after filtering, but safe)
  M <- M[rowSums(M) > 0, colSums(M) > 0, drop = FALSE]

  attr(M, "n_pro_total") <- nrow(hc)        # stash diagnostics for the caller to log
  attr(M, "n_pro_kept")  <- length(keep)
  M
}

# Two-tailed empirical p-value + z-score from an observed value and null vector.
# p uses the (1 + k)/(N + 1) form: the observed network counts as one draw under
# the null, so p is never exactly 0 and the floor is 1/(N+1).
z_p <- function(obs, null) {
  mu <- mean(null); s <- sd(null)
  list(obs = obs, null_mean = mu, null_sd = s, z = (obs - mu) / s,
       p = (1 + sum(abs(null - mu) >= abs(obs - mu))) / (length(null) + 1))
}
