#!/usr/bin/env Rscript
# ============================================================================
# analyse_modularity.R  —  MODULARITY (Beckett's Q) of the bacterial-community
# x prophage-community infection network, computed straight from Netapss output.
#
# Q (Beckett 2016, DIRTLPAwb+) measures whether the network partitions into
# compartments of communities that interact more among themselves than expected.
# The optimiser is a HEURISTIC, so each computation takes the BEST OF 10 random
# restarts. CRUCIALLY the same best-of-10 is applied to the observed network AND
# to every null, so neither side is favoured (avoids inflating significance).
#
# Matrix construction is shared with analyse_nestedness.R via coevolution_common.R
# (identical logic). Three host-count thresholds are run as a sensitivity test:
# >=5 distinct hosts is the primary network; >=3 and >=10 bracket it.
#
# Significance: observed Q vs an r2dtable null (random matrices preserving the
# row and column totals), as a z-score and empirical p-value.
#
# SLOW: Beckett on the >=3 matrix is ~260 s/run; with best-of-10 x 99 nulls this
# is an overnight job. Results are written INCREMENTALLY (one threshold at a
# time) so partial output lands early.
#
# Out (Results/coevolution/):
#   network_modularity_results.csv      one row per threshold (appended)
#   modules_min<thr>.csv                per-node module membership
#   fig_modules_min<thr>.png            module-web visualisation
# ============================================================================
suppressPackageStartupMessages({
  library(bipartite); library(data.table); library(parallel)
})
source("coevolution_common.R")   # build_baccomm_procomm(), z_p()

set.seed(42)                      # reproducible restarts + nulls
RES_DIR    <- "Results"
OUT_DIR    <- file.path(RES_DIR, "coevolution")
THRESHOLDS <- c(3L, 5L, 10L)      # minimum distinct hosts per prophage community
N_NULL     <- 99                  # nulls -> p-floor 1/100 = 0.01 (z-score carries the weight)
N_RESTART  <- 10                  # best-of-N restarts, applied to BOTH observed AND nulls
N_CORES    <- max(1L, detectCores() - 1L)
OUTCSV     <- file.path(OUT_DIR, "data/network_modularity_results.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (file.exists(OUTCSV)) file.remove(OUTCSV)

# Best-of-N modularity: run Beckett N times, return the highest-Q module object.
# Identical procedure is used for the observed web and every null (fairness).
best_modules <- function(web, n = N_RESTART) {
  best <- -Inf; bm <- NULL
  for (i in seq_len(n)) {
    mo <- computeModules(web, method = "Beckett")
    if (mo@likelihood > best) { best <- mo@likelihood; bm <- mo }
  }
  bm
}
best_Q <- function(web, n = N_RESTART) best_modules(web, n)@likelihood

analyse_one <- function(thr) {
  web <- build_baccomm_procomm(RES_DIR, min_hosts = thr)
  cat(sprintf("\n=== min_hosts >= %d : %d BAC-comm x %d PRO-comm  (%d/%d prophage communities kept) ===\n",
              thr, nrow(web), ncol(web),
              attr(web, "n_pro_kept"), attr(web, "n_pro_total")))

  # observed: best of N_RESTART
  mo  <- best_modules(web); obs <- mo@likelihood
  cat(sprintf("  observed Q = %.3f (best of %d) ; running %d nulls x %d restarts on %d cores ...\n",
              obs, N_RESTART, N_NULL, N_RESTART, N_CORES))

  # nulls: r2dtable (preserve row & col sums), each scored best-of-N, in parallel
  nulls <- nullmodel(web, N = N_NULL, method = "r2dtable")
  nv    <- unlist(mclapply(nulls, best_Q, mc.cores = N_CORES))

  r    <- z_p(obs, nv)
  cons <- module2constraints(mo)   # per-node module id (rows then cols)
  cat(sprintf("  Q = %.3f  null = %.3f +/- %.3f  z = %+.2f  p = %.4f  modules = %d\n",
              r$obs, r$null_mean, r$null_sd, r$z, r$p, length(unique(cons))))

  # per-node module membership
  fwrite(data.table(node    = c(rownames(web), colnames(web)),
                    side    = c(rep("BAC_comm", nrow(web)), rep("PRO_comm", ncol(web))),
                    module  = cons, min_hosts = thr),
         file.path(OUT_DIR, sprintf("modules_min%d.csv", thr)))

  # module-web figure
  png(file.path(OUT_DIR, sprintf("fig_modules_min%d.png", thr)), 1600, 1200, res = 150)
  tryCatch(plotModuleWeb(mo, labsize = 0.4), error = function(e) plot.new())
  title(main = sprintf("Modularity (Q = %.3f, z = %+.1f, p = %.3f, %d modules)  |  min hosts >= %d",
                       r$obs, r$z, r$p, length(unique(cons)), thr))
  dev.off()

  out <- data.table(network = "BACcomm_x_prophageComm", min_hosts = thr,
                    n_row = nrow(web), n_col = ncol(web), metric = "modularity_Q",
                    observed = r$obs, null_mean = r$null_mean, null_sd = r$null_sd,
                    z = r$z, p = r$p, n_modules = length(unique(cons)),
                    n_null = N_NULL, n_restart = N_RESTART)
  fwrite(out, OUTCSV, append = file.exists(OUTCSV))   # incremental: header only on first write
  out
}

# Order: smallest/fastest matrix LAST is not ideal for "partial output early",
# but biology wants >=5 (primary) early. >=10 is fastest, >=5 next, >=3 slowest.
# Run fastest -> slowest so partial results land as soon as possible.
for (thr in c(10L, 5L, 3L)) analyse_one(thr)
cat("\nDONE modularity — see", OUTCSV, "\n")
