#!/usr/bin/env Rscript
# Nestedness (NODF) of the community infection network vs r2dtable nulls (Fig 5.15).

suppressPackageStartupMessages({
  library(bipartite); library(data.table); library(parallel)
})
source("coevolution_common.R")   # build_baccomm_procomm(), z_p()

set.seed(42)                      # reproducible nulls
RES_DIR    <- "Results"
OUT_DIR    <- file.path(RES_DIR, "coevolution")
THRESHOLDS <- c(3L, 5L, 10L)      # minimum distinct hosts per prophage community
N_NULL     <- 999                 # nulls (cheap for NODF) -> p-floor 1/1000 = 0.001
N_CORES    <- max(1L, detectCores() - 1L)
OUTCSV     <- file.path(OUT_DIR, "data/network_nestedness_results.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (file.exists(OUTCSV)) file.remove(OUTCSV)

# NODF of a web. bipartite::nested(method = "NODF2") returns the standard
# Almeida-Neto NODF on a 0-100 scale (higher = more nested).
nodf <- function(web) as.numeric(nested(web, method = "NODF2"))

analyse_one <- function(thr) {
  web <- build_baccomm_procomm(RES_DIR, min_hosts = thr)
  cat(sprintf("\n=== min_hosts >= %d : %d BAC-comm x %d PRO-comm  (%d/%d prophage communities kept) ===\n",
              thr, nrow(web), ncol(web),
              attr(web, "n_pro_kept"), attr(web, "n_pro_total")))

  obs <- nodf(web)
  cat(sprintf("  observed NODF = %.3f ; running %d nulls on %d cores ...\n",
              obs, N_NULL, N_CORES))

  # r2dtable nulls: random matrices with the SAME row & column sums as the web.
  nulls <- nullmodel(web, N = N_NULL, method = "r2dtable")
  nv    <- unlist(mclapply(nulls, nodf, mc.cores = N_CORES))   # NODF of each null, in parallel

  r <- z_p(obs, nv)
  cat(sprintf("  NODF = %.3f  null = %.3f +/- %.3f  z = %+.2f  p = %.4f\n",
              r$obs, r$null_mean, r$null_sd, r$z, r$p))

  # --- packed-matrix figure (standard nestedness visualisation) -------------
  # visweb with rows/cols sorted by marginal totals shows the nested "staircase":
  # a fully nested matrix packs into a triangle in the top-left.
  png(file.path(OUT_DIR, sprintf("fig_nested_min%d.png", thr)), 1600, 1200, res = 150)
  tryCatch(
    visweb(web, type = "nested", labsize = 0.4,
           xlab = "Prophage communities", ylab = "Bacterial communities"),
    error = function(e) plot.new())
  title(main = sprintf("Nestedness (NODF = %.1f, z = %+.1f, p = %.3f)  |  min hosts >= %d",
                       r$obs, r$z, r$p, thr))
  dev.off()

  data.table(network = "BACcomm_x_prophageComm", min_hosts = thr,
             n_row = nrow(web), n_col = ncol(web), metric = "nestedness_NODF",
             observed = r$obs, null_mean = r$null_mean, null_sd = r$null_sd,
             z = r$z, p = r$p, n_null = N_NULL)
}

# run all thresholds, append results incrementally
for (thr in THRESHOLDS) {
  out <- analyse_one(thr)
  fwrite(out, OUTCSV, append = file.exists(OUTCSV))   # header only on first write
}
cat("\nDONE nestedness — see", OUTCSV, "\n")
