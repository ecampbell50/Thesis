# DESeq2, each phage vs time-matched control at every timepoint (ashr shrinkage). Used for Fig 4.3, 4.4.

if (!("DESeq2" %in% installed.packages())) BiocManager::install("DESeq2", update = FALSE)
if (!("ashr" %in% installed.packages())) BiocManager::install("ashr", update = FALSE)
if (!("ggrepel" %in% installed.packages())) install.packages("ggrepel")
if (!("RColorBrewer" %in% installed.packages())) install.packages("RColorBrewer")

library(DESeq2)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these variables before running
# ══════════════════════════════════════════════════════════════════════════════
STRAIN <- "D68"   # "C67", "D32", or "D68"
PHAGE  <- "BON"   # "BON" or "CLY" (compared against CTLR at each timepoint)

# Self-contained paths: resolve relative to this script's location in paper_package
# (…/Figure2_host_expression/{script,data,figure}). No external paths.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT  <- normalizePath(file.path(.script_dir, ".."))
DATA  <- file.path(ROOT, "paper_package", "Figure2_host_expression", "data", "DESeq2_vs_control")
OUT_CSV_DIR <- DATA                      # per-timepoint DESeq2 CSVs written back here
dir.create(OUT_CSV_DIR, showWarnings = FALSE, recursive = TRUE)

# Input counts file
COUNTS_FILE <- file.path(DATA, paste0(STRAIN, "_DESeq2data.csv"))

# Samples to exclude (failed QC replicates) — set to character(0) if none
EXCLUDE_SAMPLES <- c(0)
#EXCLUDE_SAMPLES <- c()

# Thresholds for labelling extreme genes on volcano plots
LFC_THRESHOLD  <- 1
PADJ_THRESHOLD <- 1e-10
# ══════════════════════════════════════════════════════════════════════════════


# ── Defence proteins ──────────────────────────────────────────────────────────
# Fill in the appropriate block for your strain; the others are ignored.
if (STRAIN == "C67") {

  defence_genes <- c(
    # ── Original ──────────────────────────────────────────────────────────────
    "C67_prot_00121" = "PDC-M22_S1G1",
    "C67_prot_00122" = "PDC-M22_S1G2",
    "C67_prot_00133" = "RM_Type_IV_S1G1",
    "C67_prot_00134" = "RM_Type_IV_S1G2",
    "C67_prot_00284" = "RM_Type_I_S1G1",
    "C67_prot_00285" = "RM_Type_I_S1G2",
    "C67_prot_00286" = "RM_Type_I_S1G3",
    "C67_prot_00287" = "RM_Type_I_S1G4",
    "C67_prot_00297" = "RM_Type_III_S1G1",
    "C67_prot_00298" = "RM_Type_III_S1G2",
    "C67_prot_01021" = "RM_Type_I_S2G1",
    "C67_prot_01022" = "RM_Type_I_S2G2",
    "C67_prot_01023" = "RM_Type_I_S2G3",
    "C67_prot_01024" = "RM_Type_I_S2G4",
    "C67_prot_01211" = "PDC-S30_S1G1",
    "C67_prot_01330" = "RosmerTA_S1G1",
    "C67_prot_01331" = "RosmerTA_S1G2",
    "C67_prot_01418" = "Ogmios_S1G1",
    "C67_prot_01435" = "tmn_S1G1",
    "C67_prot_01476" = "PD-T4-6_S1G1",

    # ── New from Ptolemaea ─────────────────────────────────────────────────────
    "C67_prot_00060"  = "PDC-S07_S1G1",
    "C67_prot_00603"  = "PDC-S07_S2G1",
    "C67_prot_01356"  = "PDC-S07_S3G1",
    "C67_prot_00253"  = "Gabija_S1G1",
    "C67_prot_00763"  = "CAS_Class1-IV_S1G1",
    "C67_prot_01000"  = "DRT_other_S1G1",
    "C67_prot_01545"  = "DRT_other_S2G1",
    "C67_prot_01598"  = "Theoris_S1G1",
    "C67_prot_00210"  = "DMS_other_S2G1",
    "C67_prot_00597"  = "DMS_other_S3G1",
    "C67_prot_00619"  = "DMS_other_S4G1",
    "C67_prot_01170"  = "DMS_other_S5G1",
    "C67_prot_01206"  = "DMS_other_S6G1",
    "C67_prot_01226"  = "DMS_other_S7G1",
    "C67_prot_01617"  = "DMS_other_S8G1",
    "C67_prot_01819"  = "RM_Type_I_S3G1"
  )

  prophage_genes <- c(
    "C67_prot_02009" = "hypothetical protein",
    "C67_prot_00947" = "type B 50S ribosomal protein L31",
    "C67_prot_00948" = "Putative phage integrase",
    "C67_prot_00949" = "DUF3800 domain-containing protein",
    "C67_prot_00950" = "DNA-binding phage protein",
    "C67_prot_00951" = "DNA-binding phage protein",
    "C67_prot_00952" = "Phage protein",
    "C67_prot_00953" = "Phage protein",
    "C67_prot_00954" = "Phage protein",
    "C67_prot_00955" = "Phage protein",
    "C67_prot_00956" = "hypothetical protein",
    "C67_prot_00957" = "Phage protein",
    "C67_prot_00958" = "Pathogenicity island protein",
    "C67_prot_00959" = "DNA-binding phage protein",
    "C67_prot_00960" = "Phage membrane protein",
    "C67_prot_00961" = "Phage protein",
    "C67_prot_00962" = "Phage protein",
    "C67_prot_00963" = "Phage primase",
    "C67_prot_00964" = "Site-specific recombinase",
    "C67_prot_00965" = "Rod shape determining protein"
  )

} else if (STRAIN == "D32") {
  defence_genes <- c(
    # ── CAS_Class1-IV ─────────────────────────────────────
    "D32_prot_00742" = "CAS_Class1-IV_S1G1",

    # ── DMS_other ─────────────────────────────────────────
    "D32_prot_00009" = "DMS_other_S1G1",
    "D32_prot_00085" = "DMS_other_S2G1",
    "D32_prot_00247" = "DMS_other_S3G1",
    "D32_prot_00257" = "DMS_other_S3G2",
    "D32_prot_00379" = "DMS_other_S6G1",
    "D32_prot_00401" = "DMS_other_S7G1",
    "D32_prot_00920" = "DMS_other_S8G1",
    "D32_prot_00940" = "DMS_other_S9G1",
    "D32_prot_01852" = "DMS_other_S4G1",
    "D32_prot_01911" = "DMS_other_S5G1",

    # ── DRT_other ─────────────────────────────────────────
    "D32_prot_01263" = "DRT_other_S1G1",
    "D32_prot_01431" = "DRT_other_S2G1",

    # ── Gabija ────────────────────────────────────────────
    "D32_prot_00042" = "Gabija_S1G1",

    # ── Ogmios ────────────────────────────────────────────
    "D32_prot_01386" = "Ogmios_S1G1",

    # ── PD-T4-6 ───────────────────────────────────────────
    "D32_prot_01328" = "PD-T4-6_S1G1",

    # ── PDC-M22 ───────────────────────────────────────────
    "D32_prot_00172" = "PDC-M22_S1G1",
    "D32_prot_00173" = "PDC-M22_S1G2",

    # ── PDC-S07 ───────────────────────────────────────────
    "D32_prot_00292" = "PDC-S07_S1G1",
    "D32_prot_00395" = "PDC-S07_S2G1",
    "D32_prot_01175" = "PDC-S07_S3G1",

    # ── PDC-S30 ───────────────────────────────────────────
    "D32_prot_00925" = "PDC-S30_S1G1",

    # ── RM_Type_I ─────────────────────────────────────────
    "D32_prot_00008" = "RM_Type_I_S1G1",
    "D32_prot_00010" = "RM_Type_I_S1G2",
    "D32_prot_00011" = "RM_Type_I_S1G3",
    "D32_prot_00882" = "RM_Type_I_S2G1",
    "D32_prot_00883" = "RM_Type_I_S2G2",
    "D32_prot_00884" = "RM_Type_I_S2G3",
    "D32_prot_01452" = "RM_Type_I_S3G1",
    "D32_prot_01453" = "RM_Type_I_S3G2",
    "D32_prot_01454" = "RM_Type_I_S3G3",
    "D32_prot_01455" = "RM_Type_I_S3G4",

    # ── RM_Type_III ───────────────────────────────────────
    "D32_prot_01970" = "RM_Type_III_S1G1",
    "D32_prot_01971" = "RM_Type_III_S1G2",

    # ── RM_Type_IV ────────────────────────────────────────
    "D32_prot_00160" = "RM_Type_IV_S1G1",
    "D32_prot_00161" = "RM_Type_IV_S1G2",

    # ── RosmerTA ──────────────────────────────────────────
    "D32_prot_01200" = "RosmerTA_S1G1",
    "D32_prot_01201" = "RosmerTA_S1G2",

    # ── Theoris ───────────────────────────────────────────
    "D32_prot_01867" = "Theoris_S1G1",

    # ── tmn ───────────────────────────────────────────────
    "D32_prot_01369" = "tmn_S1G1"

  )
  prophage_genes_D32 <- c(
    "D32_prot_01483" = "Rod shape determining protein",
    "D32_prot_01484" = "Site-specific recombinase",
    "D32_prot_01485" = "Phage primase",
    "D32_prot_01486" = "Age protein",
    "D32_prot_01487" = "Phage protein",
    "D32_prot_01488" = "Phage membrane protein",
    "D32_prot_01489" = "DNA-binding phage protein",
    "D32_prot_01490" = "Pathogenicity island protein",
    "D32_prot_01491" = "Phage protein",
    "D32_prot_01492" = "hypothetical protein",
    "D32_prot_01493" = "Phage protein",
    "D32_prot_01494" = "Phage protein",
    "D32_prot_01495" = "Phage protein",
    "D32_prot_01496" = "Phage protein",
    "D32_prot_01497" = "DNA-binding phage protein",
    "D32_prot_01498" = "DNA-binding phage protein",
    "D32_prot_01499" = "DUF3800 domain-containing protein",
    "D32_prot_01500" = "Putative phage integrase",
    "D32_prot_01501" = "type B 50S ribosomal protein L31"
  )

} else if (STRAIN == "D68") {

  defence_genes <- c(
    # ── Abi2DF ────────────────────────────────────────────
    "D68_prot_01256" = "Abi2DF_S1G1",
    "D68_prot_01327" = "Abi2DF_S2G1",

    # ── CAS_Class1-IV ─────────────────────────────────────
    "D68_prot_00729" = "CAS_Class1-IV_S1G1",

    # ── DMS_other ─────────────────────────────────────────
    "D68_prot_00232" = "DMS_other_S1G1",
    "D68_prot_00254" = "DMS_other_S2G1",
    "D68_prot_00291" = "DMS_other_S4G1",
    "D68_prot_00367" = "DMS_other_S5G1",
    "D68_prot_00528" = "DMS_other_S7G1",
    "D68_prot_00548" = "DMS_other_S8G1",
    "D68_prot_00919" = "DMS_other_S9G1",
    "D68_prot_00929" = "DMS_other_S9G2",
    "D68_prot_01833" = "DMS_other_S3G1",
    "D68_prot_01915" = "DMS_other_S6G1",

    # ── DRT_other ─────────────────────────────────────────
    "D68_prot_01201" = "DRT_other_S2G1",
    "D68_prot_01696" = "DRT_other_S1G1",

    # ── Gabija ────────────────────────────────────────────
    "D68_prot_00324" = "Gabija_S1G1",

    # ── Ogmios ────────────────────────────────────────────
    "D68_prot_01102" = "Ogmios_S1G1",

    # ── PD-T4-6 ───────────────────────────────────────────
    "D68_prot_01160" = "PD-T4-6_S1G1",

    # ── PDC-S07 ───────────────────────────────────────────
    "D68_prot_00237" = "PDC-S07_S1G1",
    "D68_prot_00238" = "PDC-S07_S1G2",
    "D68_prot_01046" = "PDC-S07_S2G1",

    # ── PDC-S30 ───────────────────────────────────────────
    "D68_prot_00533" = "PDC-S30_S1G1",

    # ── RM_Type_I ─────────────────────────────────────────
    "D68_prot_00289" = "RM_Type_I_S1G1",
    "D68_prot_00292" = "RM_Type_I_S1G2",
    "D68_prot_00293" = "RM_Type_I_S1G3",
    "D68_prot_00490" = "RM_Type_I_S2G1",
    "D68_prot_00491" = "RM_Type_I_S2G2",
    "D68_prot_00492" = "RM_Type_I_S2G3",
    "D68_prot_01222" = "RM_Type_I_S3G1",
    "D68_prot_01223" = "RM_Type_I_S3G2",
    "D68_prot_01224" = "RM_Type_I_S3G3",
    "D68_prot_01225" = "RM_Type_I_S3G4",

    # ── RM_Type_III ───────────────────────────────────────
    "D68_prot_01973" = "RM_Type_III_S1G1",
    "D68_prot_01974" = "RM_Type_III_S1G2",

    # ── RM_Type_IV ────────────────────────────────────────
    "D68_prot_00442" = "RM_Type_IV_S1G1",
    "D68_prot_00443" = "RM_Type_IV_S1G2",

    # ── RosmerTA ──────────────────────────────────────────
    "D68_prot_01020" = "RosmerTA_S1G1",
    "D68_prot_01021" = "RosmerTA_S1G2",

    # ── Theoris ───────────────────────────────────────────
    "D68_prot_01871" = "Theoris_S1G1",

    # ── tmn ───────────────────────────────────────────────
    "D68_prot_01119" = "tmn_S1G1"

  )
  prophage_genes_D68 <- c(
    "D68_prot_01332" = "Rod shape determining protein",
    "D68_prot_01333" = "Site-specific recombinase",
    "D68_prot_01334" = "Phage primase",
    "D68_prot_01335" = "Phage protein",
    "D68_prot_01336" = "Phage protein",
    "D68_prot_01337" = "Phage membrane protein",
    "D68_prot_01338" = "DNA-binding phage protein",
    "D68_prot_01339" = "Pathogenicity island protein",
    "D68_prot_01340" = "Phage protein",
    "D68_prot_01341" = "hypothetical protein",
    "D68_prot_01342" = "Phage protein",
    "D68_prot_01343" = "Phage protein",
    "D68_prot_01344" = "Phage protein",
    "D68_prot_01345" = "Phage protein",
    "D68_prot_01346" = "DNA-binding phage protein",
    "D68_prot_01347" = "DNA-binding phage protein",
    "D68_prot_01348" = "DUF3800 domain-containing protein",
    "D68_prot_01349" = "Putative phage integrase",
    "D68_prot_01350" = "type B 50S ribosomal protein L31"
  )

}

# Normalise prophage variable name across strains
if (STRAIN == "D32") prophage_genes <- prophage_genes_D32
if (STRAIN == "D68") prophage_genes <- prophage_genes_D68


# ── Data loading ──────────────────────────────────────────────────────────────
counts <- read.csv(COUNTS_FILE, row.names = 1)
colnames(counts) <- sub("^X", "", colnames(counts))

# Select columns for the chosen PHAGE and the CTLR at all timepoints
counts_sub <- counts[, grepl(paste0(PHAGE, "|CTLR"), colnames(counts))]
counts_sub <- counts_sub[, !colnames(counts_sub) %in% EXCLUDE_SAMPLES]

# ── Parse sample metadata from column names ──────────────────────────────────
# Column format: {time}_{strain}_{treatment}_{replicate}
col_parts <- strsplit(colnames(counts_sub), "_")
colData_sub <- data.frame(
  row.names = colnames(counts_sub),
  time      = factor(sapply(col_parts, `[`, 1)),
  treatment = factor(sapply(col_parts, `[`, 3), levels = c("CTLR", PHAGE))
)

# Combined group factor for pairwise contrasts at each timepoint
colData_sub$group <- factor(paste0(colData_sub$time, "_", colData_sub$treatment))

# Derive timepoints (sorted numerically)
all_times     <- sort(unique(as.numeric(as.character(colData_sub$time))))
all_times_chr <- as.character(all_times)


# ── DESeq2 ────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData   = colData_sub,
  design    = ~ group
)

keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
dds  <- DESeq(dds)
resultsNames(dds)


# ── Results for each timepoint: PHAGE vs CTLR ────────────────────────────────
results_list <- list()
for (tp in all_times_chr) {
  res <- lfcShrink(dds,
                   contrast = c("group",
                                paste0(tp, "_", PHAGE),
                                paste0(tp, "_CTLR")),
                   type = "ashr")
  results_list[[tp]] <- res
  write.csv(as.data.frame(res),
            file = file.path(OUT_CSV_DIR, paste0("DESeq2_", STRAIN, "_", PHAGE, "_vs_CTLR_T", tp, "min.csv")))
  cat("\n--- T=", tp, "min: ", PHAGE, " vs CTLR ---\n")
  print(summary(res, alpha = 0.05))
}


# ── VST ───────────────────────────────────────────────────────────────────────
vsd <- vst(dds, blind = TRUE)


# ── PCA ───────────────────────────────────────────────────────────────────────
time_colours <- c("0" = "green", "2" = "yellow", "10" = "orange",
                  "20" = "red",  "30" = "purple", "50" = "blue")
treatment_shapes <- c("CTLR" = 1, "BON" = 16, "CLY" = 17)

pca_data    <- plotPCA(vsd, intgroup = c("time", "treatment"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                  colour = time, shape = treatment,
                                  label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf, box.padding = 0.5) +
  scale_colour_manual(values = time_colours) +
  scale_shape_manual(values = treatment_shapes) +
  labs(x      = paste0("PC1: ", percent_var[1], "% variance"),
       y      = paste0("PC2: ", percent_var[2], "% variance"),
       colour = "Timepoint (min)",
       shape  = "Treatment",
       title  = paste0(STRAIN, " — ", PHAGE, " vs CTLR PCA")) +
  theme_bw()

# png(paste0("PCA_", STRAIN, "_", PHAGE, "_vs_CTLR.png"), width = 800, height = 500)
# print(pca_plot)
# dev.off()
print(pca_plot)


# ── Sample distance heatmap ───────────────────────────────────────────────────
sampleDists      <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
heatmap_colours  <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# png(paste0("heatmap_", STRAIN, "_", PHAGE, "_vs_CTLR.png"), width = 800, height = 700)
print(pheatmap(sampleDistMatrix,
               clustering_distance_rows = sampleDists,
               clustering_distance_cols = sampleDists,
               col  = heatmap_colours,
               main = paste0(STRAIN, " — ", PHAGE, " vs CTLR Sample Distances")))
# dev.off()


# ── Faceted volcano (all timepoints in one figure) ───────────────────────────

# Combine all timepoint results into a single data frame
all_results <- bind_rows(
  lapply(all_times_chr, function(tp) {
    as.data.frame(results_list[[tp]]) %>%
      rownames_to_column("gene") %>%
      filter(!is.na(padj)) %>%
      mutate(timepoint = tp)
  })
) %>%
  mutate(
    timepoint     = factor(timepoint, levels = all_times_chr),
    sig           = padj < 0.05 & abs(log2FoldChange) > 1,
    is_defence    = gene %in% names(defence_genes),
    defence_label = ifelse(is_defence, defence_genes[gene], NA_character_),
    is_prophage    = gene %in% names(prophage_genes),
    prophage_label = ifelse(is_prophage, prophage_genes[gene], NA_character_),
    is_extreme    = sig & !is_defence & !is_prophage &
                    (abs(log2FoldChange) > LFC_THRESHOLD | padj < PADJ_THRESHOLD),
    category      = case_when(
      is_defence  ~ "defence",
      is_prophage ~ "prophage",
      sig         ~ "significant",
      TRUE        ~ "not_significant"
    )
  )

# Per-facet direction annotations
annot_right <- data.frame(
  timepoint = factor(all_times_chr, levels = all_times_chr),
  label     = paste0("Higher in ", PHAGE)
)
annot_left <- data.frame(
  timepoint = factor(all_times_chr, levels = all_times_chr),
  label     = "Higher in CTLR"
)

volcano_facet <- ggplot(all_results, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = filter(all_results, category == "not_significant"),
             colour = "grey70", alpha = 0.5, size = 1) +
  geom_point(data = filter(all_results, category == "significant"),
             colour = "red", alpha = 0.7, size = 1) +
  geom_point(data = filter(all_results, is_defence),
             colour = "dodgerblue", size = 2.5, shape = 18) +
  geom_point(data = filter(all_results, is_prophage),
             colour = "forestgreen", size = 2.5, shape = 18) +
  geom_text_repel(data = filter(all_results, is_extreme),
                  aes(label = gene),
                  colour = "black", size = 2,
                  max.overlaps = 15, box.padding = 0.3,
                  fontface = "italic") +
  geom_text_repel(data = filter(all_results, is_defence & sig),
                  aes(label = defence_label),
                  colour = "dodgerblue", size = 1.5,
                  max.overlaps = Inf, box.padding = 0.5,
                  fontface = "bold") +
  geom_text_repel(data = filter(all_results, is_prophage & sig),
                  aes(label = prophage_label),
                  colour = "forestgreen", size = 1.5,
                  max.overlaps = Inf, box.padding = 0.5,
                  fontface = "bold") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  geom_text(data = annot_left,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.5, inherit.aes = FALSE,
            colour = "grey40", size = 2.5, fontface = "italic") +
  geom_text(data = annot_right,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.5, inherit.aes = FALSE,
            colour = "grey40", size = 2.5, fontface = "italic") +
  facet_wrap(~ timepoint, nrow = 2,
             labeller = labeller(
               timepoint = function(x) paste0("T=", x, "min: ", PHAGE, " vs CTLR"))) +
  labs(title    = paste0(STRAIN, " — ", PHAGE, " vs CTLR at Each Timepoint"),
       subtitle = paste0("red = sig. DE (padj<0.05, |LFC|>1)  |  ",
                         "blue = defence genes (n=", length(defence_genes), ")  |  ",
                         "green = prophage genes (n=", length(prophage_genes), ")  |  ",
                         "labelled = |LFC|>", LFC_THRESHOLD, " or padj<", PADJ_THRESHOLD),
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  theme_bw(base_size = 12) +
  theme(plot.subtitle = element_text(size = 9, colour = "grey40"),
        aspect.ratio = 1)   # makes each facet panel square

# ── View at high quality (opens in Preview/default image viewer) ──────────────
tmp <- tempfile(fileext = ".png")
ggsave(tmp, volcano_facet, width = 20, height = 14, dpi = 300)
browseURL(tmp)

# ── Save permanently (uncomment when ready) ───────────────────────────────────
# ggsave(paste0("Volcano_faceted_", STRAIN, "_", PHAGE, "_vs_CTLR.png"),
#        volcano_facet, width = 20, height = 14, dpi = 300)
