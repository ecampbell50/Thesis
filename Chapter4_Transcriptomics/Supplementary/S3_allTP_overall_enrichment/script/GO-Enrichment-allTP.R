suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(GO.db)
  library(AnnotationDbi)
})

# Input: Rscript GO-Enrichment-allTP.R STRAIN PHAGE
# Unions DEGs across all post-infection timepoints (T2,10,20,30,50) and runs
# a single GO enrichment — "which pathways are broadly active during infection?"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  STRAIN <- args[1]; PHAGE <- args[2]
} else {
  STRAIN <- "C67"; PHAGE <- "CLY"
}

TIMES     <- c(2, 10, 20, 30, 50)
# Self-contained: resolve relative to this script's location in paper_package
# (…/S3_allTP_overall_enrichment/{script,data}). Inputs live in data/_inputs/.
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
  error = function(e) getwd())
ROOT      <- normalizePath(file.path(.script_dir, ".."))
DATA      <- file.path(ROOT, "data")
func_file <- file.path(DATA, "_inputs", paste0("functional_annotation_", STRAIN, "_with_names.csv"))

out_dir <- file.path(DATA, STRAIN, PHAGE, paste0(STRAIN, "_", PHAGE, "_allTP_GO_results"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Filters
padj_thr <- 0.05
lfc_thr  <- 1

# Helpers
clean_term_string <- function(x) {
  x <- str_trim(x)
  x[x == ""] <- NA_character_
  x
}

parse_ratio <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || val == "") return(NA_real_)
    parts <- strsplit(val, "/")[[1]]
    if (length(parts) != 2) return(NA_real_)
    num <- as.numeric(parts[1]); den <- as.numeric(parts[2])
    if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
    num / den
  })
}

safe_min_positive <- function(x) {
  x <- x[is.finite(x) & !is.na(x) & x > 0]
  if (length(x) == 0) return(1e-10)
  min(x)
}

# Load annotation
if (!file.exists(func_file)) stop("Functional annotation file not found: ", func_file)
anno <- read_csv(func_file, show_col_types = FALSE)

required_anno_cols <- c("GeneID", "GOs", "GO_names")
missing_anno_cols  <- setdiff(required_anno_cols, colnames(anno))
if (length(missing_anno_cols) > 0) {
  stop("Annotation file missing columns: ", paste(missing_anno_cols, collapse = ", "))
}

anno <- anno %>% filter(!str_starts(GeneID, "PEN"), !str_starts(GeneID, "ZUD"))

# Read all timepoint DESeq2 files and collect genes
all_genes_tested <- character(0)
all_sig_genes    <- character(0)

for (tp in TIMES) {
  f <- file.path(DATA, "_inputs", "DESeq2_vs_control", paste0("DESeq2_", STRAIN, "_", PHAGE, "_vs_CTLR_T", tp, "min.csv"))
  if (!file.exists(f)) { message("Missing file, skipping: ", f); next }

  df <- read_csv(f, show_col_types = FALSE)
  if ("...1" %in% colnames(df)) colnames(df)[colnames(df) == "...1"] <- "gene"
  if (!"gene" %in% colnames(df)) next

  df <- df %>% filter(!str_starts(gene, "PEN"), !str_starts(gene, "ZUD"))

  all_genes_tested <- union(all_genes_tested, unique(df$gene))

  if (all(c("padj", "log2FoldChange") %in% colnames(df))) {
    sig <- df %>%
      filter(!is.na(padj), padj < padj_thr,
             !is.na(log2FoldChange), abs(log2FoldChange) >= lfc_thr) %>%
      pull(gene)
    all_sig_genes <- union(all_sig_genes, sig)
  }
}

genes_bg  <- intersect(all_genes_tested, unique(anno$GeneID))
genes_sig <- intersect(all_sig_genes,    genes_bg)

if (length(genes_bg)  == 0) stop("No background genes after intersection with annotation.")
if (length(genes_sig) == 0) stop("No significant genes across any timepoint.")

message(STRAIN, " ", PHAGE, " — background: ", length(genes_bg),
        " genes | significant (union): ", length(genes_sig), " genes")

# Build TERM2GENE and TERM2NAME
term2gene_go <- anno %>%
  dplyr::select(GeneID, GOs) %>%
  filter(!is.na(GOs), GOs != "") %>%
  separate_rows(GOs, sep = ";") %>%
  mutate(GOs = clean_term_string(GOs)) %>%
  filter(!is.na(GOs)) %>%
  transmute(term = GOs, gene = GeneID) %>%
  distinct()

term2name_go <- anno %>%
  dplyr::select(GOs, GO_names) %>%
  filter(!is.na(GOs), GOs != "", !is.na(GO_names), GO_names != "") %>%
  separate_rows(GOs, sep = ";") %>%
  separate_rows(GO_names, sep = "; ") %>%
  mutate(GOs = clean_term_string(GOs), GO_names = clean_term_string(GO_names)) %>%
  filter(!is.na(GOs), !is.na(GO_names)) %>%
  transmute(term = GOs, name = GO_names) %>%
  distinct(term, .keep_all = TRUE)

go_ontology <- AnnotationDbi::select(
  GO.db, keys = unique(term2gene_go$term), columns = "ONTOLOGY", keytype = "GOID"
) %>%
  as_tibble() %>%
  filter(!is.na(GOID), !is.na(ONTOLOGY)) %>%
  distinct(GOID, ONTOLOGY)
colnames(go_ontology) <- c("term", "ONTOLOGY")

term2gene_go <- term2gene_go %>% left_join(go_ontology, by = "term")

topN              <- 20
min_count_for_plot <- 1

make_plots <- function(df, title_prefix, file_prefix) {
  if (nrow(df) == 0) return(invisible(NULL))
  if (!"Count" %in% colnames(df) && "geneID" %in% colnames(df)) {
    df$Count <- lengths(strsplit(as.character(df$geneID), "/"))
  }

  plot_df <- df %>%
    mutate(
      Description = str_trim(Description),
      Description = ifelse(is.na(Description) | Description == "", ID, Description)
    ) %>%
    group_by(Description) %>%
    summarise(
      Count      = max(as.numeric(Count), na.rm = TRUE),
      GeneRatio  = max(parse_ratio(GeneRatio), na.rm = TRUE),
      p.adjust   = min(as.numeric(p.adjust), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(Count >= min_count_for_plot) %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = topN)

  write_csv(plot_df, file.path(out_dir, paste0(file_prefix, "_plot_table_top.csv")))
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df$p.adjust[!is.finite(plot_df$p.adjust) | is.na(plot_df$p.adjust)] <- 1
  min_pos <- safe_min_positive(plot_df$p.adjust)
  plot_df$p.plot <- ifelse(plot_df$p.adjust <= 0 | is.na(plot_df$p.adjust), min_pos, plot_df$p.adjust)
  plot_df$Term   <- factor(str_wrap(plot_df$Description, 45),
                            levels = rev(unique(str_wrap(plot_df$Description, 45))))

  p_dot <- ggplot(plot_df, aes(GeneRatio, Term, size = Count, color = p.plot)) +
    geom_point(position = position_jitter(height = 0.08, width = 0)) +
    scale_color_gradient(low = "#D73027", high = "#4575B4", trans = "log10",
                         name = "p.adjust", labels = label_scientific(digits = 1)) +
    scale_size_continuous(name = "Count") +
    labs(title = paste0(title_prefix, " - dotplot (all timepoints)"), x = "GeneRatio", y = NULL) +
    theme_bw(base_size = 13)
  ggsave(file.path(out_dir, paste0(file_prefix, "_dotplot.png")), p_dot, width = 10, height = 7, dpi = 300)

  p_bar <- ggplot(plot_df, aes(Count, Term, fill = p.plot)) +
    geom_col() +
    scale_fill_gradient(low = "#4575B4", high = "#D73027", trans = "log10",
                        name = "p.adjust", labels = label_scientific(digits = 1)) +
    labs(title = paste0(title_prefix, " - barplot (all timepoints)"), x = "Count", y = NULL) +
    theme_bw(base_size = 13)
  ggsave(file.path(out_dir, paste0(file_prefix, "_barplot.png")), p_bar, width = 10, height = 7, dpi = 300)
}

# All ontologies together
ego_all    <- enricher(gene = genes_sig, universe = genes_bg,
                       TERM2GENE = term2gene_go %>% dplyr::select(term, gene),
                       TERM2NAME = term2name_go,
                       pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
ego_all_df <- as.data.frame(ego_all)
write_csv(ego_all_df, file.path(out_dir, "GO_ALL_results.csv"))
make_plots(ego_all_df, "GO enrichment (all ontologies)", "GO_ALL")

# Per-ontology
run_ontology <- function(ontology_code, label) {
  t2g <- term2gene_go %>% filter(ONTOLOGY == ontology_code) %>% dplyr::select(term, gene)
  t2n <- term2name_go %>% filter(term %in% t2g$term)
  if (nrow(t2g) == 0) { message("No terms for ", label); return(invisible(NULL)) }

  ego    <- enricher(gene = genes_sig, universe = genes_bg,
                     TERM2GENE = t2g, TERM2NAME = t2n,
                     pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  ego_df <- as.data.frame(ego)
  write_csv(ego_df, file.path(out_dir, paste0("GO_", label, "_results.csv")))
  make_plots(ego_df, paste0("GO ", label), paste0("GO_", label))
  message("Finished GO ", label)
}

run_ontology("BP", "BP")
run_ontology("MF", "MF")
run_ontology("CC", "CC")

cat("Done. Files saved in:\n", out_dir, "\n")
