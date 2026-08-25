# ====================================================================
# Ravindu: Novogene 2026 Mice Analysis
# Original Tumor 
# Clustering, DEG, Pathway enrichment Expression of known Biomarkers
# ====================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(openxlsx)
  library(pheatmap)
  library(tibble)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(stringr)
  library(forcats)
  library(grid)
  library(gridExtra)
  library(mMCPcounter)
})

# ===================
# OUTPUT DIRECTORIES
# ===================
outdir <- "original_tumor_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

qc_dir <- file.path(outdir, "01_QC")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

deg_dir <- file.path(outdir, "02_DEG_analysis")
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)

deg_table_dir <- file.path(deg_dir, "01_DEG_tables")
dir.create(deg_table_dir, showWarnings = FALSE, recursive = TRUE)

deg_heatmap_dir <- file.path(deg_dir, "02_DEG_heatmaps")
dir.create(deg_heatmap_dir, showWarnings = FALSE, recursive = TRUE)

expr_plot_dir <- file.path(outdir, "03_Expression_boxplots")
dir.create(expr_plot_dir, showWarnings = FALSE, recursive = TRUE)

scatter_root_dir <- file.path(outdir, "04_Correlation_scatterplots")
dir.create(scatter_root_dir, showWarnings = FALSE, recursive = TRUE)

corr_bar_dir <- file.path(outdir, "05_Correlation_barplots")
dir.create(corr_bar_dir, showWarnings = FALSE, recursive = TRUE)

pearson_bar_dir <- file.path(outdir, "05b_Pearson_correlation_barplots")
dir.create(pearson_bar_dir, showWarnings = FALSE, recursive = TRUE)

enrich_dir <- file.path(outdir, "06_Pathway_analysis")
dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

go_bp_dir <- file.path(enrich_dir, "01_GO_BP_ORA")
dir.create(go_bp_dir, showWarnings = FALSE, recursive = TRUE)

hallmark_gsea_dir <- file.path(enrich_dir, "02_Hallmark_GSEA")
dir.create(hallmark_gsea_dir, showWarnings = FALSE, recursive = TRUE)

immune_gsea_dir <- file.path(enrich_dir, "03_ImmuneSigDB_GSEA")
dir.create(immune_gsea_dir, showWarnings = FALSE, recursive = TRUE)

mmcp_root_dir <- file.path(outdir, "07_mMCP_counter")
dir.create(mmcp_root_dir, showWarnings = FALSE, recursive = TRUE)

mmcp_input_dir <- file.path(mmcp_root_dir, "01_input")
mmcp_scores_dir <- file.path(mmcp_root_dir, "02_scores")
mmcp_group_dir <- file.path(mmcp_root_dir, "03_group_comparisons")
mmcp_plot_dir <- file.path(mmcp_root_dir, "04_plots")

dir.create(mmcp_input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(mmcp_scores_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(mmcp_group_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(mmcp_plot_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# LABELS / THEMES / COLORS
# =========================

tumor_label <- "original tumor"
comparison_label <- "KO original tumor vs WT original tumor"

theme_nature <- function(base_size = 14, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.8),
      strip.text = element_text(face = "bold", colour = "black")
    )
}
theme_set(theme_nature())

group_cols_box <- c("WT" = "#F8766D", "KO" = "#00BFC4")
group_cols_qc <- c("WT" = "#B22222", "KO" = "#0B3C8C")
corr_bar_cols <- c("Positive" = "#E64B35FF", "Negative" = "#4DBBD5FF")
pearson_bar_cols <- c("Positive" = "#E64B35FF", "Negative" = "#4DBBD5FF")

theme_enrich_medium <- function() {
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 7, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7, face = "bold")
  )
}

theme_enrich_small <- function() {
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 7, colour = "black"),
    axis.text.y = element_text(size = 4, colour = "black"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    strip.text = element_text(size = 7, face = "bold")
  )
}

short_sample_name <- function(x) {
  paste0(substr(x, 1, 7), "_", substr(x, nchar(x) - 1, nchar(x)))
}

# ========
# INPUT
# ========

count_file <- "../counts/Galaxy Column join on original.tabular"

raw <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# =============================================
# MERGE TECHNICAL REPLICATE LANES FOR T2552WT
# =============================================

t2552wt_cols <- c(
  "T2552WT_MKRN260030404-1A_23M5KWLT4_L1",
  "T2552WT_MKRN260030404-1A_23M72JLT4_L4",
  "T2552WT_MKRN260030404-1A_23M53NLT4_L4"
)

if (!all(t2552wt_cols %in% colnames(raw))) {
  missing_cols <- t2552wt_cols[!t2552wt_cols %in% colnames(raw)]
  stop(paste("Missing T2552WT columns:", paste(missing_cols, collapse = ", ")))
}

raw[["T2552WT_MKRN260030404-1A_23M5KWLT4_LX"]] <- rowSums(raw[, t2552wt_cols], na.rm = TRUE)
raw <- raw[, !colnames(raw) %in% t2552wt_cols]

gene_col <- 1
gene_id_raw <- as.character(raw[[gene_col]])
gene_id <- sub("\\..*$", "", gene_id_raw)

count_df <- raw[, -gene_col, drop = FALSE]
count_mat <- as.matrix(count_df)
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- gene_id

sample_names <- colnames(count_mat)

group <- ifelse(
  grepl("KO", sample_names, ignore.case = TRUE), "KO",
  ifelse(grepl("WT", sample_names, ignore.case = TRUE), "WT", NA)
)

group <- factor(group, levels = c("WT", "KO"))

sample_info <- data.frame(
  sample = sample_names,
  group = group,
  sample_type = tumor_label,
  sample_short = short_sample_name(sample_names),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_names

# =========================
# FILTER / DESEQ2 OBJECT
# =========================

keep_rows <- rowSums(count_mat, na.rm = TRUE) > 1
count_mat <- count_mat[keep_rows, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = sample_info,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

norm_mat <- counts(dds, normalized = TRUE)
vst_obj <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_obj)

write.csv(
  data.frame(gene_id = rownames(norm_mat), norm_mat, check.names = FALSE),
  file = file.path(outdir, "normalized_counts_original_tumor.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = rownames(vst_mat), vst_mat, check.names = FALSE),
  file = file.path(outdir, "vst_expression_matrix_original_tumor.csv"),
  row.names = FALSE
)

# =============
# ANNOTATION
# =============

anno <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = rownames(vst_mat),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
)

anno <- anno[!is.na(anno$ENSEMBL), ]
anno <- anno[!duplicated(anno$ENSEMBL), ]

row_map <- data.frame(
  ENSEMBL = rownames(vst_mat),
  stringsAsFactors = FALSE
) %>%
  left_join(anno, by = c("ENSEMBL"))

row_map$SYMBOL_UPPER <- toupper(row_map$SYMBOL)

# =========
# PCA QC
# =========

run_pca_plot <- function(mat, sample_info, remove_prefixes = NULL, title_text, file_prefix) {
  mat_use <- mat
  if (!is.null(remove_prefixes) && length(remove_prefixes) > 0) {
    remove_mask <- Reduce(
      `|`,
      lapply(remove_prefixes, function(pref) grepl(paste0("^", pref), colnames(mat_use)))
    )
    removed_samples <- colnames(mat_use)[remove_mask]
    kept_samples <- colnames(mat_use)[!remove_mask]
    mat_use <- mat_use[, kept_samples, drop = FALSE]
  } else {
    removed_samples <- character(0)
    kept_samples <- colnames(mat_use)
  }
  
  pca_obj <- prcomp(t(mat_use), center = TRUE, scale. = FALSE)
  percent_var <- (pca_obj$sdev^2 / sum(pca_obj$sdev^2)) * 100
  
  pca_df <- data.frame(
    sample = rownames(pca_obj$x),
    sample_short = short_sample_name(rownames(pca_obj$x)),
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    group = sample_info[rownames(pca_obj$x), "group"],
    sample_type = tumor_label,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    pca_df,
    file = file.path(qc_dir, paste0(file_prefix, "_coordinates.csv")),
    row.names = FALSE
  )
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
    geom_point(size = 4, alpha = 0.95) +
    stat_ellipse(
      aes(group = group),
      type = "norm",
      geom = "polygon",
      alpha = 0.15,
      colour = NA,
      show.legend = FALSE
    ) +
    scale_color_manual(values = group_cols_qc) +
    scale_fill_manual(values = group_cols_qc) +
    geom_text(
      aes(label = sample_short),
      vjust = -0.8,
      size = 2.8,
      colour = "black",
      check_overlap = TRUE
    ) +
    labs(
      title = paste0(title_text, " (", tumor_label, ")"),
      x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
      y = paste0("PC2 (", round(percent_var[2], 1), "%)")
    ) +
    theme_nature(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14),
      legend.position = "right"
    )
  
  ggsave(
    file.path(qc_dir, paste0(file_prefix, ".png")),
    p_pca,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  list(plot = p_pca, pca_df = pca_df, removed_samples = removed_samples, kept_samples = kept_samples)
}

pca_all <- run_pca_plot(
  mat = vst_mat,
  sample_info = sample_info,
  remove_prefixes = NULL,
  title_text = "PCA of mouse RNA-seq samples",
  file_prefix = "PCA_PC1_vs_PC2_all_samples_original_tumor"
)

pca_filtered <- run_pca_plot(
  mat = vst_mat,
  sample_info = sample_info,
  remove_prefixes = c("T2498WT", "T2552WT"),
  title_text = "PCA of mouse RNA-seq filtered",
  file_prefix = "PCA_PC1_vs_PC2_filtered_samples_original_tumor"
)

# ======================================
# SAMPLE-TO-SAMPLE CORRELATION HEATMAPS
# ======================================

sample_cor_mat <- cor(vst_mat, method = "spearman", use = "pairwise.complete.obs")
sample_short_labels <- short_sample_name(colnames(sample_cor_mat))

ann_full <- data.frame(group = sample_info$group, row.names = rownames(sample_info))

pheatmap(
  sample_cor_mat,
  annotation_col = ann_full,
  annotation_row = ann_full,
  labels_col = sample_short_labels,
  labels_row = sample_short_labels,
  fontsize_col = 9,
  fontsize_row = 9,
  main = paste0("Sample-to-sample Spearman correlation (", tumor_label, ")"),
  filename = file.path(qc_dir, "Sample_correlation_heatmap_original_tumor.png"),
  width = 9,
  height = 8
)

filtered_samples <- pca_filtered$kept_samples
sample_cor_mat_filtered <- cor(
  vst_mat[, filtered_samples, drop = FALSE],
  method = "spearman",
  use = "pairwise.complete.obs"
)

sample_short_labels_filtered <- short_sample_name(colnames(sample_cor_mat_filtered))
ann_filtered <- data.frame(
  group = sample_info[filtered_samples, "group", drop = TRUE],
  row.names = filtered_samples
)

pheatmap(
  sample_cor_mat_filtered,
  annotation_col = ann_filtered,
  annotation_row = ann_filtered,
  labels_col = sample_short_labels_filtered,
  labels_row = sample_short_labels_filtered,
  fontsize_col = 9,
  fontsize_row = 9,
  main = paste0("Sample-to-sample Spearman correlation filtered (", tumor_label, ")"),
  filename = file.path(qc_dir, "Sample_correlation_heatmap_filtered_original_tumor.png"),
  width = 9,
  height = 8
)

# =============================
# SAMPLE CLUSTERING DENDROGRAM
# =============================

hc <- hclust(dist(t(vst_mat), method = "euclidean"), method = "complete")
hc$labels <- short_sample_name(hc$labels)

png(file.path(qc_dir, "Sample_clustering_dendrogram_original_tumor.png"),
    width = 2400, height = 1800, res = 300)
plot(
  hc,
  main = paste0("Sample clustering (rlog/vst, Euclidean) (", tumor_label, ")"),
  xlab = "",
  sub = "",
  cex = 0.9
)
dev.off()

# ==============
# DEG ANALYSIS
# ==============

res <- results(dds, contrast = c("group", "KO", "WT"), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  rownames_to_column("ENSEMBL") %>%
  left_join(row_map %>% select(ENSEMBL, SYMBOL, ENTREZID, GENENAME), by = "ENSEMBL") %>%
  arrange(padj, desc(abs(log2FoldChange)))

res_df$regulated <- case_when(
  !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0.58 ~ "Up_in_KO",
  !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.58 ~ "Down_in_KO",
  TRUE ~ "Not_significant"
)

deg_sig <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.58)
deg_top50 <- res_df %>% filter(!is.na(padj)) %>% arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n = 50)
deg_top100 <- res_df %>% filter(!is.na(padj)) %>% arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n = 100)

deg_summary <- data.frame(
  Comparison = "KO_original_tumor_vs_WT_original_tumor",
  Total_DE_Genes = nrow(deg_sig),
  Upregulated = sum(deg_sig$log2FoldChange > 0.58, na.rm = TRUE),
  Downregulated = sum(deg_sig$log2FoldChange < -0.58, na.rm = TRUE),
  padj_Cutoff = 0.05,
  log2FC_Cutoff = 0.58,
  stringsAsFactors = FALSE
)

write.csv(res_df, file.path(deg_table_dir, "DESeq2_all_results_KO_vs_WT_original_tumor.csv"), row.names = FALSE)
write.csv(deg_sig, file.path(deg_table_dir, "DESeq2_significant_DEGs_KO_vs_WT_original_tumor.csv"), row.names = FALSE)
write.csv(deg_top50, file.path(deg_table_dir, "DESeq2_top50_ranked_genes_KO_vs_WT_original_tumor.csv"), row.names = FALSE)
write.csv(deg_top100, file.path(deg_table_dir, "DESeq2_top100_ranked_genes_KO_vs_WT_original_tumor.csv"), row.names = FALSE)
write.csv(deg_summary, file.path(deg_table_dir, "DEG_summary_KO_vs_WT_original_tumor.csv"), row.names = FALSE)

deg_summary_tbl <- gridExtra::tableGrob(
  deg_summary,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(base_size = 9)
)

png(file.path(deg_dir, "DEG_summary_table_original_tumor.png"), width = 3600, height = 800, res = 300)
grid.newpage()
grid.draw(deg_summary_tbl)
dev.off()

# ==============
# DEG HEATMAPS
# ==============

plot_deg_heatmap <- function(gene_table, mat, sample_info, filename, title_text, scale_rows = TRUE) {
  genes_use <- intersect(gene_table$ENSEMBL, rownames(mat))
  if (length(genes_use) < 2) return(NULL)
  
  hm_mat <- mat[genes_use, , drop = FALSE]
  gene_symbols <- gene_table$SYMBOL[match(rownames(hm_mat), gene_table$ENSEMBL)]
  
  rownames(hm_mat) <- gene_symbols
  rownames(hm_mat)[is.na(rownames(hm_mat)) | rownames(hm_mat) == ""] <- genes_use[is.na(rownames(hm_mat)) | rownames(hm_mat) == ""]
  rownames(hm_mat) <- make.unique(rownames(hm_mat))
  
  ann_col <- data.frame(
    group = sample_info[colnames(hm_mat), "group"],
    sample_type = tumor_label,
    row.names = colnames(hm_mat)
  )
  
  pheatmap(
    hm_mat,
    scale = if (scale_rows) "row" else "none",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = ann_col,
    labels_col = short_sample_name(colnames(hm_mat)),
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = ifelse(nrow(hm_mat) > 80, 5, 7),
    fontsize_col = 9,
    angle_col = 45,
    main = paste0(title_text, " (", tumor_label, ")"),
    filename = filename,
    width = 10,
    height = 12
  )
}

plot_deg_heatmap(deg_sig, vst_mat, sample_info,
                 file.path(deg_heatmap_dir, "Heatmap_all_significant_DEGs_KO_vs_WT_original_tumor.png"),
                 "Heatmap of all significant DEGs")
plot_deg_heatmap(deg_top50, vst_mat, sample_info,
                 file.path(deg_heatmap_dir, "Heatmap_top50_DEGs_KO_vs_WT_original_tumor.png"),
                 "Heatmap of top 50 DEGs")
plot_deg_heatmap(deg_top100, vst_mat, sample_info,
                 file.path(deg_heatmap_dir, "Heatmap_top100_DEGs_KO_vs_WT_original_tumor.png"),
                 "Heatmap of top 100 DEGs")

# =================
# BIOMARKER PANEL
# =================

biomarkers <- c(
  "KLK3", "AR", "NKX3-1", "TMPRSS2", "AMACR", "ERG", "MKI67", "CDH1", "VIM",
  "AURKA", "BRCA1", "BRCA2", "MYC", "TP53", "PTEN", "RB1", "TTF1", "INSM1",
  "NKX2-1", "ACP3", "CHGA", "CHGB", "TFRC", "NCAM1", "SCG2", "SYP", "FOLH1"
)

biomarkers_upper <- toupper(biomarkers)
biomarker_presence <- data.frame(
  gene = biomarkers,
  found_in_mouse_annotation = biomarkers_upper %in% row_map$SYMBOL_UPPER,
  stringsAsFactors = FALSE
)

write.csv(
  biomarker_presence,
  file = file.path(outdir, "biomarker_presence_check_original_tumor.csv"),
  row.names = FALSE
)

available_biomarkers <- biomarker_presence$gene[biomarker_presence$found_in_mouse_annotation]

get_symbol_rows <- function(symbol) {
  idx <- which(row_map$SYMBOL_UPPER == toupper(symbol))
  if (length(idx) == 0) return(integer(0))
  match(row_map$ENSEMBL[idx], rownames(vst_mat))
}

if (!("FOLH1" %in% available_biomarkers)) {
  stop("FOLH1 was not found in the mouse annotation.")
}

gene_expr_df <- lapply(available_biomarkers, function(g) {
  rows <- get_symbol_rows(g)
  x <- as.numeric(colMeans(vst_mat[rows, , drop = FALSE], na.rm = TRUE))
  data.frame(
    sample = colnames(vst_mat),
    group = sample_info[colnames(vst_mat), "group"],
    sample_type = tumor_label,
    gene = g,
    expr = x,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(gene_expr_df, file = file.path(outdir, "biomarker_expression_long_VST_original_tumor.csv"), row.names = FALSE)

plot_group_box <- function(df, gene, ylab = "VST expression") {
  wilcox <- tryCatch(wilcox.test(expr ~ group, data = df), error = function(e) NULL)
  p_label <- if (!is.null(wilcox)) paste0("Wilcoxon p = ", signif(wilcox$p.value, 3)) else "Wilcoxon p = NA"
  y_range <- range(df$expr, na.rm = TRUE)
  y_pad <- ifelse(diff(y_range) == 0, 0.5, 0.15 * diff(y_range))
  y_pos <- y_range[2] + 0.5 * y_pad
  
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.9, linewidth = 0.9) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.55, color = "black") +
    annotate("text", x = 1.5, y = y_pos, label = p_label, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = group_cols_box) +
    labs(
      title = paste0(gene, " (", tumor_label, ")"),
      x = NULL,
      y = ylab
    ) +
    coord_cartesian(ylim = c(y_range[1], y_range[2] + y_pad)) +
    theme_nature(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 15),
      legend.position = "none"
    )
  
  ggsave(file.path(expr_plot_dir, paste0(gene, "_KO_vs_WT_boxplot_original_tumor.png")),
         p, width = 6, height = 5.5, dpi = 300)
  
  if (!is.null(wilcox)) {
    data.frame(gene = gene, p.value = wilcox$p.value)
  } else {
    data.frame(gene = gene, p.value = NA_real_)
  }
}

group_test_results <- lapply(available_biomarkers, function(g) {
  df_g <- gene_expr_df %>% filter(gene == g)
  plot_group_box(df_g, g)
}) %>% bind_rows()

group_test_results$p.adj <- p.adjust(group_test_results$p.value, method = "BH")
write.csv(group_test_results, file = file.path(outdir, "KO_vs_WT_wilcox_biomarkers_original_tumor.csv"), row.names = FALSE)

group_summary <- gene_expr_df %>%
  group_by(gene, group) %>%
  summarise(
    n = sum(is.finite(expr)),
    mean_expr = mean(expr, na.rm = TRUE),
    median_expr = median(expr, na.rm = TRUE),
    sd_expr = sd(expr, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(group_summary, file = file.path(outdir, "biomarker_group_summary_original_tumor.csv"), row.names = FALSE)

# ====================
# FOLH1 CORRELATIONS
# ====================

calc_group_correlations <- function(group_name) {
  samples_in_group <- rownames(sample_info)[sample_info$group == group_name]
  folh1_rows <- get_symbol_rows("FOLH1")
  folh1_expr_group <- as.numeric(colMeans(vst_mat[folh1_rows, samples_in_group, drop = FALSE], na.rm = TRUE))
  
  res <- lapply(available_biomarkers[available_biomarkers != "FOLH1"], function(g) {
    rows <- get_symbol_rows(g)
    x <- as.numeric(colMeans(vst_mat[rows, samples_in_group, drop = FALSE], na.rm = TRUE))
    ok <- is.finite(x) & is.finite(folh1_expr_group)
    
    if (sum(ok) < 3) {
      return(data.frame(
        gene = g,
        group = group_name,
        cor = NA_real_,
        p.value = NA_real_,
        n = sum(ok),
        stringsAsFactors = FALSE
      ))
    }
    
    ct <- suppressWarnings(cor.test(x[ok], folh1_expr_group[ok], method = "spearman", exact = FALSE))
    data.frame(
      gene = g,
      group = group_name,
      cor = unname(ct$estimate),
      p.value = ct$p.value,
      n = sum(ok),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  res$p.adj <- p.adjust(res$p.value, method = "BH")
  res %>% arrange(desc(cor))
}

cor_results_wt <- calc_group_correlations("WT")
cor_results_ko <- calc_group_correlations("KO")
cor_results_all <- bind_rows(cor_results_wt, cor_results_ko)

write.csv(cor_results_wt, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_WT_original_tumor.csv"), row.names = FALSE)
write.csv(cor_results_ko, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_KO_original_tumor.csv"), row.names = FALSE)
write.csv(cor_results_all, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_by_group_original_tumor.csv"), row.names = FALSE)

plot_correlation_bar <- function(cor_df, group_name) {
  plot_df <- cor_df %>%
    filter(is.finite(cor)) %>%
    mutate(direction = ifelse(cor >= 0, "Positive", "Negative"), gene = factor(gene, levels = rev(gene)))
  
  p <- ggplot(plot_df, aes(x = cor, y = gene, fill = direction)) +
    geom_col(width = 0.76, colour = "white", linewidth = 0.5) +
    geom_vline(xintercept = 0, linewidth = 1, colour = "black") +
    scale_fill_manual(values = corr_bar_cols) +
    labs(
      title = paste0("Correlation of FOLH1 with known biomarkers (", group_name, ", ", tumor_label, ")"),
      x = "Spearman correlation with FOLH1",
      y = NULL
    ) +
    theme_classic(base_size = 18, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title = element_text(face = "bold", colour = "black", size = 18),
      axis.text = element_text(colour = "black", size = 15),
      axis.line = element_line(linewidth = 0.9, colour = "black"),
      axis.ticks = element_line(linewidth = 0.8, colour = "black"),
      legend.position = "none"
    )
  
  ggsave(file.path(corr_bar_dir, paste0("FOLH1_biomarker_correlation_barplot_", group_name, "_original_tumor.png")),
         p, width = 10, height = 8, dpi = 300)
}

plot_correlation_bar(cor_results_wt, "WT")
plot_correlation_bar(cor_results_ko, "KO")

# =============================
# PEARSON CORRELATION BARPLOTS
# =============================

calc_group_pearson <- function(group_name) {
  samples_in_group <- rownames(sample_info)[sample_info$group == group_name]
  folh1_rows <- get_symbol_rows("FOLH1")
  folh1_expr_group <- as.numeric(colMeans(vst_mat[folh1_rows, samples_in_group, drop = FALSE], na.rm = TRUE))
  
  res <- lapply(available_biomarkers[available_biomarkers != "FOLH1"], function(g) {
    rows <- get_symbol_rows(g)
    x <- as.numeric(colMeans(vst_mat[rows, samples_in_group, drop = FALSE], na.rm = TRUE))
    ok <- is.finite(x) & is.finite(folh1_expr_group)
    
    if (sum(ok) < 3) {
      return(data.frame(
        gene = g,
        group = group_name,
        cor = NA_real_,
        p.value = NA_real_,
        n = sum(ok),
        stringsAsFactors = FALSE
      ))
    }
    
    ct <- suppressWarnings(cor.test(x[ok], folh1_expr_group[ok], method = "pearson"))
    data.frame(
      gene = g,
      group = group_name,
      cor = unname(ct$estimate),
      p.value = ct$p.value,
      n = sum(ok),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  res$p.adj <- p.adjust(res$p.value, method = "BH")
  res %>% arrange(desc(cor))
}

pearson_results_wt <- calc_group_pearson("WT")
pearson_results_ko <- calc_group_pearson("KO")
pearson_results_all <- bind_rows(pearson_results_wt, pearson_results_ko)

write.csv(pearson_results_wt, file = file.path(outdir, "FOLH1_vs_biomarkers_pearson_WT_original_tumor.csv"), row.names = FALSE)
write.csv(pearson_results_ko, file = file.path(outdir, "FOLH1_vs_biomarkers_pearson_KO_original_tumor.csv"), row.names = FALSE)
write.csv(pearson_results_all, file = file.path(outdir, "FOLH1_vs_biomarkers_pearson_by_group_original_tumor.csv"), row.names = FALSE)

plot_pearson_bar <- function(cor_df, group_name) {
  plot_df <- cor_df %>%
    filter(is.finite(cor)) %>%
    mutate(direction = ifelse(cor >= 0, "Positive", "Negative"), gene = factor(gene, levels = rev(gene)))
  
  p <- ggplot(plot_df, aes(x = cor, y = gene, fill = direction)) +
    geom_col(width = 0.76, colour = "white", linewidth = 0.5) +
    geom_vline(xintercept = 0, linewidth = 1, colour = "black") +
    scale_fill_manual(values = pearson_bar_cols) +
    labs(
      title = paste0("Pearson correlation of FOLH1 with biomarkers (", group_name, ", ", tumor_label, ")"),
      x = "Pearson correlation with FOLH1",
      y = NULL
    ) +
    theme_classic(base_size = 18, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title = element_text(face = "bold", colour = "black", size = 18),
      axis.text = element_text(colour = "black", size = 15),
      axis.line = element_line(linewidth = 0.9, colour = "black"),
      axis.ticks = element_line(linewidth = 0.8, colour = "black"),
      legend.position = "none"
    )
  
  ggsave(file.path(pearson_bar_dir, paste0("FOLH1_biomarker_pearson_barplot_", group_name, "_original_tumor.png")),
         p, width = 10, height = 8, dpi = 300)
}

plot_pearson_bar(pearson_results_wt, "WT")
plot_pearson_bar(pearson_results_ko, "KO")

# ===============================================
# COMBINED WT+KO SCATTERPLOTS FOR EACH BIOMARKER
# ===============================================
scatter_fun_group <- function(g, save_dir = scatter_root_dir) {
  samples_use <- rownames(sample_info)
  
  folh1_rows <- get_symbol_rows("FOLH1")
  folh1_expr <- as.numeric(colMeans(vst_mat[folh1_rows, samples_use, drop = FALSE], na.rm = TRUE))
  
  rows <- get_symbol_rows(g)
  biomarker_expr <- as.numeric(colMeans(vst_mat[rows, samples_use, drop = FALSE], na.rm = TRUE))
  
  df <- data.frame(
    sample = samples_use,
    FOLH1 = folh1_expr,
    biomarker = biomarker_expr,
    group = sample_info[samples_use, "group"],
    sample_type = tumor_label,
    stringsAsFactors = FALSE
  )
  
  df <- df[is.finite(df$FOLH1) & is.finite(df$biomarker), , drop = FALSE]
  if (nrow(df) < 3) return(NULL)
  
  ct_all <- suppressWarnings(cor.test(df$FOLH1, df$biomarker, method = "spearman", exact = FALSE))
  lab_all <- paste0(
    "All samples: Spearman rho = ",
    round(unname(ct_all$estimate), 3),
    ", p = ",
    signif(ct_all$p.value, 3)
  )
  
  p <- ggplot(df, aes(x = FOLH1, y = biomarker, color = group)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_smooth(aes(group = group, fill = group), method = "lm", se = FALSE, linewidth = 0.85) +
    scale_color_manual(values = group_cols_box) +
    scale_fill_manual(values = group_cols_box) +
    labs(
      title = paste0("FOLH1 vs ", g, " (", tumor_label, ")"),
      subtitle = lab_all,
      x = "FOLH1 (VST)",
      y = paste0(g, " (VST)")
    ) +
    theme_nature(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  ggsave(
    file.path(save_dir, paste0("FOLH1_vs_", g, "_scatter_combined_WT_KO_original_tumor.png")),
    p,
    width = 6.4,
    height = 5.4,
    dpi = 300
  )
  
  invisible(df)
}

for (g in available_biomarkers[available_biomarkers != "FOLH1"]) {
  scatter_fun_group(g, scatter_root_dir)
}

heat_genes <- available_biomarkers
heat_rows <- unlist(lapply(heat_genes, get_symbol_rows))
heat_rows <- unique(heat_rows)
heat_mat <- vst_mat[heat_rows, , drop = FALSE]

heat_anno_col <- data.frame(group = sample_info$group)
rownames(heat_anno_col) <- rownames(sample_info)

rownames(heat_mat) <- row_map$SYMBOL[match(rownames(heat_mat), row_map$ENSEMBL)]
rownames(heat_mat)[is.na(rownames(heat_mat)) | rownames(heat_mat) == ""] <- heat_rows[is.na(rownames(heat_mat)) | rownames(heat_mat) == ""]

pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = heat_anno_col,
  labels_col = short_sample_name(colnames(heat_mat)),
  show_colnames = TRUE,
  main = paste0("Mouse biomarker panel and FOLH1 (", tumor_label, ")"),
  filename = file.path(outdir, "biomarker_heatmap_original_tumor.png"),
  width = 10,
  height = 8
)
# =================
# ENRICHMENT PREP
# =================

res_ranked <- res_df %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(rank_metric = log2FoldChange) %>%
  arrange(desc(rank_metric))

universe_entrez <- res_df %>%
  filter(!is.na(ENTREZID), !is.na(baseMean)) %>%
  pull(ENTREZID) %>%
  unique()

rank_df <- res_ranked %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(ENTREZID) %>%
  summarise(rank_metric = rank_metric[which.max(abs(rank_metric))], .groups = "drop") %>%
  arrange(desc(rank_metric))

gene_list <- rank_df$rank_metric
names(gene_list) <- rank_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

sig_entrez_up <- deg_sig %>%
  filter(log2FoldChange > 0.58, !is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

sig_entrez_down <- deg_sig %>%
  filter(log2FoldChange < -0.58, !is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

# ===========
# GO BP ORA
# ===========

ego_up <- if (length(sig_entrez_up) > 0) {
  enrichGO(
    gene = sig_entrez_up,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
} else NULL

ego_down <- if (length(sig_entrez_down) > 0) {
  enrichGO(
    gene = sig_entrez_down,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
} else NULL

if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
  write.csv(as.data.frame(ego_up), file.path(go_bp_dir, "GO_BP_ORA_up_in_KO_original_tumor.csv"), row.names = FALSE)
  png(file.path(go_bp_dir, "GO_BP_ORA_up_in_KO_dotplot_original_tumor.png"), width = 2400, height = 1800, res = 300)
  print(dotplot(ego_up, showCategory = 20) + ggtitle(paste0("GO BP enriched in upregulated genes (", tumor_label, ")")) + theme_enrich_medium())
  dev.off()
}

if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
  write.csv(as.data.frame(ego_down), file.path(go_bp_dir, "GO_BP_ORA_down_in_KO_original_tumor.csv"), row.names = FALSE)
  png(file.path(go_bp_dir, "GO_BP_ORA_down_in_KO_dotplot_original_tumor.png"), width = 2400, height = 1800, res = 300)
  print(dotplot(ego_down, showCategory = 20) + ggtitle(paste0("GO BP enriched in downregulated genes (", tumor_label, ")")) + theme_enrich_medium())
  dev.off()
}

# =============
# MSIGDB SETS
# =============

msig_hallmark_mm <- msigdbr(species = "Mus musculus", category = "H")
hallmark_term2gene <- msig_hallmark_mm %>% select(gs_name, entrez_gene)

msig_immune_mm <- msigdbr(species = "Mus musculus", category = "C7")
immune_term2gene <- msig_immune_mm %>% select(gs_name, entrez_gene)

# ==============
# HALLMARK GSEA
# ==============

gsea_hallmark <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
  write.csv(as.data.frame(gsea_hallmark), file.path(hallmark_gsea_dir, "GSEA_Hallmark_ranked_original_tumor.csv"), row.names = FALSE)
  png(file.path(hallmark_gsea_dir, "GSEA_Hallmark_dotplot_original_tumor.png"), width = 2400, height = 1800, res = 300)
  print(dotplot(gsea_hallmark, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign) + ggtitle(paste0("Hallmark GSEA ranked analysis (", tumor_label, ")")) + theme_enrich_small())
  dev.off()
}

# ==================
# IMMUNESIGDB GSEA
# ==================

gsea_immune <- GSEA(
  geneList = gene_list,
  TERM2GENE = immune_term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (!is.null(gsea_immune) && nrow(as.data.frame(gsea_immune)) > 0) {
  write.csv(as.data.frame(gsea_immune), file.path(immune_gsea_dir, "GSEA_ImmuneSigDB_ranked_original_tumor.csv"), row.names = FALSE)
  png(file.path(immune_gsea_dir, "GSEA_ImmuneSigDB_dotplot_original_tumor.png"), width = 2400, height = 1800, res = 300)
  print(dotplot(gsea_immune, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign) + ggtitle(paste0("ImmuneSigDB GSEA ranked analysis (", tumor_label, ")")) + theme_enrich_small())
  dev.off()
}

# ==================
# mMCP-Counter PREP
# ==================

mmcp_map <- row_map %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  select(ENSEMBL, SYMBOL) %>%
  distinct()

mmcp_expr <- data.frame(ENSEMBL = rownames(norm_mat), norm_mat, check.names = FALSE) %>%
  left_join(mmcp_map, by = "ENSEMBL") %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  select(-ENSEMBL) %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

mmcp_mat <- as.data.frame(mmcp_expr)
rownames(mmcp_mat) <- mmcp_mat$SYMBOL
mmcp_mat$SYMBOL <- NULL
mmcp_mat <- as.matrix(mmcp_mat)

write.csv(
  data.frame(gene_symbol = rownames(mmcp_mat), mmcp_mat, check.names = FALSE),
  file = file.path(mmcp_input_dir, "mMCP_counter_input_symbol_matrix.csv"),
  row.names = FALSE
)

# =================
# RUN mMCP-counter
# ================

if (!"mMCPcounter" %in% rownames(installed.packages())) {
  stop("Package 'mMCPcounter' is not installed. Run: remotes::install_github('cit-bioinfo/mMCP-counter')")
}

mmcp_res <- mMCPcounter.estimate(
  exp = mmcp_mat,
  features = "Gene.Symbol",
  genomeVersion = "GCRm39"
)

mmcp_scores <- as.data.frame(mmcp_res)
mmcp_scores$cell_type <- rownames(mmcp_scores)

mmcp_scores_long <- mmcp_scores %>%
  pivot_longer(
    cols = -cell_type,
    names_to = "sample",
    values_to = "score"
  )

write.csv(mmcp_scores_long, file = file.path(mmcp_scores_dir, "mMCP_counter_scores_long.csv"), row.names = FALSE)

mmcp_wide <- mmcp_scores_long %>%
  pivot_wider(names_from = sample, values_from = score)

write.csv(mmcp_wide, file = file.path(mmcp_scores_dir, "mMCP_counter_scores_wide.csv"), row.names = FALSE)

# ================================
# mMCP-Counter GROUP COMPARISONS
# ================================

mmcp_scores2 <- mmcp_scores_long %>%
  left_join(sample_info, by = "sample")

mmcp_group_summary <- mmcp_scores2 %>%
  group_by(cell_type, group) %>%
  summarise(
    n = sum(is.finite(score)),
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(mmcp_group_summary, file = file.path(mmcp_group_dir, "mMCP_group_summary.csv"), row.names = FALSE)

pairwise_groups <- list(
  c("WT", "KO")
)

mmcp_pairwise_results <- mmcp_scores2 %>%
  group_by(cell_type) %>%
  summarise(
    p.value = tryCatch(wilcox.test(score ~ group)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    comparison = "WT_vs_KO",
    group1 = "WT",
    group2 = "KO",
    p.adj = p.adjust(p.value, method = "BH")
  )

write.csv(mmcp_pairwise_results, file = file.path(mmcp_group_dir, "mMCP_pairwise_tests.csv"), row.names = FALSE)

# ====================
# mMCP-Counter PLOTS
# ===================

plot_mmcp_box <- function(cell_type_name) {
  df <- mmcp_scores2 %>% filter(cell_type == cell_type_name)
  if (nrow(df) == 0) return(NULL)
  
  wtko_test <- tryCatch(
    wilcox.test(score ~ group, data = df),
    error = function(e) NULL
  )
  
  p_lab <- if (!is.null(wtko_test)) {
    paste0("Wilcoxon p = ", signif(wtko_test$p.value, 3))
  } else {
    "Wilcoxon p = NA"
  }
  
  y_range <- range(df$score, na.rm = TRUE)
  y_pad <- ifelse(diff(y_range) == 0, 0.5, 0.15 * diff(y_range))
  y_pos <- y_range[2] + 0.45 * y_pad
  
  p <- ggplot(df, aes(x = group, y = score, fill = group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.9, linewidth = 0.85) +
    geom_jitter(width = 0.12, size = 2.0, alpha = 0.65, colour = "black") +
    annotate("text", x = 1.5, y = y_pos, label = p_lab, size = 4.2, fontface = "bold") +
    scale_fill_manual(values = group_cols_box) +
    labs(
      title = paste0("mMCP-counter: ", cell_type_name, " (", tumor_label, ")"),
      x = NULL,
      y = "Score"
    ) +
    coord_cartesian(ylim = c(y_range[1], y_range[2] + y_pad)) +
    theme_nature(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
  
  ggsave(
    file.path(mmcp_plot_dir, paste0("mMCP_boxplot_", gsub("[^A-Za-z0-9]+", "_", cell_type_name), ".png")),
    p,
    width = 7.4,
    height = 5.8,
    dpi = 300
  )
}

mmcp_cell_types <- sort(unique(mmcp_scores2$cell_type))
for (ct in mmcp_cell_types) plot_mmcp_box(ct)

mmcp_heat <- mmcp_scores2 %>%
  select(cell_type, sample, score) %>%
  pivot_wider(names_from = sample, values_from = score) %>%
  as.data.frame()

rownames(mmcp_heat) <- mmcp_heat$cell_type
mmcp_heat$cell_type <- NULL
mmcp_heat <- as.matrix(mmcp_heat)

mmcp_ann_col <- data.frame(group = sample_info[colnames(mmcp_heat), "group"])
rownames(mmcp_ann_col) <- colnames(mmcp_heat)

pheatmap(
  mmcp_heat,
  scale = "row",
  annotation_col = mmcp_ann_col,
  main = paste0("mMCP-counter scores (", tumor_label, ")"),
  fontsize_col = 8,
  fontsize_row = 8,
  filename = file.path(mmcp_plot_dir, "mMCP_heatmap_all_samples.png"),
  width = 10,
  height = 8
)

mmcp_pca_obj <- prcomp(t(mmcp_heat), scale. = TRUE)
mmcp_pct <- (mmcp_pca_obj$sdev^2 / sum(mmcp_pca_obj$sdev^2)) * 100
mmcp_pca <- data.frame(
  sample = rownames(mmcp_pca_obj$x),
  PC1 = mmcp_pca_obj$x[,1],
  PC2 = mmcp_pca_obj$x[,2],
  group = sample_info[rownames(mmcp_pca_obj$x), "group"],
  stringsAsFactors = FALSE
)

p_mmcp_pca <- ggplot(mmcp_pca, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  scale_color_manual(values = group_cols_box) +
  labs(
    title = paste0("PCA of mMCP-counter scores (", tumor_label, ")"),
    x = paste0("PC1 (", round(mmcp_pct[1], 1), "%)"),
    y = paste0("PC2 (", round(mmcp_pct[2], 1), "%)")
  ) +
  theme_nature(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right")

ggsave(file.path(mmcp_plot_dir, "mMCP_PCA.png"), p_mmcp_pca, width = 7, height = 6, dpi = 300)

# =================
# EXCEL WORKBOOK
# =================

wb <- createWorkbook()

addWorksheet(wb, "Sample_info")
writeData(wb, "Sample_info", sample_info)

addWorksheet(wb, "PCA_QC_all")
writeData(wb, "PCA_QC_all", pca_all$pca_df)

addWorksheet(wb, "PCA_QC_filtered")
writeData(wb, "PCA_QC_filtered", pca_filtered$pca_df)

addWorksheet(wb, "DEG_summary")
writeData(wb, "DEG_summary", deg_summary)

addWorksheet(wb, "DESeq2_all_results")
writeData(wb, "DESeq2_all_results", res_df)

addWorksheet(wb, "DESeq2_sig_DEGs")
writeData(wb, "DESeq2_sig_DEGs", deg_sig)

addWorksheet(wb, "DESeq2_top50")
writeData(wb, "DESeq2_top50", deg_top50)

addWorksheet(wb, "DESeq2_top100")
writeData(wb, "DESeq2_top100", deg_top100)

addWorksheet(wb, "Biomarker_presence")
writeData(wb, "Biomarker_presence", biomarker_presence)

addWorksheet(wb, "Group_summary")
writeData(wb, "Group_summary", group_summary)

addWorksheet(wb, "Wilcox_results")
writeData(wb, "Wilcox_results", group_test_results)

addWorksheet(wb, "Spearman_WT")
writeData(wb, "Spearman_WT", cor_results_wt)

addWorksheet(wb, "Spearman_KO")
writeData(wb, "Spearman_KO", cor_results_ko)

addWorksheet(wb, "Pearson_WT")
writeData(wb, "Pearson_WT", pearson_results_wt)

addWorksheet(wb, "Pearson_KO")
writeData(wb, "Pearson_KO", pearson_results_ko)

addWorksheet(wb, "Expression_long")
writeData(wb, "Expression_long", gene_expr_df)

if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
  addWorksheet(wb, "GO_BP_up")
  writeData(wb, "GO_BP_up", as.data.frame(ego_up))
}
if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
  addWorksheet(wb, "GO_BP_down")
  writeData(wb, "GO_BP_down", as.data.frame(ego_down))
}
if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
  addWorksheet(wb, "GSEA_Hallmark")
  writeData(wb, "GSEA_Hallmark", as.data.frame(gsea_hallmark))
}
if (!is.null(gsea_immune) && nrow(as.data.frame(gsea_immune)) > 0) {
  addWorksheet(wb, "GSEA_ImmuneSigDB")
  writeData(wb, "GSEA_ImmuneSigDB", as.data.frame(gsea_immune))
}

addWorksheet(wb, "mMCP_scores_long")
writeData(wb, "mMCP_scores_long", mmcp_scores_long)

addWorksheet(wb, "mMCP_scores_wide")
writeData(wb, "mMCP_scores_wide", mmcp_wide)

addWorksheet(wb, "mMCP_group_summary")
writeData(wb, "mMCP_group_summary", mmcp_group_summary)

addWorksheet(wb, "mMCP_pairwise")
writeData(wb, "mMCP_pairwise", mmcp_pairwise_results)

saveWorkbook(
  wb,
  file = file.path(outdir, "Overall_mouse_FOLH1_biomarker_and_DEG_results_original_tumor.xlsx"),
  overwrite = TRUE
)
