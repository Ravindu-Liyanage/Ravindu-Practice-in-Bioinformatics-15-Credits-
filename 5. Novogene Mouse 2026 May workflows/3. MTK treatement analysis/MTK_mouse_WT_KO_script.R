# ====================================================================
# Ravindu: Novogene 2026 Mice Analysis
# MTYK Treatement 
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
  library(grid)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(rstatix)
})

# OUTPUT DIRECTORIES
outdir <- "mouse_FOLH1_biomarker_analysis_MTK_treated"
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

scatter_control_dir <- file.path(scatter_root_dir, "control")
scatter_treated_dir <- file.path(scatter_root_dir, "treated")
dir.create(scatter_control_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(scatter_treated_dir, showWarnings = FALSE, recursive = TRUE)

corr_bar_dir <- file.path(outdir, "05_Correlation_barplots")
dir.create(corr_bar_dir, showWarnings = FALSE, recursive = TRUE)

enrich_dir <- file.path(outdir, "06_Pathway_analysis")
dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

go_bp_dir <- file.path(enrich_dir, "01_GO_BP")
dir.create(go_bp_dir, showWarnings = FALSE, recursive = TRUE)

hallmark_dir <- file.path(enrich_dir, "02_Hallmark")
dir.create(hallmark_dir, showWarnings = FALSE, recursive = TRUE)

gsea_dir <- file.path(enrich_dir, "03_GSEA_ranked")
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

immune_dir <- file.path(enrich_dir, "04_ImmuneSigDB")
dir.create(immune_dir, showWarnings = FALSE, recursive = TRUE)

# LABELS / THEME

plot_label <- "MTK treated"

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

# ==========
# INPUT
# ==========

count_file <- "../counts/Galaxy Column join on MTK.tabular"

raw <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_col <- 1
gene_id_raw <- as.character(raw[[gene_col]])
gene_id <- sub("\\..*$", "", gene_id_raw)

count_df <- raw[, -gene_col, drop = FALSE]
count_mat <- as.matrix(count_df)
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- gene_id

sample_names <- colnames(count_mat)

condition <- ifelse(
  substr(sample_names, 1, 1) == "C", "control",
  ifelse(substr(sample_names, 1, 1) == "M", "treated", NA)
)

condition <- factor(
  condition,
  levels = c("treated", "control")
)

sample_info <- data.frame(
  sample = sample_names,
  sample_short = substr(sample_names, 1, 7),
  sample_short10 = substr(sample_names, 1, 10),
  condition = condition,
  group = condition,
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_names

# ==================
# FILTER / DESEQ2
# ==================

keep_rows <- rowSums(count_mat, na.rm = TRUE) > 1
count_mat <- count_mat[keep_rows, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = sample_info,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

norm_mat <- counts(dds, normalized = TRUE)
vst_obj <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_obj)

write.csv(
  data.frame(gene_id = rownames(norm_mat), norm_mat, check.names = FALSE),
  file = file.path(outdir, "normalized_counts_MTK_treated.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = rownames(vst_mat), vst_mat, check.names = FALSE),
  file = file.path(outdir, "vst_expression_matrix_MTK_treated.csv"),
  row.names = FALSE
)

# =========
# PCA QC
# =========

pca_obj <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
percent_var <- (pca_obj$sdev^2 / sum(pca_obj$sdev^2)) * 100

pca_df <- data.frame(
  sample = rownames(pca_obj$x),
  sample_short = substr(rownames(pca_obj$x), 1, 7),
  sample_short10 = substr(rownames(pca_obj$x), 1, 10),
  PC1 = pca_obj$x[, 1],
  PC2 = pca_obj$x[, 2],
  condition = sample_info[rownames(pca_obj$x), "condition"],
  group = sample_info[rownames(pca_obj$x), "group"],
  stringsAsFactors = FALSE
)

write.csv(
  pca_df,
  file = file.path(qc_dir, "PCA_coordinates_PC1_PC2_MTK_treated.csv"),
  row.names = FALSE
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, fill = condition)) +
  geom_point(size = 4, alpha = 0.95) +
  stat_ellipse(
    aes(group = condition),
    type = "norm",
    geom = "polygon",
    alpha = 0.12,
    colour = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("control" = "#F8766D", "treated" = "#00BFC4")) +
  scale_fill_manual(values = c("control" = "#F8766D", "treated" = "#00BFC4")) +
  geom_text(
    aes(label = sample_short),
    vjust = -0.8,
    size = 2.0,
    colour = "black",
    check_overlap = TRUE
  ) +
  labs(
    title = paste0("PCA of mouse RNA-seq samples (", plot_label, ")"),
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_nature(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.position = "right"
  )

ggsave(
  file.path(qc_dir, "PCA_PC1_vs_PC2_MTK_treated.png"),
  p_pca,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

sample_cor_mat <- cor(vst_mat, method = "spearman", use = "pairwise.complete.obs")

pheatmap(
  sample_cor_mat,
  annotation_col = data.frame(condition = sample_info$condition, row.names = rownames(sample_info)),
  annotation_row = data.frame(condition = sample_info$condition, row.names = rownames(sample_info)),
  labels_col = sample_info[colnames(sample_cor_mat), "sample_short10"],
  labels_row = sample_info[rownames(sample_cor_mat), "sample_short10"],
  main = paste0("Sample-to-sample Spearman correlation (", plot_label, ")"),
  fontsize_col = 8,
  fontsize_row = 8,
  filename = file.path(qc_dir, "Sample_correlation_heatmap_MTK_treated.png"),
  width = 8,
  height = 7
)

hc <- hclust(dist(t(vst_mat), method = "euclidean"), method = "complete")
png(file.path(qc_dir, "Sample_clustering_dendrogram_MTK_treated.png"), width = 2600, height = 2000, res = 300)
plot(
  hc,
  labels = sample_info[hc$labels, "sample_short10"],
  main = paste0("Sample clustering (rlog/vst, Euclidean) (", plot_label, ")"),
  xlab = "",
  sub = "",
  cex = 0.8
)
dev.off()

# ============
# ANNOTATION
# ============

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
  left_join(anno, by = "ENSEMBL")

row_map$SYMBOL_UPPER <- toupper(row_map$SYMBOL)

# ============
# BIOMARKERS
# ============

biomarkers <- c(
  "KLK3", "AR", "NKX3-1", "TMPRSS2", "AMACR", "ERG", "MKI67", "CDH1", "VIM",
  "AURKA", "BRCA1", "BRCA2", "MYC", "TP53", "PTEN", "RB1", "TTF1", "INSM1",
  "NKX2-1", "ACP3", "CHGA", "CHGB", "TFRC", "NCAM1", "SCG2", "SYP", "FOLH1"
)

biomarker_presence <- data.frame(
  gene = biomarkers,
  found_in_mouse_annotation = toupper(biomarkers) %in% row_map$SYMBOL_UPPER,
  stringsAsFactors = FALSE
)

write.csv(
  biomarker_presence,
  file = file.path(outdir, "biomarker_presence_check_MTK_treated.csv"),
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

# ======================
# EXPRESSION LONG TABLE
# ======================

gene_expr_df <- lapply(available_biomarkers, function(g) {
  rows <- get_symbol_rows(g)
  x <- as.numeric(colMeans(vst_mat[rows, , drop = FALSE], na.rm = TRUE))
  data.frame(
    sample = colnames(vst_mat),
    sample_short = sample_info[colnames(vst_mat), "sample_short"],
    sample_short10 = sample_info[colnames(vst_mat), "sample_short10"],
    condition = sample_info[colnames(vst_mat), "condition"],
    group = sample_info[colnames(vst_mat), "group"],
    gene = g,
    expr = x,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(
  gene_expr_df,
  file = file.path(outdir, "biomarker_expression_long_VST_MTK_treated.csv"),
  row.names = FALSE
)

# =========
# BOXPLOTS
# =========

plot_group_box <- function(df, gene, ylab = "VST expression") {
  
  wilcox_res <- tryCatch(
    wilcox.test(expr ~ condition, data = df),
    error = function(e) NULL
  )
  
  p_label <- if (!is.null(wilcox_res)) {
    paste0("Wilcoxon p = ", signif(wilcox_res$p.value, 3))
  } else {
    "Wilcoxon p = NA"
  }
  
  y_range <- range(df$expr, na.rm = TRUE)
  y_diff <- diff(y_range)
  if (y_diff == 0) y_diff <- 0.5
  y_pad <- 0.18 * y_diff
  y_pos <- y_range[2] + 0.45 * y_pad
  
  p <- ggplot(df, aes(x = condition, y = expr, fill = condition)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.9, linewidth = 0.9) +
    geom_jitter(width = 0.12, size = 2.0, alpha = 0.55, color = "black") +
    annotate(
      "text",
      x = 1.5,
      y = y_pos,
      label = p_label,
      size = 4.3,
      fontface = "bold"
    ) +
    scale_fill_manual(values = c("control" = "#F8766D", "treated" = "#00BFC4")) +
    labs(
      title = paste0(gene, " (", plot_label, ")"),
      subtitle = "Wilcoxon rank-sum test across 2 groups",
      x = NULL,
      y = ylab
    ) +
    coord_cartesian(ylim = c(y_range[1], y_range[2] + y_pad)) +
    theme_nature(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 15),
      legend.position = "none"
    )
  
  ggsave(
    file.path(expr_plot_dir, paste0(gene, "_group_boxplot_MTK_treated.png")),
    p,
    width = 6.4,
    height = 5.6,
    dpi = 300
  )
  
  wilcox_df <- if (!is.null(wilcox_res)) {
    data.frame(
      gene = gene,
      test = "Wilcoxon rank-sum",
      p.value = wilcox_res$p.value,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      gene = gene,
      test = "Wilcoxon rank-sum",
      p.value = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  wilcox_df
}

group_test_results <- lapply(available_biomarkers, function(g) {
  df_g <- gene_expr_df %>% dplyr::filter(gene == g)
  plot_group_box(df_g, g)
}) %>% bind_rows()

group_test_results$p.adj <- p.adjust(group_test_results$p.value, method = "BH")

write.csv(
  group_test_results,
  file = file.path(outdir, "group_wilcox_biomarkers_MTK_treated.csv"),
  row.names = FALSE
)

group_summary <- gene_expr_df %>%
  group_by(gene, condition) %>%
  summarise(
    n = sum(is.finite(expr)),
    mean_expr = mean(expr, na.rm = TRUE),
    median_expr = median(expr, na.rm = TRUE),
    sd_expr = sd(expr, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  group_summary,
  file = file.path(outdir, "biomarker_group_summary_MTK_treated.csv"),
  row.names = FALSE
)

# ==============
# CORRELATIONS
# ==============

calc_group_correlations <- function(group_name) {
  samples_in_group <- rownames(sample_info)[sample_info$condition == group_name]
  
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

cor_results_control <- calc_group_correlations("control")
cor_results_treated <- calc_group_correlations("treated")
cor_results_all <- bind_rows(cor_results_control, cor_results_treated)

write.csv(cor_results_control, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_control_MTK_treated.csv"), row.names = FALSE)
write.csv(cor_results_treated, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_treated_MTK_treated.csv"), row.names = FALSE)
write.csv(cor_results_all, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_by_group_MTK_treated.csv"), row.names = FALSE)

plot_correlation_bar <- function(cor_df, group_name) {
  plot_df <- cor_df %>%
    filter(is.finite(cor)) %>%
    mutate(
      direction = ifelse(cor >= 0, "Positive", "Negative"),
      gene = factor(gene, levels = rev(gene))
    )
  
  p <- ggplot(plot_df, aes(x = cor, y = gene, fill = direction)) +
    geom_col(width = 0.76, colour = "white", linewidth = 0.5) +
    geom_vline(xintercept = 0, linewidth = 1, colour = "black") +
    scale_fill_manual(values = c("Positive" = "#E64B35FF", "Negative" = "#4DBBD5FF")) +
    labs(
      title = paste0("Correlation of FOLH1 with known biomarkers (", group_name, ", ", plot_label, ")"),
      x = "Spearman correlation with FOLH1",
      y = NULL
    ) +
    theme_classic(base_size = 18, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold", colour = "black", size = 18),
      axis.text = element_text(colour = "black", size = 15),
      axis.line = element_line(linewidth = 0.9, colour = "black"),
      axis.ticks = element_line(linewidth = 0.8, colour = "black"),
      legend.position = "none"
    )
  
  ggsave(
    file.path(corr_bar_dir, paste0("FOLH1_biomarker_correlation_barplot_", group_name, "_MTK_treated.png")),
    p,
    width = 10,
    height = 8,
    dpi = 300
  )
}

plot_correlation_bar(cor_results_control, "control")
plot_correlation_bar(cor_results_treated, "treated")

scatter_fun_group <- function(g, group_name, save_dir) {
  samples_in_group <- rownames(sample_info)[sample_info$condition == group_name]
  
  folh1_rows <- get_symbol_rows("FOLH1")
  folh1_expr_group <- as.numeric(colMeans(vst_mat[folh1_rows, samples_in_group, drop = FALSE], na.rm = TRUE))
  
  rows <- get_symbol_rows(g)
  biomarker_expr <- as.numeric(colMeans(vst_mat[rows, samples_in_group, drop = FALSE], na.rm = TRUE))
  
  df <- data.frame(
    sample = samples_in_group,
    sample_short = sample_info[samples_in_group, "sample_short"],
    sample_short10 = sample_info[samples_in_group, "sample_short10"],
    FOLH1 = folh1_expr_group,
    biomarker = biomarker_expr,
    group = group_name,
    stringsAsFactors = FALSE
  )
  
  ok <- is.finite(df$FOLH1) & is.finite(df$biomarker)
  df <- df[ok, , drop = FALSE]
  
  if (nrow(df) < 3) return(NULL)
  
  ct <- suppressWarnings(cor.test(df$FOLH1, df$biomarker, method = "spearman", exact = FALSE))
  lab <- paste0("Spearman rho = ", round(unname(ct$estimate), 3), ", p = ", signif(ct$p.value, 3))
  
  point_col <- ifelse(group_name == "control", "#F8766D", "#00BFC4")
  
  p <- ggplot(df, aes(x = FOLH1, y = biomarker)) +
    geom_point(size = 2.8, alpha = 0.9, color = point_col) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.85) +
    labs(
      title = paste0("FOLH1 vs ", g, " (", group_name, ", ", plot_label, ")"),
      subtitle = lab,
      x = "FOLH1 (VST)",
      y = paste0(g, " (VST)")
    ) +
    theme_nature(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave(
    file.path(save_dir, paste0("FOLH1_vs_", g, "_scatter_", group_name, "_MTK_treated.png")),
    p,
    width = 6.2,
    height = 5.2,
    dpi = 300
  )
  
  invisible(df)
}

for (g in available_biomarkers[available_biomarkers != "FOLH1"]) {
  scatter_fun_group(g, "control", scatter_control_dir)
  scatter_fun_group(g, "treated", scatter_treated_dir)
}

# ===================
# BIOMARKER HEATMAP
# ===================

heat_genes <- available_biomarkers
heat_rows <- unlist(lapply(heat_genes, get_symbol_rows))
heat_rows <- unique(heat_rows)
heat_mat <- vst_mat[heat_rows, , drop = FALSE]

heat_anno_col <- data.frame(condition = sample_info$condition)
rownames(heat_anno_col) <- rownames(sample_info)

rownames(heat_mat) <- row_map$SYMBOL[match(rownames(heat_mat), row_map$ENSEMBL)]
rownames(heat_mat)[is.na(rownames(heat_mat)) | rownames(heat_mat) == ""] <- heat_rows[is.na(rownames(heat_mat)) | rownames(heat_mat) == ""]

pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = heat_anno_col,
  show_colnames = TRUE,
  labels_col = sample_info[colnames(heat_mat), "sample_short10"],
  main = paste0("Mouse biomarker panel and FOLH1 (", plot_label, ")"),
  fontsize_col = 8,
  filename = file.path(outdir, "biomarker_heatmap_MTK_treated.png"),
  width = 10,
  height = 8
)

# =============
# DEG ANALYSIS
# =============

res <- results(
  dds,
  contrast = c("condition", "treated", "control"),
  alpha = 0.05
)

res_df <- as.data.frame(res) %>%
  rownames_to_column("ENSEMBL") %>%
  left_join(
    row_map %>% select(ENSEMBL, SYMBOL, ENTREZID, GENENAME),
    by = "ENSEMBL"
  ) %>%
  arrange(padj, desc(abs(log2FoldChange)))

res_df$regulated <- case_when(
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange > 0.58 ~ "Up_in_treated",
  
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    res_df$log2FoldChange < -0.58 ~ "Down_in_treated",
  
  TRUE ~ "Not_significant"
)


deg_sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0.58)

deg_top50 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 50)

deg_top100 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 100)

deg_summary <- data.frame(
  Comparison = "treated_vs_control",
  Total_DE_Genes = nrow(deg_sig),
  Upregulated = sum(deg_sig$log2FoldChange > 0.58, na.rm = TRUE),
  Downregulated = sum(deg_sig$log2FoldChange < -0.58, na.rm = TRUE),
  padj_Cutoff = 0.05,
  log2FC_Cutoff = 0.58,
  stringsAsFactors = FALSE
)

write.csv(res_df, file.path(deg_table_dir, "DESeq2_all_results_treated_vs_control_MTK_treated.csv"), row.names = FALSE)
write.csv(deg_sig, file.path(deg_table_dir, "DESeq2_significant_DEGs_treated_vs_control_MTK_treated.csv"), row.names = FALSE)
write.csv(deg_top50, file.path(deg_table_dir, "DESeq2_top50_ranked_genes_treated_vs_control_MTK_treated.csv"), row.names = FALSE)
write.csv(deg_top100, file.path(deg_table_dir, "DESeq2_top100_ranked_genes_treated_vs_control_MTK_treated.csv"), row.names = FALSE)
write.csv(deg_summary, file.path(deg_table_dir, "DEG_summary_treated_vs_control_MTK_treated.csv"), row.names = FALSE)

plot_deg_heatmap <- function(gene_table, mat, sample_info, filename, title_text, scale_rows = TRUE) {
  genes_use <- intersect(gene_table$ENSEMBL, rownames(mat))
  if (length(genes_use) < 2) return(NULL)
  
  hm_mat <- mat[genes_use, , drop = FALSE]
  gene_symbols <- gene_table$SYMBOL[match(rownames(hm_mat), gene_table$ENSEMBL)]
  gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- rownames(hm_mat)
  rownames(hm_mat) <- make.unique(gene_symbols)
  
  ann_col <- data.frame(
    condition = sample_info[colnames(hm_mat), "condition"],
    row.names = colnames(hm_mat)
  )
  
  pheatmap(
    hm_mat,
    scale = if (scale_rows) "row" else "none",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = ann_col,
    labels_col = sample_info[colnames(hm_mat), "sample_short10"],
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = ifelse(nrow(hm_mat) > 80, 5, 7),
    fontsize_col = 8,
    main = paste0(title_text, " (", plot_label, ")"),
    filename = filename,
    width = 10,
    height = 12
  )
}

plot_deg_heatmap(
  gene_table = deg_sig,
  mat = vst_mat,
  sample_info = sample_info,
  filename = file.path(deg_heatmap_dir, "Heatmap_all_significant_DEGs_treated_vs_control_MTK_treated.png"),
  title_text = "Heatmap of all significant DEGs"
)

plot_deg_heatmap(
  gene_table = deg_top50,
  mat = vst_mat,
  sample_info = sample_info,
  filename = file.path(deg_heatmap_dir, "Heatmap_top50_DEGs_treated_vs_control_MTK_treated.png"),
  title_text = "Heatmap of top 50 DEGs"
)

plot_deg_heatmap(
  gene_table = deg_top100,
  mat = vst_mat,
  sample_info = sample_info,
  filename = file.path(deg_heatmap_dir, "Heatmap_top100_DEGs_treated_vs_control_MTK_treated.png"),
  title_text = "Heatmap of top 100 DEGs"
)

# =================
# ENRICHMENT PREP
# =================

res_ranked <- res_df %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(rank_metric = log2FoldChange) %>%
  arrange(desc(rank_metric))

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
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
} else NULL

ego_down <- if (length(sig_entrez_down) > 0) {
  enrichGO(
    gene = sig_entrez_down,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
} else NULL

if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
  write.csv(as.data.frame(ego_up), file.path(go_bp_dir, "GO_BP_ORA_up_in_treated_MTK_treated.csv"), row.names = FALSE)
  png(file.path(go_bp_dir, "GO_BP_ORA_up_in_treated_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(ego_up, showCategory = 20) +
      ggtitle(paste0("GO BP enriched in upregulated genes (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
  write.csv(as.data.frame(ego_down), file.path(go_bp_dir, "GO_BP_ORA_down_in_treated_MTK_treated.csv"), row.names = FALSE)
  png(file.path(go_bp_dir, "GO_BP_ORA_down_in_treated_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(ego_down, showCategory = 20) +
      ggtitle(paste0("GO BP enriched in downregulated genes (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

# ========
# MSIGDB
# ========

msig_hallmark_mm <- msigdbr(species = "Mus musculus", category = "H")
hallmark_term2gene <- msig_hallmark_mm %>% select(gs_name, entrez_gene)

msig_immune_mm <- msigdbr(species = "Mus musculus", category = "C7")
immune_term2gene <- msig_immune_mm %>% select(gs_name, entrez_gene)

# =======================
# HALLMARK ORA + GSEA
# =======================

hallmark_up <- if (length(sig_entrez_up) > 0) {
  enricher(
    gene = sig_entrez_up,
    TERM2GENE = hallmark_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
} else NULL

hallmark_down <- if (length(sig_entrez_down) > 0) {
  enricher(
    gene = sig_entrez_down,
    TERM2GENE = hallmark_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
} else NULL

gsea_hallmark <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (!is.null(hallmark_up) && nrow(as.data.frame(hallmark_up)) > 0) {
  write.csv(as.data.frame(hallmark_up), file.path(hallmark_dir, "Hallmark_ORA_up_in_treated_MTK_treated.csv"), row.names = FALSE)
  png(file.path(hallmark_dir, "Hallmark_ORA_up_in_treated_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(hallmark_up, showCategory = 20) +
      ggtitle(paste0("Hallmark enriched in upregulated genes (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

if (!is.null(hallmark_down) && nrow(as.data.frame(hallmark_down)) > 0) {
  write.csv(as.data.frame(hallmark_down), file.path(hallmark_dir, "Hallmark_ORA_down_in_treated_MTK_treated.csv"), row.names = FALSE)
  png(file.path(hallmark_dir, "Hallmark_ORA_down_in_treated_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(hallmark_down, showCategory = 20) +
      ggtitle(paste0("Hallmark enriched in downregulated genes (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
  write.csv(as.data.frame(gsea_hallmark), file.path(gsea_dir, "GSEA_Hallmark_ranked_MTK_treated.csv"), row.names = FALSE)
  png(file.path(gsea_dir, "GSEA_Hallmark_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(gsea_hallmark, showCategory = 20, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(paste0("Hallmark GSEA ranked analysis (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

# ========================
# IMMUNESIGDB ORA + GSEA
# ========================

immune_up <- if (length(sig_entrez_up) > 0) {
  enricher(
    gene = sig_entrez_up,
    TERM2GENE = immune_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
} else NULL

immune_down <- if (length(sig_entrez_down) > 0) {
  enricher(
    gene = sig_entrez_down,
    TERM2GENE = immune_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
} else NULL

gsea_immune <- GSEA(
  geneList = gene_list,
  TERM2GENE = immune_term2gene,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (!is.null(immune_up) && nrow(as.data.frame(immune_up)) > 0) {
  write.csv(as.data.frame(immune_up), file.path(immune_dir, "ImmuneSigDB_ORA_up_in_treated_MTK_treated.csv"), row.names = FALSE)
  png(file.path(immune_dir, "ImmuneSigDB_ORA_up_in_treated_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(immune_up, showCategory = 20) +
      ggtitle(paste0("ImmuneSigDB enriched in upregulated genes (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

if (!is.null(immune_down) && nrow(as.data.frame(immune_down)) > 0) {
  write.csv(as.data.frame(immune_down), file.path(immune_dir, "ImmuneSigDB_ORA_down_in_treated_MTK_treated.csv"), row.names = FALSE)
  png(file.path(immune_dir, "ImmuneSigDB_ORA_down_in_treated_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(immune_down, showCategory = 20) +
      ggtitle(paste0("ImmuneSigDB enriched in downregulated genes (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

if (!is.null(gsea_immune) && nrow(as.data.frame(gsea_immune)) > 0) {
  write.csv(as.data.frame(gsea_immune), file.path(gsea_dir, "GSEA_ImmuneSigDB_ranked_MTK_treated.csv"), row.names = FALSE)
  png(file.path(gsea_dir, "GSEA_ImmuneSigDB_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(gsea_immune, showCategory = 20, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(paste0("ImmuneSigDB GSEA ranked analysis (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

# ================
# GO BP GSEA
# ================

gsea_go_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (!is.null(gsea_go_bp) && nrow(as.data.frame(gsea_go_bp)) > 0) {
  write.csv(as.data.frame(gsea_go_bp), file.path(gsea_dir, "GSEA_GO_BP_ranked_MTK_treated.csv"), row.names = FALSE)
  png(file.path(gsea_dir, "GSEA_GO_BP_dotplot_MTK_treated.png"), width = 2400, height = 1800, res = 300)
  print(
    dotplot(gsea_go_bp, showCategory = 20, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(paste0("GO BP GSEA ranked analysis (", plot_label, ")")) +
      theme(axis.text.y = element_text(size = 6))
  )
  dev.off()
}

# ================================
# SIMPLE DEG SUMMARY TABLE IMAGE
# ================================

library(gridExtra)
library(grid)

deg_summary_tbl <- tableGrob(
  deg_summary,
  rows = NULL,
  theme = ttheme_minimal(base_size = 11)
)

png(file.path(deg_dir, "DEG_summary_table_MTK_treated.png"), width = 2600, height = 700, res = 300)
grid.newpage()
grid.draw(deg_summary_tbl)
dev.off()

# ===================
# EXCEL WORKBOOK
# ===================

wb <- createWorkbook()

addWorksheet(wb, "Sample_info")
writeData(wb, "Sample_info", sample_info)

addWorksheet(wb, "PCA_QC")
writeData(wb, "PCA_QC", pca_df)

addWorksheet(wb, "Biomarker_presence")
writeData(wb, "Biomarker_presence", biomarker_presence)

addWorksheet(wb, "Wilcox_results")
writeData(wb, "Wilcox_results", group_test_results)

addWorksheet(wb, "Group_summary")
writeData(wb, "Group_summary", group_summary)

addWorksheet(wb, "Expression_long")
writeData(wb, "Expression_long", gene_expr_df)

addWorksheet(wb, "Spearman_control")
writeData(wb, "Spearman_control", cor_results_control)

addWorksheet(wb, "Spearman_treated")
writeData(wb, "Spearman_treated", cor_results_treated)

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

if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
  addWorksheet(wb, "GO_BP_up")
  writeData(wb, "GO_BP_up", as.data.frame(ego_up))
}
if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
  addWorksheet(wb, "GO_BP_down")
  writeData(wb, "GO_BP_down", as.data.frame(ego_down))
}
if (!is.null(hallmark_up) && nrow(as.data.frame(hallmark_up)) > 0) {
  addWorksheet(wb, "Hallmark_up")
  writeData(wb, "Hallmark_up", as.data.frame(hallmark_up))
}
if (!is.null(hallmark_down) && nrow(as.data.frame(hallmark_down)) > 0) {
  addWorksheet(wb, "Hallmark_down")
  writeData(wb, "Hallmark_down", as.data.frame(hallmark_down))
}
if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
  addWorksheet(wb, "GSEA_Hallmark")
  writeData(wb, "GSEA_Hallmark", as.data.frame(gsea_hallmark))
}
if (!is.null(gsea_go_bp) && nrow(as.data.frame(gsea_go_bp)) > 0) {
  addWorksheet(wb, "GSEA_GO_BP")
  writeData(wb, "GSEA_GO_BP", as.data.frame(gsea_go_bp))
}
if (!is.null(immune_up) && nrow(as.data.frame(immune_up)) > 0) {
  addWorksheet(wb, "ImmuneSigDB_up")
  writeData(wb, "ImmuneSigDB_up", as.data.frame(immune_up))
}
if (!is.null(immune_down) && nrow(as.data.frame(immune_down)) > 0) {
  addWorksheet(wb, "ImmuneSigDB_down")
  writeData(wb, "ImmuneSigDB_down", as.data.frame(immune_down))
}
if (!is.null(gsea_immune) && nrow(as.data.frame(gsea_immune)) > 0) {
  addWorksheet(wb, "GSEA_ImmuneSigDB")
  writeData(wb, "GSEA_ImmuneSigDB", as.data.frame(gsea_immune))
}

saveWorkbook(
  wb,
  file = file.path(outdir, "Overall_mouse_FOLH1_biomarker_and_DEG_results_MTK_treated.xlsx"),
  overwrite = TRUE
)