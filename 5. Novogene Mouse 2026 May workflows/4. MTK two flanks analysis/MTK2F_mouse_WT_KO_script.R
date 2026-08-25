# ====================================================================
# Ravindu: Novogene 2026 Mice Analysis
# MTK Treatment (Two Flanks) 
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
  library(gridExtra)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(rstatix)
})

# OUTPUT DIRECTORIES

outdir <- "mouse_MTK_two_flanks_biomarker_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

qc_dir <- file.path(outdir, "01_QC")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

deg_root_dir <- file.path(outdir, "02_DEG_analysis")
dir.create(deg_root_dir, showWarnings = FALSE, recursive = TRUE)

expr_plot_dir <- file.path(outdir, "03_Expression_boxplots")
dir.create(expr_plot_dir, showWarnings = FALSE, recursive = TRUE)

scatter_root_dir <- file.path(outdir, "04_Correlation_scatterplots")
dir.create(scatter_root_dir, showWarnings = FALSE, recursive = TRUE)

scatter_wt_control_dir <- file.path(scatter_root_dir, "WT_control")
scatter_wt_mtk_dir <- file.path(scatter_root_dir, "WT_MTK")
scatter_ko_control_dir <- file.path(scatter_root_dir, "KO_control")
scatter_ko_mtk_dir <- file.path(scatter_root_dir, "KO_MTK")

dir.create(scatter_wt_control_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(scatter_wt_mtk_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(scatter_ko_control_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(scatter_ko_mtk_dir, showWarnings = FALSE, recursive = TRUE)

corr_bar_dir <- file.path(outdir, "05_Correlation_barplots")
dir.create(corr_bar_dir, showWarnings = FALSE, recursive = TRUE)

enrich_root_dir <- file.path(outdir, "06_Pathway_analysis")
dir.create(enrich_root_dir, showWarnings = FALSE, recursive = TRUE)

plot_label <- "MTK two flanks"

# THEME

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

group_colors <- c(
  "WT_control" = "#F8766D",
  "WT_MTK" = "#7CAE00",
  "KO_control" = "#00BFC4",
  "KO_MTK" = "#C77CFF"
)

# =========
# INPUT
# =========

count_file <- file.path("../counts", "Galaxy Column join on MTK two flanks.tabular")
raw <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gene_id_raw <- as.character(raw[[1]])
gene_id <- sub("\\..*$", "", gene_id_raw)

count_df <- raw[, -1, drop = FALSE]
count_mat <- as.matrix(count_df)
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- gene_id

sample_names <- colnames(count_mat)
sample_short7 <- substr(sample_names, 1, 7)
sample_short10 <- substr(sample_names, 1, 10)

wtko <- ifelse(
  grepl("WT", sample_names, ignore.case = TRUE), "WT",
  ifelse(grepl("KO", sample_names, ignore.case = TRUE), "KO", NA)
)

cond <- ifelse(
  substr(sample_names, 7, 7) == "C", "control",
  ifelse(substr(sample_names, 7, 7) == "M", "MTK", NA)
)

group <- factor(
  paste(wtko, cond, sep = "_"),
  levels = c("WT_control", "WT_MTK", "KO_control", "KO_MTK")
)

sample_info <- data.frame(
  sample = sample_names,
  sample_short = sample_short7,
  sample_short10 = sample_short10,
  wtko = wtko,
  condition = cond,
  group = group,
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_names

# =================
# FILTER / DESEQ2
# =================

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
  file = file.path(outdir, "normalized_counts_MTK_two_flanks.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = rownames(vst_mat), vst_mat, check.names = FALSE),
  file = file.path(outdir, "vst_expression_matrix_MTK_two_flanks.csv"),
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
  wtko = sample_info[rownames(pca_obj$x), "wtko"],
  condition = sample_info[rownames(pca_obj$x), "condition"],
  group = sample_info[rownames(pca_obj$x), "group"],
  stringsAsFactors = FALSE
)

write.csv(
  pca_df,
  file = file.path(qc_dir, "PCA_coordinates_PC1_PC2_MTK_two_flanks.csv"),
  row.names = FALSE
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
  geom_point(size = 4, alpha = 0.95) +
  stat_ellipse(
    aes(group = group),
    type = "norm",
    geom = "polygon",
    alpha = 0.12,
    colour = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
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
  file.path(qc_dir, "PCA_PC1_vs_PC2_MTK_two_flanks.png"),
  p_pca,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

sample_cor_mat <- cor(vst_mat, method = "spearman", use = "pairwise.complete.obs")

pheatmap(
  sample_cor_mat,
  annotation_col = data.frame(group = sample_info$group, row.names = rownames(sample_info)),
  annotation_row = data.frame(group = sample_info$group, row.names = rownames(sample_info)),
  labels_col = sample_info[colnames(sample_cor_mat), "sample_short10"],
  labels_row = sample_info[rownames(sample_cor_mat), "sample_short10"],
  main = paste0("Sample-to-sample Spearman correlation (", plot_label, ")"),
  fontsize_col = 8,
  fontsize_row = 8,
  filename = file.path(qc_dir, "Sample_correlation_heatmap_MTK_two_flanks.png"),
  width = 8,
  height = 7
)

hc <- hclust(dist(t(vst_mat), method = "euclidean"), method = "complete")
png(file.path(qc_dir, "Sample_clustering_dendrogram_MTK_two_flanks.png"), width = 2600, height = 2000, res = 300)
plot(
  hc,
  labels = sample_info[hc$labels, "sample_short10"],
  main = paste0("Sample clustering (vst, Euclidean) (", plot_label, ")"),
  xlab = "",
  sub = "",
  cex = 0.8
)
dev.off()

# ===========
# ANNOTATION
# ===========

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
  file = file.path(outdir, "biomarker_presence_check_MTK_two_flanks.csv"),
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

# =================
# EXPRESSION TABLE
# =================

gene_expr_df <- lapply(available_biomarkers, function(g) {
  rows <- get_symbol_rows(g)
  x <- as.numeric(colMeans(vst_mat[rows, , drop = FALSE], na.rm = TRUE))
  data.frame(
    sample = colnames(vst_mat),
    sample_short = sample_info[colnames(vst_mat), "sample_short"],
    sample_short10 = sample_info[colnames(vst_mat), "sample_short10"],
    wtko = sample_info[colnames(vst_mat), "wtko"],
    condition = sample_info[colnames(vst_mat), "condition"],
    group = sample_info[colnames(vst_mat), "group"],
    gene = g,
    expr = x,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(
  gene_expr_df,
  file = file.path(outdir, "biomarker_expression_long_VST_MTK_two_flanks.csv"),
  row.names = FALSE
)

# ==========================
# BOXPLOTS + KRUSKAL + DUNN
# ==========================

plot_group_box <- function(df, gene, ylab = "VST expression") {
  kruskal_res <- tryCatch(kruskal.test(expr ~ group, data = df), error = function(e) NULL)
  p_label <- if (!is.null(kruskal_res)) paste0("p = ", signif(kruskal_res$p.value, 3)) else "Kruskal-Wallis p = NA"
  
  y_range <- range(df$expr, na.rm = TRUE)
  y_diff <- diff(y_range)
  if (y_diff == 0) y_diff <- 0.5
  y_pad <- 0.18 * y_diff
  
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.9, linewidth = 0.9) +
    geom_jitter(width = 0.12, size = 2.0, alpha = 0.55, color = "black") +
    annotate("text", x = 2.5, y = y_range[2] + 0.45 * y_pad, label = p_label, size = 4.3, fontface = "bold") +
    scale_fill_manual(values = group_colors) +
    labs(
      title = paste0(gene, " (", plot_label, ")"),
      subtitle = "Kruskal-Wallis test across 4 groups",
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
    file.path(expr_plot_dir, paste0(gene, "_group_boxplot_MTK_two_flanks.png")),
    p,
    width = 7.4,
    height = 5.9,
    dpi = 300
  )
  
  important_comparisons <- tibble::tribble(
    ~group1,        ~group2,
    "WT_control",   "WT_MTK",
    "KO_control",   "KO_MTK",
    "WT_control",   "KO_control",
    "WT_MTK",       "KO_MTK"
  )
  
  dunn_res <- tryCatch(
    {
      rstatix::dunn_test(
        data = df,
        formula = expr ~ group,
        p.adjust.method = "BH",
        detailed = TRUE
      ) %>%
        dplyr::semi_join(important_comparisons, by = c("group1", "group2")) %>%
        dplyr::mutate(
          gene = gene,
          comparison = paste(group1, "vs", group2)
        ) %>%
        dplyr::select(
          gene, comparison, group1, group2, estimate, statistic, p, p.adj,
          p.signif, p.adj.signif, method
        )
    },
    error = function(e) {
      data.frame(
        gene = gene,
        comparison = NA_character_,
        group1 = NA_character_,
        group2 = NA_character_,
        estimate = NA_real_,
        statistic = NA_real_,
        p = NA_real_,
        p.adj = NA_real_,
        p.signif = NA_character_,
        p.adj.signif = NA_character_,
        method = "Dunn test",
        stringsAsFactors = FALSE
      )
    }
  )
  
  list(
    kruskal = data.frame(
      gene = gene,
      test = "Kruskal-Wallis",
      p.value = if (!is.null(kruskal_res)) kruskal_res$p.value else NA_real_,
      stringsAsFactors = FALSE
    ),
    dunn = dunn_res
  )
}

boxplot_results <- lapply(available_biomarkers, function(g) {
  df_g <- gene_expr_df %>% dplyr::filter(gene == g)
  plot_group_box(df_g, g)
})

group_test_results <- dplyr::bind_rows(lapply(boxplot_results, `[[`, "kruskal"))
dunn_posthoc_results <- dplyr::bind_rows(lapply(boxplot_results, `[[`, "dunn"))

write.csv(
  group_test_results,
  file = file.path(outdir, "group_kruskal_biomarkers_MTK_two_flanks.csv"),
  row.names = FALSE
)

write.csv(
  dunn_posthoc_results,
  file = file.path(outdir, "group_dunn_posthoc_biomarkers_MTK_two_flanks.csv"),
  row.names = FALSE
)

group_summary <- gene_expr_df %>%
  group_by(gene, group, wtko, condition) %>%
  summarise(
    n = sum(is.finite(expr)),
    mean_expr = mean(expr, na.rm = TRUE),
    median_expr = median(expr, na.rm = TRUE),
    sd_expr = sd(expr, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  group_summary,
  file = file.path(outdir, "biomarker_group_summary_MTK_two_flanks.csv"),
  row.names = FALSE
)

# ==============
# CORRELATIONS
# ==============

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

cor_results_wt_control <- calc_group_correlations("WT_control")
cor_results_wt_mtk <- calc_group_correlations("WT_MTK")
cor_results_ko_control <- calc_group_correlations("KO_control")
cor_results_ko_mtk <- calc_group_correlations("KO_MTK")

cor_results_all <- bind_rows(
  cor_results_wt_control,
  cor_results_wt_mtk,
  cor_results_ko_control,
  cor_results_ko_mtk
)

write.csv(cor_results_wt_control, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_WT_control_MTK_two_flanks.csv"), row.names = FALSE)
write.csv(cor_results_wt_mtk, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_WT_MTK_MTK_two_flanks.csv"), row.names = FALSE)
write.csv(cor_results_ko_control, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_KO_control_MTK_two_flanks.csv"), row.names = FALSE)
write.csv(cor_results_ko_mtk, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_KO_MTK_MTK_two_flanks.csv"), row.names = FALSE)
write.csv(cor_results_all, file = file.path(outdir, "FOLH1_vs_biomarkers_spearman_by_group_MTK_two_flanks.csv"), row.names = FALSE)

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
    file.path(corr_bar_dir, paste0("FOLH1_biomarker_correlation_barplot_", group_name, "_MTK_two_flanks.png")),
    p,
    width = 10,
    height = 8,
    dpi = 300
  )
}

plot_correlation_bar(cor_results_wt_control, "WT_control")
plot_correlation_bar(cor_results_wt_mtk, "WT_MTK")
plot_correlation_bar(cor_results_ko_control, "KO_control")
plot_correlation_bar(cor_results_ko_mtk, "KO_MTK")

# ===============
# SCATTERPLOTS
# ===============

scatter_fun_group <- function(g, group_name, save_dir) {
  samples_in_group <- rownames(sample_info)[sample_info$group == group_name]
  
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
  
  point_col <- switch(
    group_name,
    "WT_control" = "#F8766D",
    "WT_MTK" = "#7CAE00",
    "KO_control" = "#00BFC4",
    "KO_MTK" = "#C77CFF"
  )
  
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
    file.path(save_dir, paste0("FOLH1_vs_", g, "_scatter_", group_name, "_MTK_two_flanks.png")),
    p,
    width = 6.2,
    height = 5.2,
    dpi = 300
  )
  
  invisible(df)
}

for (g in available_biomarkers[available_biomarkers != "FOLH1"]) {
  scatter_fun_group(g, "WT_control", scatter_wt_control_dir)
  scatter_fun_group(g, "WT_MTK", scatter_wt_mtk_dir)
  scatter_fun_group(g, "KO_control", scatter_ko_control_dir)
  scatter_fun_group(g, "KO_MTK", scatter_ko_mtk_dir)
}

# ==================
# BIOMARKER HEATMAP
# ==================

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
  show_colnames = TRUE,
  labels_col = sample_info[colnames(heat_mat), "sample_short10"],
  main = paste0("Mouse biomarker panel and FOLH1 (", plot_label, ")"),
  fontsize_col = 8,
  filename = file.path(outdir, "biomarker_heatmap_MTK_two_flanks.png"),
  width = 10,
  height = 8
)

# =================
# HELPER FUNCTIONS
# =================

plot_deg_heatmap <- function(gene_table, mat, sample_info, filename, title_text, scale_rows = TRUE) {
  genes_use <- intersect(gene_table$ENSEMBL, rownames(mat))
  if (length(genes_use) < 2) return(NULL)
  
  hm_mat <- mat[genes_use, , drop = FALSE]
  gene_symbols <- gene_table$SYMBOL[match(rownames(hm_mat), gene_table$ENSEMBL)]
  gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- rownames(hm_mat)
  rownames(hm_mat) <- make.unique(gene_symbols)
  
  ann_col <- data.frame(
    group = sample_info[colnames(hm_mat), "group"],
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
    fontsize_row = ifelse(nrow(hm_mat) > 80, 3.5, 6),
    fontsize_col = 8,
    main = paste0(title_text, " (", plot_label, ")"),
    filename = filename,
    width = 10,
    height = 12
  )
}

write_enrichment_plot <- function(enrich_obj, file_csv, file_png, plot_title, width = 9, height = 6) {
  if (!is.null(enrich_obj) && nrow(as.data.frame(enrich_obj)) > 0) {
    write.csv(as.data.frame(enrich_obj), file_csv, row.names = FALSE)
    p <- dotplot(enrich_obj, showCategory = 20) +
      ggtitle(plot_title) +
      theme(axis.text.y = element_text(size = 5))
    ggsave(file_png, p, width = width, height = height, dpi = 300)
  }
}

write_gsea_plot <- function(gsea_obj, file_csv, file_png, plot_title, width = 10, height = 6) {
  if (!is.null(gsea_obj) && nrow(as.data.frame(gsea_obj)) > 0) {
    write.csv(as.data.frame(gsea_obj), file_csv, row.names = FALSE)
    p <- dotplot(gsea_obj, showCategory = 20, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(plot_title) +
      theme(axis.text.y = element_text(size = 5))
    ggsave(file_png, p, width = width, height = height, dpi = 300)
  }
}

run_deg_and_pathway_analysis <- function(dds, vst_mat, sample_info, row_map, comparison_name, numerator, denominator,
                                         out_deg_root, out_enrich_root, plot_label,
                                         lfc_cutoff = 0.58, padj_cutoff = 0.05) {
  
  comp_deg_dir <- file.path(out_deg_root, comparison_name)
  dir.create(comp_deg_dir, showWarnings = FALSE, recursive = TRUE)
  
  deg_table_dir <- file.path(comp_deg_dir, "01_DEG_tables")
  deg_heatmap_dir <- file.path(comp_deg_dir, "02_DEG_heatmaps")
  dir.create(deg_table_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(deg_heatmap_dir, showWarnings = FALSE, recursive = TRUE)
  
  comp_enrich_dir <- file.path(out_enrich_root, comparison_name)
  dir.create(comp_enrich_dir, showWarnings = FALSE, recursive = TRUE)
  
  go_bp_dir <- file.path(comp_enrich_dir, "01_GO_BP")
  hallmark_dir <- file.path(comp_enrich_dir, "02_Hallmark")
  gsea_dir <- file.path(comp_enrich_dir, "03_GSEA_ranked")
  immune_dir <- file.path(comp_enrich_dir, "04_ImmuneSigDB")
  
  dir.create(go_bp_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(hallmark_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(immune_dir, showWarnings = FALSE, recursive = TRUE)
  
  res <- results(dds, contrast = c("group", numerator, denominator), alpha = padj_cutoff)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(row_map %>% select(ENSEMBL, SYMBOL, ENTREZID, GENENAME), by = "ENSEMBL") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  res_df$regulated <- case_when(
    !is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange > lfc_cutoff ~ paste0("Up_in_", numerator),
    !is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff ~ paste0("Down_in_", numerator),
    TRUE ~ "Not_significant"
  )
  
  deg_sig <- res_df %>%
    filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  
  deg_top50 <- res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 50)
  
  deg_top100 <- res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 100)
  
  deg_summary <- data.frame(
    Comparison = comparison_name,
    Numerator = numerator,
    Denominator = denominator,
    Total_DE_Genes = nrow(deg_sig),
    Upregulated = sum(deg_sig$log2FoldChange > lfc_cutoff, na.rm = TRUE),
    Downregulated = sum(deg_sig$log2FoldChange < -lfc_cutoff, na.rm = TRUE),
    padj_Cutoff = padj_cutoff,
    log2FC_Cutoff = lfc_cutoff,
    stringsAsFactors = FALSE
  )
  
  write.csv(res_df, file.path(deg_table_dir, paste0("DESeq2_all_results_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(deg_sig, file.path(deg_table_dir, paste0("DESeq2_significant_DEGs_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(deg_top50, file.path(deg_table_dir, paste0("DESeq2_top50_ranked_genes_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(deg_top100, file.path(deg_table_dir, paste0("DESeq2_top100_ranked_genes_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(deg_summary, file.path(deg_table_dir, paste0("DEG_summary_", comparison_name, ".csv")), row.names = FALSE)
  
  contrast_samples <- rownames(sample_info)[sample_info$group %in% c(numerator, denominator)]
  hm_sample_info <- sample_info[contrast_samples, , drop = FALSE]
  hm_mat_use <- vst_mat[, contrast_samples, drop = FALSE]
  
  plot_deg_heatmap(
    gene_table = deg_sig,
    mat = hm_mat_use,
    sample_info = hm_sample_info,
    filename = file.path(deg_heatmap_dir, paste0("Heatmap_all_significant_DEGs_", comparison_name, ".png")),
    title_text = paste0("Heatmap of all significant DEGs: ", comparison_name)
  )
  
  plot_deg_heatmap(
    gene_table = deg_top50,
    mat = hm_mat_use,
    sample_info = hm_sample_info,
    filename = file.path(deg_heatmap_dir, paste0("Heatmap_top50_DEGs_", comparison_name, ".png")),
    title_text = paste0("Heatmap of top 50 DEGs: ", comparison_name)
  )
  
  plot_deg_heatmap(
    gene_table = deg_top100,
    mat = hm_mat_use,
    sample_info = hm_sample_info,
    filename = file.path(deg_heatmap_dir, paste0("Heatmap_top100_DEGs_", comparison_name, ".png")),
    title_text = paste0("Heatmap of top 100 DEGs: ", comparison_name)
  )
  
  rank_df <- res_df %>%
    filter(!is.na(ENTREZID), !is.na(log2FoldChange)) %>%
    group_by(ENTREZID) %>%
    summarise(rank_metric = log2FoldChange[which.max(abs(log2FoldChange))], .groups = "drop") %>%
    arrange(desc(rank_metric))
  
  gene_list <- rank_df$rank_metric
  names(gene_list) <- rank_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  sig_entrez_up <- deg_sig %>%
    filter(log2FoldChange > lfc_cutoff, !is.na(ENTREZID)) %>%
    pull(ENTREZID) %>%
    unique()
  
  sig_entrez_down <- deg_sig %>%
    filter(log2FoldChange < -lfc_cutoff, !is.na(ENTREZID)) %>%
    pull(ENTREZID) %>%
    unique()
  
  msig_hallmark_mm <- msigdbr( species = "Mus musculus", category = "H")
  hallmark_term2gene <- msig_hallmark_mm %>% select(gs_name, entrez_gene)
  
  msig_immune_mm <- msigdbr( species = "Mus musculus", category = "C7")
  immune_term2gene <- msig_immune_mm %>% select(gs_name, entrez_gene)
  
  ego_up <- if (length(sig_entrez_up) > 0) {
    enrichGO(
      gene = sig_entrez_up,
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
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  } else NULL
  
  write_enrichment_plot(
    ego_up,
    file.path(go_bp_dir, paste0("GO_BP_ORA_up_in_", numerator, "_", comparison_name, ".csv")),
    file.path(go_bp_dir, paste0("GO_BP_ORA_up_in_", numerator, "_dotplot_", comparison_name, ".png")),
    paste0("GO BP enriched in upregulated genes: ", comparison_name, " (", plot_label, ")")
  )
  
  write_enrichment_plot(
    ego_down,
    file.path(go_bp_dir, paste0("GO_BP_ORA_down_in_", numerator, "_", comparison_name, ".csv")),
    file.path(go_bp_dir, paste0("GO_BP_ORA_down_in_", numerator, "_dotplot_", comparison_name, ".png")),
    paste0("GO BP enriched in downregulated genes: ", comparison_name, " (", plot_label, ")")
  )
  
  hallmark_up <- if (length(sig_entrez_up) > 0) {
    enricher(
      gene = sig_entrez_up,
      TERM2GENE = hallmark_term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  } else NULL
  
  hallmark_down <- if (length(sig_entrez_down) > 0) {
    enricher(
      gene = sig_entrez_down,
      TERM2GENE = hallmark_term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  } else NULL
  
  gsea_hallmark <- GSEA(
    geneList = gene_list,
    TERM2GENE = hallmark_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
  write_enrichment_plot(
    hallmark_up,
    file.path(hallmark_dir, paste0("Hallmark_ORA_up_in_", numerator, "_", comparison_name, ".csv")),
    file.path(hallmark_dir, paste0("Hallmark_ORA_up_in_", numerator, "_dotplot_", comparison_name, ".png")),
    paste0("Hallmark enriched in upregulated genes: ", comparison_name, " (", plot_label, ")")
  )
  
  write_enrichment_plot(
    hallmark_down,
    file.path(hallmark_dir, paste0("Hallmark_ORA_down_in_", numerator, "_", comparison_name, ".csv")),
    file.path(hallmark_dir, paste0("Hallmark_ORA_down_in_", numerator, "_dotplot_", comparison_name, ".png")),
    paste0("Hallmark enriched in downregulated genes: ", comparison_name, " (", plot_label, ")")
  )
  
  write_gsea_plot(
    gsea_hallmark,
    file.path(gsea_dir, paste0("GSEA_Hallmark_ranked_", comparison_name, ".csv")),
    file.path(gsea_dir, paste0("GSEA_Hallmark_dotplot_", comparison_name, ".png")),
    paste0("Hallmark GSEA ranked analysis: ", comparison_name, " (", plot_label, ")")
  )
  
  gsea_go_bp <- gseGO(
    geneList = gene_list,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
  write_gsea_plot(
    gsea_go_bp,
    file.path(gsea_dir, paste0("GSEA_GO_BP_ranked_", comparison_name, ".csv")),
    file.path(gsea_dir, paste0("GSEA_GO_BP_dotplot_", comparison_name, ".png")),
    paste0("GO BP GSEA ranked analysis: ", comparison_name, " (", plot_label, ")")
  )
  
  immune_up <- if (length(sig_entrez_up) > 0) {
    enricher(
      gene = sig_entrez_up,
      TERM2GENE = immune_term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  } else NULL
  
  immune_down <- if (length(sig_entrez_down) > 0) {
    enricher(
      gene = sig_entrez_down,
      TERM2GENE = immune_term2gene,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  } else NULL
  
  gsea_immune <- GSEA(
    geneList = gene_list,
    TERM2GENE = immune_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    verbose = FALSE
  )
  
  write_enrichment_plot(
    immune_up,
    file.path(immune_dir, paste0("ImmuneSigDB_ORA_up_in_", numerator, "_", comparison_name, ".csv")),
    file.path(immune_dir, paste0("ImmuneSigDB_ORA_up_in_", numerator, "_dotplot_", comparison_name, ".png")),
    paste0("ImmuneSigDB enriched in upregulated genes: ", comparison_name, " (", plot_label, ")")
  )
  
  write_enrichment_plot(
    immune_down,
    file.path(immune_dir, paste0("ImmuneSigDB_ORA_down_in_", numerator, "_", comparison_name, ".csv")),
    file.path(immune_dir, paste0("ImmuneSigDB_ORA_down_in_", numerator, "_dotplot_", comparison_name, ".png")),
    paste0("ImmuneSigDB enriched in downregulated genes: ", comparison_name, " (", plot_label, ")")
  )
  
  write_gsea_plot(
    gsea_immune,
    file.path(gsea_dir, paste0("GSEA_ImmuneSigDB_ranked_", comparison_name, ".csv")),
    file.path(gsea_dir, paste0("GSEA_ImmuneSigDB_dotplot_", comparison_name, ".png")),
    paste0("ImmuneSigDB GSEA ranked analysis: ", comparison_name, " (", plot_label, ")")
  )
  
  list(
    comparison_name = comparison_name,
    numerator = numerator,
    denominator = denominator,
    deg_summary = deg_summary,
    res_df = res_df,
    deg_sig = deg_sig,
    deg_top50 = deg_top50,
    deg_top100 = deg_top100,
    ego_up = ego_up,
    ego_down = ego_down,
    hallmark_up = hallmark_up,
    hallmark_down = hallmark_down,
    gsea_hallmark = gsea_hallmark,
    gsea_go_bp = gsea_go_bp,
    immune_up = immune_up,
    immune_down = immune_down,
    gsea_immune = gsea_immune
  )
}

# ==============================
# RUN FOUR REQUIRED COMPARISONS
# ==============================

comparison_definitions <- list(
  list(
    name = "WT_MTK_vs_WT_control",
    numerator = "WT_MTK",
    denominator = "WT_control"
  ),
  
  list(
    name = "KO_MTK_vs_KO_control",
    numerator = "KO_MTK",
    denominator = "KO_control"
  ),
  
  list(
    name = "KO_MTK_vs_WT_MTK",
    numerator = "KO_MTK",
    denominator = "WT_MTK"
  ),
  
  list(
    name = "KO_control_vs_WT_control",
    numerator = "KO_control",
    denominator = "WT_control"
  )
)

comparison_results <- lapply(comparison_definitions, function(comp) {
  run_deg_and_pathway_analysis(
    dds = dds,
    vst_mat = vst_mat,
    sample_info = sample_info,
    row_map = row_map,
    comparison_name = comp$name,
    numerator = comp$numerator,
    denominator = comp$denominator,
    out_deg_root = deg_root_dir,
    out_enrich_root = enrich_root_dir,
    plot_label = plot_label
  )
})

names(comparison_results) <- vapply(comparison_definitions, function(x) x$name, character(1))

all_deg_summaries <- bind_rows(lapply(comparison_results, `[[`, "deg_summary"))
write.csv(
  all_deg_summaries,
  file = file.path(outdir, "DEG_summary_all_comparisons_MTK_two_flanks.csv"),
  row.names = FALSE
)

# ===============
# EXCEL OUTPUT
# ===============

wb <- createWorkbook()

addWorksheet(wb, "Sample_info")
writeData(wb, "Sample_info", sample_info)

addWorksheet(wb, "Biomarker_presence")
writeData(wb, "Biomarker_presence", biomarker_presence)

addWorksheet(wb, "PCA_QC")
writeData(wb, "PCA_QC", pca_df)

addWorksheet(wb, "Group_summary")
writeData(wb, "Group_summary", group_summary)

addWorksheet(wb, "Kruskal_results")
writeData(wb, "Kruskal_results", group_test_results)

addWorksheet(wb, "Dunn_posthoc")
writeData(wb, "Dunn_posthoc", dunn_posthoc_results)

addWorksheet(wb, "Spearman_WT_control")
writeData(wb, "Spearman_WT_control", cor_results_wt_control)

addWorksheet(wb, "Spearman_WT_MTK")
writeData(wb, "Spearman_WT_MTK", cor_results_wt_mtk)

addWorksheet(wb, "Spearman_KO_control")
writeData(wb, "Spearman_KO_control", cor_results_ko_control)

addWorksheet(wb, "Spearman_KO_MTK")
writeData(wb, "Spearman_KO_MTK", cor_results_ko_mtk)

addWorksheet(wb, "Expression_long")
writeData(wb, "Expression_long", gene_expr_df)

addWorksheet(wb, "DEG_summary_all")
writeData(wb, "DEG_summary_all", all_deg_summaries)

for (nm in names(comparison_results)) {
  res_obj <- comparison_results[[nm]]
  sheet_base <- substr(gsub("[^A-Za-z0-9]", "_", nm), 1, 25)
  
  addWorksheet(wb, paste0(sheet_base, "_summary"))
  writeData(wb, paste0(sheet_base, "_summary"), res_obj$deg_summary)
  
  addWorksheet(wb, paste0(sheet_base, "_all"))
  writeData(wb, paste0(sheet_base, "_all"), res_obj$res_df)
  
  addWorksheet(wb, paste0(sheet_base, "_sig"))
  writeData(wb, paste0(sheet_base, "_sig"), res_obj$deg_sig)
  
  addWorksheet(wb, paste0(sheet_base, "_top50"))
  writeData(wb, paste0(sheet_base, "_top50"), res_obj$deg_top50)
  
  addWorksheet(wb, paste0(sheet_base, "_top100"))
  writeData(wb, paste0(sheet_base, "_top100"), res_obj$deg_top100)
  
  if (!is.null(res_obj$ego_up) && nrow(as.data.frame(res_obj$ego_up)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_GO_up"))
    writeData(wb, paste0(sheet_base, "_GO_up"), as.data.frame(res_obj$ego_up))
  }
  if (!is.null(res_obj$ego_down) && nrow(as.data.frame(res_obj$ego_down)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_GO_down"))
    writeData(wb, paste0(sheet_base, "_GO_down"), as.data.frame(res_obj$ego_down))
  }
  if (!is.null(res_obj$hallmark_up) && nrow(as.data.frame(res_obj$hallmark_up)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_Hallmark_up"))
    writeData(wb, paste0(sheet_base, "_Hallmark_up"), as.data.frame(res_obj$hallmark_up))
  }
  if (!is.null(res_obj$hallmark_down) && nrow(as.data.frame(res_obj$hallmark_down)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_Hallmark_down"))
    writeData(wb, paste0(sheet_base, "_Hallmark_down"), as.data.frame(res_obj$hallmark_down))
  }
  if (!is.null(res_obj$gsea_hallmark) && nrow(as.data.frame(res_obj$gsea_hallmark)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_GSEA_H"))
    writeData(wb, paste0(sheet_base, "_GSEA_H"), as.data.frame(res_obj$gsea_hallmark))
  }
  if (!is.null(res_obj$gsea_go_bp) && nrow(as.data.frame(res_obj$gsea_go_bp)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_GSEA_GO"))
    writeData(wb, paste0(sheet_base, "_GSEA_GO"), as.data.frame(res_obj$gsea_go_bp))
  }
  if (!is.null(res_obj$immune_up) && nrow(as.data.frame(res_obj$immune_up)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_Immune_up"))
    writeData(wb, paste0(sheet_base, "_Immune_up"), as.data.frame(res_obj$immune_up))
  }
  if (!is.null(res_obj$immune_down) && nrow(as.data.frame(res_obj$immune_down)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_Immune_down"))
    writeData(wb, paste0(sheet_base, "_Immune_down"), as.data.frame(res_obj$immune_down))
  }
  if (!is.null(res_obj$gsea_immune) && nrow(as.data.frame(res_obj$gsea_immune)) > 0) {
    addWorksheet(wb, paste0(sheet_base, "_GSEA_I"))
    writeData(wb, paste0(sheet_base, "_GSEA_I"), as.data.frame(res_obj$gsea_immune))
  }
}

saveWorkbook(
  wb,
  file = file.path(outdir, "Overall_mouse_MTK_two_flanks_biomarker_results.xlsx"),
  overwrite = TRUE
)