# ====================================================================
# Ravindu: Novogene 2026 Mice Analysis
# ISA Treatment 
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

# GLOBAL SAFETY
options(stringsAsFactors = FALSE)

# Helps in some RStudio / system setups for bitmap devices
try({
  options(bitmapType = "cairo")
}, silent = TRUE)

close_all_graphics_devices <- function() {
  while (!is.null(dev.list())) {
    try(dev.off(), silent = TRUE)
  }
  invisible(NULL)
}

safe_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

safe_png_plot <- function(file, plot_expr, width = 2400, height = 1800, res = 300) {
  safe_dir(dirname(file))
  close_all_graphics_devices()
  
  opened <- FALSE
  
  tryCatch({
    png(filename = file, width = width, height = height, res = res)
    opened <- TRUE
    eval.parent(substitute(plot_expr))
  }, error = function(e) {
    message("Skipping plot ", basename(file), ": ", conditionMessage(e))
  }, finally = {
    if (opened) {
      try(dev.off(), silent = TRUE)
    }
  })
  
  invisible(file.exists(file))
}

safe_pheatmap_file <- function(...) {
  tryCatch({
    pheatmap(...)
  }, error = function(e) {
    message("Skipping pheatmap output: ", conditionMessage(e))
  })
}

# OUTPUT DIRECTORIES

outdir <- "mouse_ISA_biomarker_analysis_sex_stratified"
safe_dir(outdir)

qc_dir <- file.path(outdir, "01_QC")
deg_root_dir <- file.path(outdir, "02_DEG_analysis")
expr_plot_dir <- file.path(outdir, "03_Expression_boxplots")
scatter_root_dir <- file.path(outdir, "04_Correlation_scatterplots")
corr_bar_dir <- file.path(outdir, "05_Correlation_barplots")
enrich_root_dir <- file.path(outdir, "06_Pathway_analysis")

safe_dir(qc_dir)
safe_dir(deg_root_dir)
safe_dir(expr_plot_dir)
safe_dir(scatter_root_dir)
safe_dir(corr_bar_dir)
safe_dir(enrich_root_dir)

plot_label <- "ISA sex-stratified"

# ==========
# THEME
# ==========

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

# =========
# INPUT
# =========

count_file <- file.path("../counts", "Galaxy Column join on ISA.tabular")

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

treatment_code <- substr(sample_names, 6, 6)
treatment <- ifelse(
  treatment_code == "C", "control",
  ifelse(treatment_code == "I", "ISA", NA)
)

female_prefixes <- c("T1729C", "T1728C", "T1730I", "T1760I")
sex <- ifelse(substr(sample_names, 1, 6) %in% female_prefixes, "female", "male")

sample_info <- data.frame(
  sample = sample_names,
  sample_short = sample_short7,
  sample_short10 = sample_short10,
  treatment = treatment,
  sex = sex,
  sex_treatment = paste(sex, treatment, sep = "_"),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_names

keep_rows <- rowSums(count_mat, na.rm = TRUE) > 1
count_mat <- count_mat[keep_rows, , drop = FALSE]

# ============================================
# GLOBAL DESEQ2 FOR NORMALIZATION / VST / QC
# ============================================
sample_info_global <- sample_info
sample_info_global$group <- factor(sample_info_global$treatment, levels = c("control", "ISA"))

dds_all <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = sample_info_global,
  design = ~ group
)

dds_all <- dds_all[rowSums(counts(dds_all)) > 1, ]
dds_all <- DESeq(dds_all)

norm_mat <- counts(dds_all, normalized = TRUE)
vst_obj <- vst(dds_all, blind = TRUE)
vst_mat <- assay(vst_obj)

write.csv(
  data.frame(gene_id = rownames(norm_mat), norm_mat, check.names = FALSE),
  file = file.path(outdir, "normalized_counts_ISA_sex_stratified.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = rownames(vst_mat), vst_mat, check.names = FALSE),
  file = file.path(outdir, "vst_expression_matrix_ISA_sex_stratified.csv"),
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
  sex = sample_info[rownames(pca_obj$x), "sex"],
  treatment = sample_info[rownames(pca_obj$x), "treatment"],
  sex_treatment = sample_info[rownames(pca_obj$x), "sex_treatment"],
  stringsAsFactors = FALSE
)

write.csv(
  pca_df,
  file = file.path(qc_dir, "PCA_coordinates_ISA_sex_stratified.csv"),
  row.names = FALSE
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment, fill = treatment, shape = sex)) +
  geom_point(size = 4, alpha = 0.95) +
  stat_ellipse(
    aes(group = treatment),
    type = "norm",
    geom = "polygon",
    alpha = 0.12,
    colour = NA,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("control" = "#F8766D", "ISA" = "#7CAE00")) +
  scale_fill_manual(values = c("control" = "#F8766D", "ISA" = "#7CAE00")) +
  geom_text(aes(label = sample_short), vjust = -0.8, size = 2.2, check_overlap = TRUE) +
  labs(
    title = paste0("PCA (", plot_label, ")"),
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_nature(base_size = 14)

ggsave(
  file.path(qc_dir, "PCA_ISA_sex_stratified.png"),
  p_pca,
  width = 7.5,
  height = 6.5,
  dpi = 300
)

sample_cor_mat <- cor(vst_mat, method = "pearson", use = "pairwise.complete.obs")

safe_pheatmap_file(
  sample_cor_mat,
  annotation_col = data.frame(
    sex = sample_info$sex,
    treatment = sample_info$treatment,
    row.names = rownames(sample_info)
  ),
  annotation_row = data.frame(
    sex = sample_info$sex,
    treatment = sample_info$treatment,
    row.names = rownames(sample_info)
  ),
  labels_col = sample_info[colnames(sample_cor_mat), "sample_short"],
  labels_row = sample_info[rownames(sample_cor_mat), "sample_short"],
  main = paste0("Sample-to-sample correlation (", plot_label, ")"),
  filename = file.path(qc_dir, "Sample_correlation_heatmap_ISA_sex_stratified.png"),
  width = 8,
  height = 7
)

hc <- hclust(dist(t(vst_mat), method = "euclidean"), method = "complete")

safe_png_plot(
  file.path(qc_dir, "Sample_clustering_dendrogram_ISA_sex_stratified.png"),
  plot(
    hc,
    labels = sample_info[hc$labels, "sample_short"],
    main = paste0("Sample clustering (VST, Euclidean) (", plot_label, ")"),
    xlab = "",
    sub = "",
    cex = 0.85
  )
)

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
# =============

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
  file = file.path(outdir, "biomarker_presence_ISA_sex_stratified.csv"),
  row.names = FALSE
)

available_biomarkers <- biomarker_presence$gene[biomarker_presence$found_in_mouse_annotation]

get_symbol_rows <- function(symbol) {
  idx <- which(row_map$SYMBOL_UPPER == toupper(symbol))
  if (length(idx) == 0) return(integer(0))
  match(row_map$ENSEMBL[idx], rownames(vst_mat))
}

if (!("FOLH1" %in% available_biomarkers)) {
  stop("FOLH1 not found in annotation for ISA set.")
}

# =======================
# EXPRESSION LONG TABLE
# =======================

gene_expr_df <- lapply(available_biomarkers, function(g) {
  rows <- get_symbol_rows(g)
  x <- as.numeric(colMeans(vst_mat[rows, , drop = FALSE], na.rm = TRUE))
  data.frame(
    sample = colnames(vst_mat),
    sample_short = sample_info[colnames(vst_mat), "sample_short"],
    sample_short10 = sample_info[colnames(vst_mat), "sample_short10"],
    sex = sample_info[colnames(vst_mat), "sex"],
    treatment = sample_info[colnames(vst_mat), "treatment"],
    sex_treatment = sample_info[colnames(vst_mat), "sex_treatment"],
    gene = g,
    expr = x,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(
  gene_expr_df,
  file = file.path(outdir, "biomarker_expression_long_ISA_sex_stratified.csv"),
  row.names = FALSE
)

# ==========
# BOXPLOTS
# ==========

plot_group_box_simple <- function(df, gene) {
  df2 <- df %>%
    mutate(group = factor(treatment, levels = c("control", "ISA")))
  
  test_res <- tryCatch(wilcox.test(expr ~ group, data = df2), error = function(e) NULL)
  p_label <- if (!is.null(test_res)) paste0("Wilcoxon p = ", signif(test_res$p.value, 3)) else "p = NA"
  
  y_range <- range(df2$expr, na.rm = TRUE)
  y_pad <- ifelse(diff(y_range) == 0, 0.5, 0.15 * diff(y_range))
  
  p <- ggplot(df2, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2, colour = "black", alpha = 0.6) +
    annotate("text", x = 1.5, y = y_range[2] + 0.4 * y_pad, label = p_label, fontface = "bold") +
    scale_fill_manual(values = c("control" = "#F8766D", "ISA" = "#7CAE00")) +
    labs(title = paste0(gene, " (ISA pooled)"), x = NULL, y = "VST expression") +
    theme_nature(base_size = 14) +
    theme(legend.position = "none")
  
  ggsave(
    file.path(expr_plot_dir, paste0(gene, "_control_vs_ISA_boxplot.png")),
    p,
    width = 6.5,
    height = 5,
    dpi = 300
  )
  
  if (!is.null(test_res)) {
    data.frame(gene = gene, p.value = test_res$p.value)
  } else {
    data.frame(gene = gene, p.value = NA_real_)
  }
}

group_test_results <- lapply(available_biomarkers, function(g) {
  df_g <- gene_expr_df %>% filter(gene == g)
  plot_group_box_simple(df_g, g)
}) %>% bind_rows()

group_test_results$p.adj <- p.adjust(group_test_results$p.value, method = "BH")

write.csv(
  group_test_results,
  file = file.path(outdir, "WT_control_vs_ISA_wilcox_biomarkers_sex_stratified.csv"),
  row.names = FALSE
)

# ==============
# CORRELATIONS
# ==============

calc_cor_both <- function(sample_selector_name, selector_vector) {
  samples <- rownames(sample_info)[selector_vector]
  folh1_rows <- get_symbol_rows("FOLH1")
  folh1_expr <- as.numeric(colMeans(vst_mat[folh1_rows, samples, drop = FALSE], na.rm = TRUE))
  
  res <- lapply(setdiff(available_biomarkers, "FOLH1"), function(g) {
    rows <- get_symbol_rows(g)
    x <- as.numeric(colMeans(vst_mat[rows, samples, drop = FALSE], na.rm = TRUE))
    ok <- is.finite(x) & is.finite(folh1_expr)
    
    if (sum(ok) < 3) {
      return(data.frame(
        gene = g,
        subset_name = sample_selector_name,
        pearson = NA_real_,
        ppearson = NA_real_,
        spearman = NA_real_,
        pspearman = NA_real_,
        n = sum(ok),
        stringsAsFactors = FALSE
      ))
    }
    
    ct_p <- tryCatch(cor.test(x[ok], folh1_expr[ok], method = "pearson"), error = function(e) NULL)
    ct_s <- tryCatch(cor.test(x[ok], folh1_expr[ok], method = "spearman"), error = function(e) NULL)
    
    data.frame(
      gene = g,
      subset_name = sample_selector_name,
      pearson = if (!is.null(ct_p)) unname(ct_p$estimate) else NA_real_,
      ppearson = if (!is.null(ct_p)) ct_p$p.value else NA_real_,
      spearman = if (!is.null(ct_s)) unname(ct_s$estimate) else NA_real_,
      pspearman = if (!is.null(ct_s)) ct_s$p.value else NA_real_,
      n = sum(ok),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  res$padj_pearson <- p.adjust(res$ppearson, method = "BH")
  res$padj_spearman <- p.adjust(res$pspearman, method = "BH")
  res %>% arrange(desc(pearson))
}

cor_control <- calc_cor_both("control", sample_info$treatment == "control")
cor_ISA <- calc_cor_both("ISA", sample_info$treatment == "ISA")
cor_male_control <- calc_cor_both("male_control", sample_info$sex == "male" & sample_info$treatment == "control")
cor_male_ISA <- calc_cor_both("male_ISA", sample_info$sex == "male" & sample_info$treatment == "ISA")
cor_female_control <- calc_cor_both("female_control", sample_info$sex == "female" & sample_info$treatment == "control")
cor_female_ISA <- calc_cor_both("female_ISA", sample_info$sex == "female" & sample_info$treatment == "ISA")

cor_all <- bind_rows(
  cor_control,
  cor_ISA,
  cor_male_control,
  cor_male_ISA,
  cor_female_control,
  cor_female_ISA
)

write.csv(
  cor_all,
  file = file.path(outdir, "FOLH1_correlations_ISA_sex_stratified.csv"),
  row.names = FALSE
)

plot_correlation_bar_pearson <- function(cor_df, subset_name, file_stub) {
  plot_df <- cor_df %>%
    filter(!is.na(pearson)) %>%
    mutate(
      direction = ifelse(pearson >= 0, "Positive", "Negative"),
      gene = factor(gene, levels = rev(gene))
    )
  
  if (nrow(plot_df) == 0) return(NULL)
  
  p <- ggplot(plot_df, aes(x = pearson, y = gene, fill = direction)) +
    geom_col(width = 0.8, colour = "white") +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.7) +
    scale_fill_manual(values = c("Positive" = "#E64B35FF", "Negative" = "#4DBBD5FF")) +
    labs(
      title = paste0("Pearson correlation of FOLH1 with known biomarkers (", subset_name, ")"),
      x = "Pearson r with FOLH1",
      y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    file.path(corr_bar_dir, paste0("FOLH1_biomarker_correlation_barplot_Pearson_", file_stub, ".png")),
    p,
    width = 9,
    height = 7,
    dpi = 300
  )
}

plot_correlation_bar_pearson(cor_control, "control", "control")
plot_correlation_bar_pearson(cor_ISA, "ISA", "ISA")
plot_correlation_bar_pearson(cor_male_control, "male_control", "male_control")
plot_correlation_bar_pearson(cor_male_ISA, "male_ISA", "male_ISA")
plot_correlation_bar_pearson(cor_female_control, "female_control", "female_control")
plot_correlation_bar_pearson(cor_female_ISA, "female_ISA", "female_ISA")

# =============
# SCATTERPLOTS
# =============

scatter_fun_group <- function(g, subset_name, selector_vector, save_dir) {
  safe_dir(save_dir)
  samples <- rownames(sample_info)[selector_vector]
  
  folh1_rows <- get_symbol_rows("FOLH1")
  folh1_expr <- as.numeric(colMeans(vst_mat[folh1_rows, samples, drop = FALSE], na.rm = TRUE))
  
  rows <- get_symbol_rows(g)
  biomarker_expr <- as.numeric(colMeans(vst_mat[rows, samples, drop = FALSE], na.rm = TRUE))
  
  df <- data.frame(
    sample = samples,
    sample_short = sample_info[samples, "sample_short"],
    FOLH1 = folh1_expr,
    biomarker = biomarker_expr,
    subset_name = subset_name,
    stringsAsFactors = FALSE
  )
  
  ok <- is.finite(df$FOLH1) & is.finite(df$biomarker)
  df <- df[ok, , drop = FALSE]
  if (nrow(df) < 3) return(NULL)
  
  ct <- suppressWarnings(cor.test(df$FOLH1, df$biomarker, method = "spearman", exact = FALSE))
  lab <- paste0("Spearman rho = ", round(unname(ct$estimate), 3), ", p = ", signif(ct$p.value, 3))
  
  point_col <- dplyr::case_when(
    subset_name == "control" ~ "#F8766D",
    subset_name == "ISA" ~ "#7CAE00",
    subset_name == "male_control" ~ "#E64B35FF",
    subset_name == "male_ISA" ~ "#4DBBD5FF",
    subset_name == "female_control" ~ "#00A087FF",
    subset_name == "female_ISA" ~ "#3C5488FF",
    TRUE ~ "#333333"
  )
  
  p <- ggplot(df, aes(x = FOLH1, y = biomarker)) +
    geom_point(size = 2.8, alpha = 0.9, color = point_col) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.85) +
    labs(
      title = paste0("FOLH1 vs ", g, " (", subset_name, ")"),
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
    file.path(save_dir, paste0("FOLH1_vs_", g, "_scatter_", subset_name, ".png")),
    p,
    width = 6.2,
    height = 5.2,
    dpi = 300
  )
  
  invisible(df)
}

scatter_subsets <- list(
  control = sample_info$treatment == "control",
  ISA = sample_info$treatment == "ISA",
  male_control = sample_info$sex == "male" & sample_info$treatment == "control",
  male_ISA = sample_info$sex == "male" & sample_info$treatment == "ISA",
  female_control = sample_info$sex == "female" & sample_info$treatment == "control",
  female_ISA = sample_info$sex == "female" & sample_info$treatment == "ISA"
)

for (nm in names(scatter_subsets)) {
  save_dir <- file.path(scatter_root_dir, nm)
  for (g in available_biomarkers[available_biomarkers != "FOLH1"]) {
    scatter_fun_group(g, nm, scatter_subsets[[nm]], save_dir)
  }
}

# ===================
# BIOMARKER HEATMAP
# ===================

heat_genes <- available_biomarkers
heat_rows <- unique(unlist(lapply(heat_genes, get_symbol_rows)))
heat_mat <- vst_mat[heat_rows, , drop = FALSE]
heat_symbols <- row_map$SYMBOL[match(rownames(heat_mat), row_map$ENSEMBL)]
heat_symbols <- ifelse(is.na(heat_symbols) | heat_symbols == "", rownames(heat_mat), heat_symbols)
rownames(heat_mat) <- make.unique(as.character(heat_symbols))

heat_anno_col <- data.frame(
  sex = sample_info$sex,
  treatment = sample_info$treatment
)
rownames(heat_anno_col) <- rownames(sample_info)

safe_pheatmap_file(
  heat_mat,
  scale = "row",
  annotation_col = heat_anno_col,
  labels_col = sample_info[colnames(heat_mat), "sample_short"],
  show_colnames = TRUE,
  main = "Biomarker panel (ISA sex-stratified)",
  filename = file.path(outdir, "biomarker_heatmap_ISA_sex_stratified.png"),
  width = 10,
  height = 8
)

# =====================
# MSIGDB SAFE HELPERS
# =====================

get_msig_term2gene <- function(collection_name) {
  fetch_msig <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    cn <- colnames(df)
    
    gene_col <- dplyr::case_when(
      "entrez_gene" %in% cn ~ "entrez_gene",
      "ncbi_gene" %in% cn ~ "ncbi_gene",
      "gene_symbol" %in% cn ~ "gene_symbol",
      "gene" %in% cn ~ "gene",
      TRUE ~ NA_character_
    )
    
    term_col <- dplyr::case_when(
      "gs_name" %in% cn ~ "gs_name",
      "term" %in% cn ~ "term",
      TRUE ~ NA_character_
    )
    
    if (is.na(gene_col) || is.na(term_col)) return(NULL)
    
    out <- df[, c(term_col, gene_col), drop = FALSE]
    colnames(out) <- c("gs_name", "gene")
    out <- out[complete.cases(out), , drop = FALSE]
    out <- unique(out)
    out
  }
  
  candidates <- list(
    tryCatch(msigdbr(db_species = "MM", species = "Mus musculus", collection = collection_name), error = function(e) NULL),
    tryCatch(msigdbr(species = "Mus musculus", collection = collection_name), error = function(e) NULL),
    tryCatch(msigdbr(species = "mouse", collection = collection_name), error = function(e) NULL)
  )
  
  for (obj in candidates) {
    got <- fetch_msig(obj)
    if (!is.null(got) && nrow(got) > 0) return(got)
  }
  
  NULL
}

hallmark_term2gene <- get_msig_term2gene("MH")
immune_term2gene <- get_msig_term2gene("M7")

if (is.null(hallmark_term2gene) || nrow(hallmark_term2gene) == 0) {
  message("Hallmark mouse term2gene is empty; Hallmark analyses will be skipped.")
}

if (is.null(immune_term2gene) || nrow(immune_term2gene) == 0) {
  message("ImmuneSigDB mouse term2gene is empty; ImmuneSigDB analyses will be skipped.")
}

# ==============================
# HELPERS FOR PAIRWISE ANALYSIS
# ==============================

make_deg_heatmap <- function(gene_table, mat, sample_info_sub, filename, title_text, scale_rows = TRUE) {
  genes_use <- intersect(gene_table$ENSEMBL, rownames(mat))
  if (length(genes_use) < 2) return(NULL)
  
  hm_mat <- mat[genes_use, , drop = FALSE]
  gene_symbols <- gene_table$SYMBOL[match(rownames(hm_mat), gene_table$ENSEMBL)]
  gene_symbols <- ifelse(is.na(gene_symbols) | gene_symbols == "", rownames(hm_mat), gene_symbols)
  rownames(hm_mat) <- make.unique(as.character(gene_symbols))
  
  ann_col <- data.frame(
    group = sample_info_sub[colnames(hm_mat), "group"],
    row.names = colnames(hm_mat)
  )
  
  safe_pheatmap_file(
    hm_mat,
    scale = if (scale_rows) "row" else "none",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = ann_col,
    labels_col = sample_info_sub[colnames(hm_mat), "sample_short"],
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = ifelse(nrow(hm_mat) > 80, 3.5, 6),
    fontsize_col = 8,
    main = title_text,
    filename = filename,
    width = 10,
    height = 12
  )
}

run_deg_and_pathway_analysis <- function(
    count_mat,
    vst_mat,
    sample_info,
    row_map,
    sample_selector,
    group_values,
    comparison_name,
    numerator_label,
    denominator_label,
    deg_root_dir,
    enrich_root_dir
) {
  sample_info_sub <- sample_info[sample_selector, , drop = FALSE]
  count_sub <- count_mat[, rownames(sample_info_sub), drop = FALSE]
  vst_sub <- vst_mat[, rownames(sample_info_sub), drop = FALSE]
  
  if (length(group_values) != nrow(sample_info_sub)) {
    stop(paste0("Length mismatch in group_values for comparison: ", comparison_name))
  }
  
  sample_info_sub$group <- factor(as.character(group_values), levels = c(denominator_label, numerator_label))
  
  if (any(is.na(sample_info_sub$group))) {
    stop(paste0("NA values created in group factor for comparison: ", comparison_name))
  }
  
  if (length(unique(sample_info_sub$group)) < 2) {
    stop(paste0("Comparison ", comparison_name, " does not contain two valid groups."))
  }
  
  rownames(sample_info_sub) <- sample_info_sub$sample
  
  dds_sub <- DESeqDataSetFromMatrix(
    countData = round(count_sub),
    colData = sample_info_sub,
    design = ~ group
  )
  
  dds_sub <- dds_sub[rowSums(counts(dds_sub)) > 1, ]
  dds_sub$group <- droplevels(dds_sub$group)
  dds_sub <- DESeq(dds_sub)
  
  deg_dir <- file.path(deg_root_dir, comparison_name)
  deg_table_dir <- file.path(deg_dir, "01_DEG_tables")
  deg_heatmap_dir <- file.path(deg_dir, "02_DEG_heatmaps")
  
  enrich_dir <- file.path(enrich_root_dir, comparison_name)
  go_bp_dir <- file.path(enrich_dir, "01_GO_BP")
  hallmark_dir <- file.path(enrich_dir, "02_Hallmark")
  gsea_dir <- file.path(enrich_dir, "03_GSEA_ranked")
  immune_dir <- file.path(enrich_dir, "04_ImmuneSigDB")
  
  safe_dir(deg_dir)
  safe_dir(deg_table_dir)
  safe_dir(deg_heatmap_dir)
  safe_dir(enrich_dir)
  safe_dir(go_bp_dir)
  safe_dir(hallmark_dir)
  safe_dir(gsea_dir)
  safe_dir(immune_dir)
  
  res <- results(dds_sub, contrast = c("group", numerator_label, denominator_label), alpha = 0.05)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(row_map %>% select(ENSEMBL, SYMBOL, ENTREZID, GENENAME), by = "ENSEMBL") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  res_df$regulated <- case_when(
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0.58 ~ paste0("Up_in_", numerator_label),
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.58 ~ paste0("Down_in_", numerator_label),
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
    Comparison = comparison_name,
    Total_DE_Genes = nrow(deg_sig),
    Upregulated = sum(deg_sig$log2FoldChange > 0.58, na.rm = TRUE),
    Downregulated = sum(deg_sig$log2FoldChange < -0.58, na.rm = TRUE),
    padj_Cutoff = 0.05,
    log2FC_Cutoff = 0.58,
    stringsAsFactors = FALSE
  )
  
  write.csv(res_df, file.path(deg_table_dir, paste0("DESeq2_all_results_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(deg_sig, file.path(deg_table_dir, paste0("DESeq2_significant_DEGs_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(deg_top50, file.path(deg_table_dir, paste0("DESeq2_top50_ranked_genes_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(deg_top100, file.path(deg_table_dir, paste0("DESeq2_top100_ranked_genes_", comparison_name, ".csv")), row.names = FALSE)
  write.csv(deg_summary, file.path(deg_table_dir, paste0("DEG_summary_", comparison_name, ".csv")), row.names = FALSE)
  
  make_deg_heatmap(
    gene_table = deg_sig,
    mat = vst_sub,
    sample_info_sub = sample_info_sub,
    filename = file.path(deg_heatmap_dir, paste0("Heatmap_all_significant_DEGs_", comparison_name, ".png")),
    title_text = paste0("Heatmap of all significant DEGs (", comparison_name, ")")
  )
  
  make_deg_heatmap(
    gene_table = deg_top50,
    mat = vst_sub,
    sample_info_sub = sample_info_sub,
    filename = file.path(deg_heatmap_dir, paste0("Heatmap_top50_DEGs_", comparison_name, ".png")),
    title_text = paste0("Heatmap of top 50 DEGs (", comparison_name, ")")
  )
  
  make_deg_heatmap(
    gene_table = deg_top100,
    mat = vst_sub,
    sample_info_sub = sample_info_sub,
    filename = file.path(deg_heatmap_dir, paste0("Heatmap_top100_DEGs_", comparison_name, ".png")),
    title_text = paste0("Heatmap of top 100 DEGs (", comparison_name, ")")
  )
  
  res_ranked <- res_df %>%
    filter(!is.na(log2FoldChange)) %>%
    mutate(rank_metric = log2FoldChange)
  
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
  
  ego_up <- if (length(sig_entrez_up) > 0) {
    tryCatch(
      enrichGO(
        gene = sig_entrez_up,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      ),
      error = function(e) NULL
    )
  } else NULL
  
  ego_down <- if (length(sig_entrez_down) > 0) {
    tryCatch(
      enrichGO(
        gene = sig_entrez_down,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      ),
      error = function(e) NULL
    )
  } else NULL
  
  if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
    write.csv(as.data.frame(ego_up), file.path(go_bp_dir, paste0("GO_BP_ORA_up_in_", numerator_label, "_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(go_bp_dir, paste0("GO_BP_ORA_up_in_", numerator_label, "_dotplot_", comparison_name, ".png")),
      print(dotplot(ego_up, showCategory = 20) + ggtitle(paste0("GO BP enriched in upregulated genes (", comparison_name, ")")))
    )
  }
  
  if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
    write.csv(as.data.frame(ego_down), file.path(go_bp_dir, paste0("GO_BP_ORA_down_in_", numerator_label, "_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(go_bp_dir, paste0("GO_BP_ORA_down_in_", numerator_label, "_dotplot_", comparison_name, ".png")),
      print(dotplot(ego_down, showCategory = 20) + ggtitle(paste0("GO BP enriched in downregulated genes (", comparison_name, ")")))
    )
  }
  
  hallmark_up <- if (!is.null(hallmark_term2gene) && nrow(hallmark_term2gene) > 0 && length(sig_entrez_up) > 0) {
    tryCatch(
      enricher(
        gene = sig_entrez_up,
        TERM2GENE = hallmark_term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      ),
      error = function(e) NULL
    )
  } else NULL
  
  hallmark_down <- if (!is.null(hallmark_term2gene) && nrow(hallmark_term2gene) > 0 && length(sig_entrez_down) > 0) {
    tryCatch(
      enricher(
        gene = sig_entrez_down,
        TERM2GENE = hallmark_term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      ),
      error = function(e) NULL
    )
  } else NULL
  
  gsea_hallmark <- if (!is.null(hallmark_term2gene) && nrow(hallmark_term2gene) > 0 && length(gene_list) > 0) {
    tryCatch(
      GSEA(
        geneList = gene_list,
        TERM2GENE = hallmark_term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
  } else NULL
  
  if (!is.null(hallmark_up) && nrow(as.data.frame(hallmark_up)) > 0) {
    write.csv(as.data.frame(hallmark_up), file.path(hallmark_dir, paste0("Hallmark_ORA_up_in_", numerator_label, "_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(hallmark_dir, paste0("Hallmark_ORA_up_in_", numerator_label, "_dotplot_", comparison_name, ".png")),
      print(dotplot(hallmark_up, showCategory = 20) + ggtitle(paste0("Hallmark enriched in upregulated genes (", comparison_name, ")")))
    )
  }
  
  if (!is.null(hallmark_down) && nrow(as.data.frame(hallmark_down)) > 0) {
    write.csv(as.data.frame(hallmark_down), file.path(hallmark_dir, paste0("Hallmark_ORA_down_in_", numerator_label, "_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(hallmark_dir, paste0("Hallmark_ORA_down_in_", numerator_label, "_dotplot_", comparison_name, ".png")),
      print(dotplot(hallmark_down, showCategory = 20) + ggtitle(paste0("Hallmark enriched in downregulated genes (", comparison_name, ")")))
    )
  }
  
  if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
    write.csv(as.data.frame(gsea_hallmark), file.path(gsea_dir, paste0("GSEA_Hallmark_ranked_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(gsea_dir, paste0("GSEA_Hallmark_dotplot_", comparison_name, ".png")),
      print(dotplot(gsea_hallmark, showCategory = 20, split = ".sign") +
              facet_grid(. ~ .sign) +
              ggtitle(paste0("Hallmark GSEA ranked analysis (", comparison_name, ")")))
    )
  }
  
  immune_up <- if (!is.null(immune_term2gene) && nrow(immune_term2gene) > 0 && length(sig_entrez_up) > 0) {
    tryCatch(
      enricher(
        gene = sig_entrez_up,
        TERM2GENE = immune_term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      ),
      error = function(e) NULL
    )
  } else NULL
  
  immune_down <- if (!is.null(immune_term2gene) && nrow(immune_term2gene) > 0 && length(sig_entrez_down) > 0) {
    tryCatch(
      enricher(
        gene = sig_entrez_down,
        TERM2GENE = immune_term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      ),
      error = function(e) NULL
    )
  } else NULL
  
  gsea_immune <- if (!is.null(immune_term2gene) && nrow(immune_term2gene) > 0 && length(gene_list) > 0) {
    tryCatch(
      GSEA(
        geneList = gene_list,
        TERM2GENE = immune_term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
  } else NULL
  
  if (!is.null(immune_up) && nrow(as.data.frame(immune_up)) > 0) {
    write.csv(as.data.frame(immune_up), file.path(immune_dir, paste0("ImmuneSigDB_ORA_up_in_", numerator_label, "_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(immune_dir, paste0("ImmuneSigDB_ORA_up_in_", numerator_label, "_dotplot_", comparison_name, ".png")),
      print(dotplot(immune_up, showCategory = 20) + ggtitle(paste0("ImmuneSigDB enriched in upregulated genes (", comparison_name, ")")))
    )
  }
  
  if (!is.null(immune_down) && nrow(as.data.frame(immune_down)) > 0) {
    write.csv(as.data.frame(immune_down), file.path(immune_dir, paste0("ImmuneSigDB_ORA_down_in_", numerator_label, "_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(immune_dir, paste0("ImmuneSigDB_ORA_down_in_", numerator_label, "_dotplot_", comparison_name, ".png")),
      print(dotplot(immune_down, showCategory = 20) + ggtitle(paste0("ImmuneSigDB enriched in downregulated genes (", comparison_name, ")")))
    )
  }
  
  if (!is.null(gsea_immune) && nrow(as.data.frame(gsea_immune)) > 0) {
    write.csv(as.data.frame(gsea_immune), file.path(gsea_dir, paste0("GSEA_ImmuneSigDB_ranked_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(gsea_dir, paste0("GSEA_ImmuneSigDB_dotplot_", comparison_name, ".png")),
      print(dotplot(gsea_immune, showCategory = 20, split = ".sign") +
              facet_grid(. ~ .sign) +
              ggtitle(paste0("ImmuneSigDB GSEA ranked analysis (", comparison_name, ")")))
    )
  }
  
  gsea_go_bp <- if (length(gene_list) > 0) {
    tryCatch(
      gseGO(
        geneList = gene_list,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
  } else NULL
  
  if (!is.null(gsea_go_bp) && nrow(as.data.frame(gsea_go_bp)) > 0) {
    write.csv(as.data.frame(gsea_go_bp), file.path(gsea_dir, paste0("GSEA_GO_BP_ranked_", comparison_name, ".csv")), row.names = FALSE)
    safe_png_plot(
      file.path(gsea_dir, paste0("GSEA_GO_BP_dotplot_", comparison_name, ".png")),
      print(dotplot(gsea_go_bp, showCategory = 20, split = ".sign") +
              facet_grid(. ~ .sign) +
              ggtitle(paste0("GO BP GSEA ranked analysis (", comparison_name, ")")))
    )
  }
  
  deg_summary_tbl <- tableGrob(
    deg_summary,
    rows = NULL,
    theme = ttheme_minimal(base_size = 11)
  )
  
  safe_png_plot(
    file.path(deg_dir, paste0("DEG_summary_table_", comparison_name, ".png")),
    {
      grid.newpage()
      grid.draw(deg_summary_tbl)
    },
    width = 2600,
    height = 700,
    res = 300
  )
  
  list(
    comparison_name = comparison_name,
    sample_info = sample_info_sub,
    results = res_df,
    deg_sig = deg_sig,
    deg_top50 = deg_top50,
    deg_top100 = deg_top100,
    deg_summary = deg_summary,
    ego_up = ego_up,
    ego_down = ego_down,
    hallmark_up = hallmark_up,
    hallmark_down = hallmark_down,
    gsea_hallmark = gsea_hallmark,
    immune_up = immune_up,
    immune_down = immune_down,
    gsea_immune = gsea_immune,
    gsea_go_bp = gsea_go_bp
  )
}

# ====================
# PAIRWISE COMPARISONS
# ====================

comparison_results <- list(
  control_vs_ISA = run_deg_and_pathway_analysis(
    count_mat = count_mat,
    vst_mat = vst_mat,
    sample_info = sample_info,
    row_map = row_map,
    sample_selector = rep(TRUE, nrow(sample_info)),
    group_values = sample_info$treatment,
    comparison_name = "control_vs_ISA",
    numerator_label = "ISA",
    denominator_label = "control",
    deg_root_dir = deg_root_dir,
    enrich_root_dir = enrich_root_dir
  ),
  
  male_control_vs_male_ISA = run_deg_and_pathway_analysis(
    count_mat = count_mat,
    vst_mat = vst_mat,
    sample_info = sample_info,
    row_map = row_map,
    sample_selector = sample_info$sex == "male",
    group_values = sample_info$treatment[sample_info$sex == "male"],
    comparison_name = "male_control_vs_male_ISA",
    numerator_label = "ISA",
    denominator_label = "control",
    deg_root_dir = deg_root_dir,
    enrich_root_dir = enrich_root_dir
  ),
  
  female_control_vs_female_ISA = run_deg_and_pathway_analysis(
    count_mat = count_mat,
    vst_mat = vst_mat,
    sample_info = sample_info,
    row_map = row_map,
    sample_selector = sample_info$sex == "female",
    group_values = sample_info$treatment[sample_info$sex == "female"],
    comparison_name = "female_control_vs_female_ISA",
    numerator_label = "ISA",
    denominator_label = "control",
    deg_root_dir = deg_root_dir,
    enrich_root_dir = enrich_root_dir
  ),
  
  male_ISA_vs_female_ISA = run_deg_and_pathway_analysis(
    count_mat = count_mat,
    vst_mat = vst_mat,
    sample_info = sample_info,
    row_map = row_map,
    sample_selector = sample_info$treatment == "ISA",
    group_values = sample_info$sex[sample_info$treatment == "ISA"],
    comparison_name = "male_ISA_vs_female_ISA",
    numerator_label = "male",
    denominator_label = "female",
    deg_root_dir = deg_root_dir,
    enrich_root_dir = enrich_root_dir
  ),
  
  male_control_vs_female_control = run_deg_and_pathway_analysis(
    count_mat = count_mat,
    vst_mat = vst_mat,
    sample_info = sample_info,
    row_map = row_map,
    sample_selector = sample_info$treatment == "control",
    group_values = sample_info$sex[sample_info$treatment == "control"],
    comparison_name = "male_control_vs_female_control",
    numerator_label = "male",
    denominator_label = "female",
    deg_root_dir = deg_root_dir,
    enrich_root_dir = enrich_root_dir
  )
)

# ===============
# EXCEL WORKBOOK
# ===============

wb <- createWorkbook()

addWorksheet(wb, "Sample_info")
writeData(wb, "Sample_info", sample_info)

addWorksheet(wb, "PCA")
writeData(wb, "PCA", pca_df)

addWorksheet(wb, "Biomarker_presence")
writeData(wb, "Biomarker_presence", biomarker_presence)

addWorksheet(wb, "Expression_long")
writeData(wb, "Expression_long", gene_expr_df)

addWorksheet(wb, "Wilcox_results")
writeData(wb, "Wilcox_results", group_test_results)

addWorksheet(wb, "Correlation_all")
writeData(wb, "Correlation_all", cor_all)

safe_sheet_name <- function(x, suffix = "") {
  nm <- paste0(x, suffix)
  nm <- gsub("[\\\\/:*?\\[\\]]", "_", nm)
  substr(nm, 1, 31)
}

for (nm in names(comparison_results)) {
  comp <- comparison_results[[nm]]
  
  sh1 <- safe_sheet_name(nm, "_summary")
  addWorksheet(wb, sh1)
  writeData(wb, sh1, comp$deg_summary)
  
  sh2 <- safe_sheet_name(nm, "_all")
  addWorksheet(wb, sh2)
  writeData(wb, sh2, comp$results)
  
  sh3 <- safe_sheet_name(nm, "_sig")
  addWorksheet(wb, sh3)
  writeData(wb, sh3, comp$deg_sig)
  
  sh4 <- safe_sheet_name(nm, "_top50")
  addWorksheet(wb, sh4)
  writeData(wb, sh4, comp$deg_top50)
  
  sh5 <- safe_sheet_name(nm, "_top100")
  addWorksheet(wb, sh5)
  writeData(wb, sh5, comp$deg_top100)
  
  if (!is.null(comp$ego_up) && nrow(as.data.frame(comp$ego_up)) > 0) {
    sh <- safe_sheet_name(nm, "_GO_up")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$ego_up))
  }
  if (!is.null(comp$ego_down) && nrow(as.data.frame(comp$ego_down)) > 0) {
    sh <- safe_sheet_name(nm, "_GO_down")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$ego_down))
  }
  if (!is.null(comp$hallmark_up) && nrow(as.data.frame(comp$hallmark_up)) > 0) {
    sh <- safe_sheet_name(nm, "_Hallmark_up")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$hallmark_up))
  }
  if (!is.null(comp$hallmark_down) && nrow(as.data.frame(comp$hallmark_down)) > 0) {
    sh <- safe_sheet_name(nm, "_Hallmark_down")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$hallmark_down))
  }
  if (!is.null(comp$gsea_hallmark) && nrow(as.data.frame(comp$gsea_hallmark)) > 0) {
    sh <- safe_sheet_name(nm, "_GSEA_Hallmark")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$gsea_hallmark))
  }
  if (!is.null(comp$immune_up) && nrow(as.data.frame(comp$immune_up)) > 0) {
    sh <- safe_sheet_name(nm, "_Immune_up")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$immune_up))
  }
  if (!is.null(comp$immune_down) && nrow(as.data.frame(comp$immune_down)) > 0) {
    sh <- safe_sheet_name(nm, "_Immune_down")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$immune_down))
  }
  if (!is.null(comp$gsea_immune) && nrow(as.data.frame(comp$gsea_immune)) > 0) {
    sh <- safe_sheet_name(nm, "_GSEA_Immune")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$gsea_immune))
  }
  if (!is.null(comp$gsea_go_bp) && nrow(as.data.frame(comp$gsea_go_bp)) > 0) {
    sh <- safe_sheet_name(nm, "_GSEA_GO_BP")
    addWorksheet(wb, sh)
    writeData(wb, sh, as.data.frame(comp$gsea_go_bp))
  }
}

saveWorkbook(
  wb,
  file = file.path(outdir, "Overall_mouse_ISA_biomarker_and_DEG_results_sex_stratified.xlsx"),
  overwrite = TRUE
)

cat("ISA sex-stratified pipeline complete.\n")