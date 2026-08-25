# ====================================================================
# Novogene 2026 Mice Analysis
# Original Tumor (Male vs Female)
# Clustering, DEG, Pathway enrichment, Expression of biomarkers, mMCP
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
  library(readxl)
  library(purrr)
})

# =====================
# DIRECTORIES / LABELS
# =====================

counts_dir        <- "../counts"
sample_info_file  <- file.path(counts_dir, "Sample information 20260516.xlsx")

outdir            <- "Male_v_Female_in_KO_WT"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

qc_dir            <- file.path(outdir, "01_QC");                            dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

deg_dir           <- file.path(outdir, "02_DEG_analysis");                  dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)
deg_WT_dir        <- file.path(deg_dir, "WT");                              dir.create(deg_WT_dir, showWarnings = FALSE, recursive = TRUE)
deg_KO_dir        <- file.path(deg_dir, "KO");                              dir.create(deg_KO_dir, showWarnings = FALSE, recursive = TRUE)

deg_table_WT_dir  <- file.path(deg_WT_dir, "01_DEG_tables");                dir.create(deg_table_WT_dir, showWarnings = FALSE, recursive = TRUE)
deg_heatmap_WT_dir<- file.path(deg_WT_dir, "02_DEG_heatmaps");              dir.create(deg_heatmap_WT_dir, showWarnings = FALSE, recursive = TRUE)

deg_table_KO_dir  <- file.path(deg_KO_dir, "01_DEG_tables");                dir.create(deg_table_KO_dir, showWarnings = FALSE, recursive = TRUE)
deg_heatmap_KO_dir<- file.path(deg_KO_dir, "02_DEG_heatmaps");              dir.create(deg_heatmap_KO_dir, showWarnings = FALSE, recursive = TRUE)

expr_root_dir     <- file.path(outdir, "03_Expression_boxplots");           dir.create(expr_root_dir, showWarnings = FALSE, recursive = TRUE)
expr_WT_dir       <- file.path(expr_root_dir, "WT");                        dir.create(expr_WT_dir, showWarnings = FALSE, recursive = TRUE)
expr_KO_dir       <- file.path(expr_root_dir, "KO");                        dir.create(expr_KO_dir, showWarnings = FALSE, recursive = TRUE)

scatter_root_dir  <- file.path(outdir, "04_Correlation_scatterplots");      dir.create(scatter_root_dir, showWarnings = FALSE, recursive = TRUE)
corr_bar_dir      <- file.path(outdir, "05_Correlation_barplots");          dir.create(corr_bar_dir, showWarnings = FALSE, recursive = TRUE)
pearson_bar_dir   <- file.path(outdir, "05b_Pearson_correlation_barplots"); dir.create(pearson_bar_dir, showWarnings = FALSE, recursive = TRUE)

enrich_dir        <- file.path(outdir, "06_Pathway_analysis");              dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)
go_bp_dir         <- file.path(enrich_dir, "01_GO_BP_ORA");                 dir.create(go_bp_dir, showWarnings = FALSE, recursive = TRUE)
hallmark_gsea_dir <- file.path(enrich_dir, "02_Hallmark_GSEA");             dir.create(hallmark_gsea_dir, showWarnings = FALSE, recursive = TRUE)
immune_gsea_dir   <- file.path(enrich_dir, "03_ImmuneSigDB_GSEA");          dir.create(immune_gsea_dir, showWarnings = FALSE, recursive = TRUE)

mmcp_root_dir     <- file.path(outdir, "07_mMCP_counter");                  dir.create(mmcp_root_dir, showWarnings = FALSE, recursive = TRUE)
mmcp_input_dir    <- file.path(mmcp_root_dir, "01_input")
mmcp_scores_dir   <- file.path(mmcp_root_dir, "02_scores")
mmcp_group_dir    <- file.path(mmcp_root_dir, "03_group_comparisons")
mmcp_plot_dir     <- file.path(mmcp_root_dir, "04_plots")
dir.create(mmcp_input_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(mmcp_scores_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(mmcp_group_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(mmcp_plot_dir,   showWarnings = FALSE, recursive = TRUE)

tumor_label       <- "control samples"
comparison_label  <- "Male vs Female within genotype"

# =====================
# THEMES / COLORS
# =====================

theme_nature <- function(base_size = 14, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title      = element_text(face = "bold", colour = "black"),
      axis.text       = element_text(colour = "black"),
      axis.line       = element_line(linewidth = 0.8, colour = "black"),
      axis.ticks      = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks.length = unit(0.18, "cm"),
      legend.title    = element_blank(),
      legend.background = element_blank(),
      legend.key      = element_blank(),
      panel.border    = element_blank(),
      panel.grid      = element_blank(),
      strip.background= element_rect(fill = "grey92", colour = "black", linewidth = 0.8),
      strip.text      = element_text(face = "bold", colour = "black")
    )
}
theme_set(theme_nature())

group_cols_box   <- c("WT" = "#F8766D", "KO" = "#00BFC4")
group_cols_qc    <- c("WT" = "#B22222", "KO" = "#0B3C8C")
sex_cols         <- c("Male" = "#1f78b4", "Female" = "#e31a1c")
corr_bar_cols    <- c("Positive" = "#E64B35FF", "Negative" = "#4DBBD5FF")
pearson_bar_cols <- c("Positive" = "#F18D9E", "Negative" = "#7FB3D5")

theme_enrich_medium <- function() {
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 7, colour = "black"),
    axis.text.y = element_text(size = 6, colour = "black"),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7),
    strip.text   = element_text(size = 7, face = "bold")
  )
}

theme_enrich_small <- function() {
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 7, colour = "black"),
    axis.text.y = element_text(size = 4, colour = "black"),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7),
    strip.text   = element_text(size = 7, face = "bold")
  )
}

short_sample_name <- function(x) {
  paste0(substr(x, 1, 7), "_", substr(x, nchar(x) - 1, nchar(x)))
}

# ====================
# 1) SAMPLES TO KEEP
# ====================

target_samples <- c(
  "WT2921C_OE", "WT2910C_OE", "WT2898C_OE", "WT2965C_OE",
  "KO2912C_OE", "KO2918C_OE", "KO2960C_OE", "KO2913C_OE",
  "C2890WT", "C2917WT", "C2893WT",
  "WT2868C", "WT2881C", "WT2882C",
  "KO2868C", "KO2881C", "KO2882C",
  "T1728C", "T1729C", "T1740C", "T1742C"
)

# ======================================
# 2) READ SAMPLE INFO (HEADER ON ROW 2)
# ======================================

sample_sheet <- read_excel(sample_info_file, sheet = 1, skip = 1)
colnames(sample_sheet) <- trimws(colnames(sample_sheet))
sample_sheet <- sample_sheet %>% mutate(across(everything(), as.character))

code_col    <- grep("^code$|code", colnames(sample_sheet), ignore.case = TRUE, value = TRUE)[1]
gender_col  <- grep("gender", colnames(sample_sheet), ignore.case = TRUE, value = TRUE)[1]
samplenr_col<- grep("sample nr|sample", colnames(sample_sheet), ignore.case = TRUE, value = TRUE)[1]

if (is.na(code_col))   stop("Could not find 'Code' column in sample info.")
if (is.na(gender_col)) stop("Could not find 'Gender' column in sample info.")
if (is.na(samplenr_col)) warning("Could not clearly find 'Sample nr' column, using first match.")

sample_meta <- sample_sheet %>%
  transmute(
    Code     = .data[[code_col]],
    Gender   = .data[[gender_col]],
    SampleNr = .data[[samplenr_col]]
  )

sample_meta$Gender <- case_when(
  grepl("^m", sample_meta$Gender, ignore.case = TRUE) ~ "Male",
  grepl("^f", sample_meta$Gender, ignore.case = TRUE) ~ "Female",
  TRUE ~ sample_meta$Gender
)

# ====================================
# 3) READ AND COMBINE ALL COUNT FILES
# ====================================

count_files <- list.files(
  counts_dir,
  pattern = "\\.tabular$",
  full.names = TRUE
)

if (length(count_files) == 0) stop("No .tabular files found in ../counts.")

read_count_file <- function(f) {
  x <- read.delim(
    f,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(x) < 2) stop(paste("File has too few columns:", f))
  
  gene_col <- as.character(x[[1]])
  gene_col <- sub("\\..*$", "", gene_col)
  
  count_cols <- setdiff(colnames(x), colnames(x)[1])
  count_mat  <- x[, count_cols, drop = FALSE]
  
  for (cc in colnames(count_mat)) {
    count_mat[[cc]] <- suppressWarnings(as.numeric(count_mat[[cc]]))
  }
  
  rownames(count_mat) <- gene_col
  count_mat
}

count_list <- lapply(count_files, read_count_file)
names(count_list) <- basename(count_files)

all_genes <- unique(unlist(lapply(count_list, rownames)))

aligned_list <- lapply(count_list, function(mat) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if (length(missing_genes) > 0) {
    na_block <- matrix(
      NA_real_,
      nrow = length(missing_genes),
      ncol = ncol(mat),
      dimnames = list(missing_genes, colnames(mat))
    )
    mat <- rbind(mat, na_block)
  }
  mat[all_genes, , drop = FALSE]
})

combined_counts <- do.call(cbind, aligned_list)

# =========================================================
# 4) FILTER TO TARGET SAMPLES (grep match on column names)
# =========================================================

all_cols    <- colnames(combined_counts)
matched_cols<- character(0)
matched_map <- data.frame(
  sample_id = character(0),
  count_col = character(0),
  stringsAsFactors = FALSE
)

for (sid in target_samples) {
  hits <- grep(sid, all_cols, value = TRUE)
  if (length(hits) == 1) {
    matched_cols <- c(matched_cols, hits)
    matched_map  <- rbind(matched_map,
                          data.frame(sample_id = sid, count_col = hits, stringsAsFactors = FALSE))
  } else if (length(hits) > 1) {
    matched_cols <- c(matched_cols, hits[1])
    matched_map  <- rbind(matched_map,
                          data.frame(sample_id = sid, count_col = hits[1], stringsAsFactors = FALSE))
    warning(paste("Multiple count columns matched", sid, "- using first:", hits[1]))
  } else {
    warning(paste("No count column matched sample ID:", sid))
  }
}

if (length(matched_cols) == 0) stop("None of the target samples were found in the count matrix.")

filtered_counts <- combined_counts[, matched_cols, drop = FALSE]

new_names <- matched_map$sample_id[match(colnames(filtered_counts), matched_map$count_col)]
colnames(filtered_counts) <- new_names
filtered_counts <- filtered_counts[, !duplicated(colnames(filtered_counts)), drop = FALSE]

# =============================================
# 5) BUILD SAMPLE INFO TABLE (genotype, sex)
# =============================================

sample_names <- colnames(filtered_counts)

genotype <- ifelse(
  grepl("KO", sample_names, ignore.case = TRUE), "KO",
  ifelse(grepl("WT", sample_names, ignore.case = TRUE), "WT", "Other")
)

keep_idx <- genotype %in% c("WT", "KO")
filtered_counts <- filtered_counts[, keep_idx, drop = FALSE]
sample_names    <- colnames(filtered_counts)
genotype        <- genotype[keep_idx]
genotype        <- factor(genotype, levels = c("WT", "KO"))

sex_vec <- character(length(sample_names))
names(sex_vec) <- sample_names

for (sid in sample_names) {
  idx_exact <- which(sample_meta$Code == sid)
  if (length(idx_exact) == 1) {
    sex_vec[sid] <- sample_meta$Gender[idx_exact]
    next
  }
  idx_grep <- grep(sid, sample_meta$Code, ignore.case = TRUE)
  if (length(idx_grep) >= 1) {
    sex_vec[sid] <- sample_meta$Gender[idx_grep[1]]
  } else {
    sex_vec[sid] <- NA_character_
    warning(paste("No gender match found for sample ID:", sid))
  }
}

sex_vec <- factor(sex_vec, levels = c("Male", "Female"))

sample_info <- data.frame(
  sample       = sample_names,
  genotype     = genotype,
  sex          = sex_vec,
  sample_type  = tumor_label,
  sample_short = short_sample_name(sample_names),
  stringsAsFactors = FALSE
)
rownames(sample_info) <- sample_names

write.csv(
  sample_info,
  file = file.path(outdir, "sample_info_selected_samples.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = rownames(filtered_counts), filtered_counts, check.names = FALSE),
  file = file.path(outdir, "raw_counts_selected_samples.csv"),
  row.names = FALSE
)

# ===========================================
# 6) DESeq2 OBJECT / NORMALIZATION (global)
# ===========================================

keep_rows       <- rowSums(filtered_counts, na.rm = TRUE) > 1
filtered_counts <- filtered_counts[keep_rows, , drop = FALSE]

dds_full <- DESeqDataSetFromMatrix(
  countData = round(filtered_counts),
  colData   = sample_info,
  design    = ~ genotype + sex
)
dds_full <- dds_full[rowSums(counts(dds_full)) > 1, ]
dds_full <- DESeq(dds_full)

norm_mat <- counts(dds_full, normalized = TRUE)
vst_obj  <- vst(dds_full, blind = TRUE)
vst_mat  <- assay(vst_obj)

# =======================================
# 7) ANNOTATION (ENSEMBL -> SYMBOL etc.)
# =======================================

anno <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys     = rownames(vst_mat),
  keytype  = "ENSEMBL",
  columns  = c("SYMBOL", "ENTREZID", "GENENAME")
)

anno <- anno[!is.na(anno$ENSEMBL), ]
anno <- anno[!duplicated(anno$ENSEMBL), ]

row_map <- data.frame(
  ENSEMBL = rownames(vst_mat),
  stringsAsFactors = FALSE
) %>%
  left_join(anno, by = "ENSEMBL")

row_map$SYMBOL_UPPER <- toupper(row_map$SYMBOL)

# =============================
# 8) PCA QC (genotype + sex)
# =============================

run_pca_plot <- function(mat, sample_info, title_text, file_prefix) {
  mat_use <- mat
  
  pca_obj    <- prcomp(t(mat_use), center = TRUE, scale. = FALSE)
  percent_var<- (pca_obj$sdev^2 / sum(pca_obj$sdev^2)) * 100
  
  pca_df <- data.frame(
    sample       = rownames(pca_obj$x),
    sample_short = short_sample_name(rownames(pca_obj$x)),
    PC1          = pca_obj$x[, 1],
    PC2          = pca_obj$x[, 2],
    genotype     = sample_info[rownames(pca_obj$x), "genotype"],
    sex          = sample_info[rownames(pca_obj$x), "sex"],
    stringsAsFactors = FALSE
  )
  
  pca_df$group <- interaction(pca_df$genotype, pca_df$sex, sep = "_")
  
  group4_cols <- c(
    "WT_Male"   = "red",
    "WT_Female" = "blue",
    "KO_Male"   = "green",
    "KO_Female" = "yellow"
  )
  
  write.csv(
    pca_df,
    file = file.path(qc_dir, paste0(file_prefix, "_coordinates.csv")),
    row.names = FALSE
  )
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    stat_ellipse(
      aes(group = genotype, fill = genotype),
      type  = "norm",
      geom  = "polygon",
      alpha = 0.12,
      colour= NA,
      show.legend = FALSE
    ) +
    geom_point(size = 4, alpha = 0.95) +
    geom_text(
      aes(label = sample_short),
      vjust = -0.8,
      size  = 2.6,
      colour= "black",
      check_overlap = TRUE
    ) +
    scale_color_manual(
      name   = "Genotype × sex",
      values = group4_cols
    ) +
    scale_fill_manual(values = c("WT" = "#377eb8", "KO" = "#e41a1c")) +
    labs(
      title = paste0(title_text, " (", tumor_label, ")"),
      x     = paste0("PC1 (", round(percent_var[1], 1), "%)"),
      y     = paste0("PC2 (", round(percent_var[2], 1), "%)")
    ) +
    theme_nature(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.60, size = 11),
      axis.title = element_text(face = "bold", size = 12),
      axis.text  = element_text(size = 11),
      legend.position = "right",
      plot.margin = unit(c(0.8, 0.8, 0.8, 0.8), "cm")
    )
  
  ggsave(
    file.path(qc_dir, paste0(file_prefix, ".png")),
    p_pca,
    width = 7,
    height = 7,
    dpi = 300
  )
  
  list(plot = p_pca, pca_df = pca_df)
}

pca_all <- run_pca_plot(
  mat         = vst_mat,
  sample_info = sample_info,
  title_text  = "PCA of mouse RNA-seq samples",
  file_prefix = "PCA_PC1_vs_PC2_selected_samples"
)

# ===============================================
# 9) SAMPLE CORRELATION HEATMAP (genotype + sex)
# ===============================================

sample_cor_mat      <- cor(vst_mat, method = "spearman", use = "pairwise.complete.obs")
sample_short_labels <- short_sample_name(colnames(sample_cor_mat))

ann_full <- data.frame(
  genotype = sample_info$genotype,
  sex      = sample_info$sex,
  row.names= rownames(sample_info)
)

pheatmap(
  sample_cor_mat,
  annotation_col = ann_full,
  annotation_row = ann_full,
  labels_col     = sample_short_labels,
  labels_row     = sample_short_labels,
  fontsize_col   = 9,
  fontsize_row   = 9,
  main           = paste0("Sample-to-sample Spearman correlation (", tumor_label, ")"),
  filename       = file.path(qc_dir, "Sample_correlation_heatmap_selected_samples.png"),
  width          = 9,
  height         = 8
)

# =========================================================
# 10) SEX-STRATIFIED DEG ANALYSIS (WT and KO separately)
# =========================================================

run_sex_DEG <- function(genotype_name, counts_full, sample_info, row_map,
                        out_table_dir, out_heatmap_dir) {
  
  message("Running DESeq2: ", genotype_name, " Female vs Male")
  
  samples_use <- rownames(sample_info)[
    sample_info$genotype == genotype_name & !is.na(sample_info$sex)
  ]
  if (length(samples_use) < 4) {
    warning("Too few samples for ", genotype_name, " to run DESeq2.")
    return(NULL)
  }
  
  counts_sub  <- counts_full[, samples_use, drop = FALSE]
  coldata_sub <- sample_info[samples_use, , drop = FALSE]
  
  coldata_sub$sex <- factor(coldata_sub$sex, levels = c("Male", "Female"))
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_sub),
    colData   = coldata_sub,
    design    = ~ sex
  )
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  
  norm_sub <- counts(dds, normalized = TRUE)
  vst_sub  <- assay(vst(dds, blind = TRUE))
  
  res <- results(dds, contrast = c("sex", "Female", "Male"), alpha = 0.05)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ENSEMBL") %>%
    left_join(row_map %>% select(ENSEMBL, SYMBOL, ENTREZID, GENENAME), by = "ENSEMBL") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  res_df$regulated <- case_when(
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >  0.58 ~ "Up_in_Female",
    !is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.58 ~ "Up_in_Male",
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
    Comparison     = paste0(genotype_name, "_Female_vs_Male"),
    Total_DE_Genes = nrow(deg_sig),
    Up_in_Female   = sum(deg_sig$log2FoldChange >  0.58, na.rm = TRUE),
    Up_in_Male     = sum(deg_sig$log2FoldChange < -0.58, na.rm = TRUE),
    padj_Cutoff    = 0.05,
    log2FC_Cutoff  = 0.58,
    stringsAsFactors = FALSE
  )
  
  write.csv(res_df,
            file.path(out_table_dir, paste0("DESeq2_all_results_Female_vs_Male_", genotype_name, ".csv")),
            row.names = FALSE)
  
  write.csv(deg_sig,
            file.path(out_table_dir, paste0("DESeq2_significant_DEGs_Female_vs_Male_", genotype_name, ".csv")),
            row.names = FALSE)
  
  write.csv(deg_top50,
            file.path(out_table_dir, paste0("DESeq2_top50_ranked_genes_Female_vs_Male_", genotype_name, ".csv")),
            row.names = FALSE)
  
  write.csv(deg_top100,
            file.path(out_table_dir, paste0("DESeq2_top100_ranked_genes_Female_vs_Male_", genotype_name, ".csv")),
            row.names = FALSE)
  
  write.csv(deg_summary,
            file.path(out_table_dir, paste0("DEG_summary_Female_vs_Male_", genotype_name, ".csv")),
            row.names = FALSE)
  
  deg_summary_tbl <- gridExtra::tableGrob(
    deg_summary,
    rows  = NULL,
    theme = gridExtra::ttheme_minimal(base_size = 9)
  )
  
  png(file.path(out_table_dir, paste0("DEG_summary_table_Female_vs_Male_", genotype_name, ".png")),
      width = 3600, height = 800, res = 300)
  grid.newpage(); grid.draw(deg_summary_tbl); dev.off()
  
  plot_deg_heatmap <- function(gene_table, mat, sample_info, filename, title_text, scale_rows = TRUE) {
    genes_use <- intersect(gene_table$ENSEMBL, rownames(mat))
    if (length(genes_use) < 2) return(NULL)
    
    hm_mat <- mat[genes_use, , drop = FALSE]
    gene_symbols <- gene_table$SYMBOL[match(rownames(hm_mat), gene_table$ENSEMBL)]
    
    rownames(hm_mat) <- gene_symbols
    rownames(hm_mat)[is.na(rownames(hm_mat)) | rownames(hm_mat) == ""] <-
      genes_use[is.na(rownames(hm_mat)) | rownames(hm_mat) == ""]
    rownames(hm_mat) <- make.unique(rownames(hm_mat))
    
    ann_col <- data.frame(
      sex      = sample_info[colnames(hm_mat), "sex"],
      genotype = sample_info[colnames(hm_mat), "genotype"],
      row.names = colnames(hm_mat)
    )
    
    pheatmap(
      hm_mat,
      scale = if (scale_rows) "row" else "none",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method        = "complete",
      annotation_col           = ann_col,
      labels_col               = short_sample_name(colnames(hm_mat)),
      show_colnames            = TRUE,
      show_rownames            = TRUE,
      fontsize_row             = ifelse(nrow(hm_mat) > 80, 5, 7),
      fontsize_col             = 9,
      angle_col                = 45,
      main                     = paste0(title_text, " (", genotype_name, ", ", tumor_label, ")"),
      filename = filename,
      width    = 10,
      height   = 12
    )
  }
  
  plot_deg_heatmap(
    deg_sig,
    vst_sub,
    sample_info[colnames(vst_sub), , drop = FALSE],
    file.path(out_heatmap_dir, paste0("Heatmap_all_significant_DEGs_Female_vs_Male_", genotype_name, ".png")),
    "All significant DEGs"
  )
  
  plot_deg_heatmap(
    deg_top50,
    vst_sub,
    sample_info[colnames(vst_sub), , drop = FALSE],
    file.path(out_heatmap_dir, paste0("Heatmap_top50_DEGs_Female_vs_Male_", genotype_name, ".png")),
    "Top 50 DEGs"
  )
  
  plot_deg_heatmap(
    deg_top100,
    vst_sub,
    sample_info[colnames(vst_sub), , drop = FALSE],
    file.path(out_heatmap_dir, paste0("Heatmap_top100_DEGs_Female_vs_Male_", genotype_name, ".png")),
    "Top 100 DEGs"
  )
  
  list(
    dds        = dds,
    norm       = norm_sub,
    vst        = vst_sub,
    res_df     = res_df,
    deg_sig    = deg_sig,
    deg_top50  = deg_top50,
    deg_top100 = deg_top100,
    deg_summary= deg_summary
  )
}

deg_WT <- run_sex_DEG("WT", filtered_counts, sample_info, row_map, deg_table_WT_dir, deg_heatmap_WT_dir)
deg_KO <- run_sex_DEG("KO", filtered_counts, sample_info, row_map, deg_table_KO_dir, deg_heatmap_KO_dir)

# ======================================================
# 11) BIOMARKER PANEL WITH SEX COMPARISON PER GENOTYPE
# ======================================================

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
  file = file.path(outdir, "biomarker_presence_check_selected_samples.csv"),
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
    sample   = colnames(vst_mat),
    genotype = sample_info[colnames(vst_mat), "genotype"],
    sex      = sample_info[colnames(vst_mat), "sex"],
    gene     = g,
    expr     = x,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(
  gene_expr_df,
  file = file.path(outdir, "biomarker_expression_long_VST_selected_samples.csv"),
  row.names = FALSE
)

# =========================================================
# 12) biomarker boxplots, Female vs Male within genotype
# =========================================================

plot_sex_box <- function(df_all, gene_name, genotype_name, outdir, ylab = "VST expression") {
  df_sub <- df_all %>%
    filter(
      gene == gene_name,
      genotype == genotype_name,
      !is.na(sex)
    ) %>%
    group_by(sample, genotype, sex, gene) %>%
    summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop")
  
  message("Gene: ", gene_name, " | Genotype: ", genotype_name,
          " | n rows: ", nrow(df_sub),
          " | n unique samples: ", n_distinct(df_sub$sample))
  
  if (nrow(df_sub) < 4) {
    warning("Too few samples for ", gene_name, " in ", genotype_name, " to test.")
    return(NULL)
  }
  
  df_sub$sex <- factor(df_sub$sex, levels = c("Male", "Female"))
  
  wilcox_res <- tryCatch(wilcox.test(expr ~ sex, data = df_sub), error = function(e) NULL)
  p_label <- if (!is.null(wilcox_res)) {
    paste0("Wilcoxon p = ", signif(wilcox_res$p.value, 3))
  } else {
    "Wilcoxon p = NA"
  }
  
  y_range <- range(df_sub$expr, na.rm = TRUE)
  y_pad   <- ifelse(diff(y_range) == 0, 0.5, 0.15 * diff(y_range))
  y_pos   <- y_range[2] + 0.5 * y_pad
  
  p <- ggplot(df_sub, aes(x = sex, y = expr, fill = sex)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.9, linewidth = 0.9) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.55, color = "black") +
    annotate("text", x = 1.5, y = y_pos, label = p_label, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = sex_cols) +
    labs(
      title = paste0(gene_name, " (", genotype_name, ")"),
      x     = NULL,
      y     = ylab
    ) +
    coord_cartesian(ylim = c(y_range[1], y_range[2] + y_pad)) +
    theme_nature(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(face = "bold", size = 18),
      axis.text  = element_text(size = 15),
      legend.position = "none"
    )
  
  out_file <- file.path(outdir, paste0(gene_name, "_Female_vs_Male_", genotype_name, "_boxplot.png"))
  ggsave(out_file, p, width = 6, height = 5.5, dpi = 300)
  
  if (!is.null(wilcox_res)) {
    data.frame(gene = gene_name, genotype = genotype_name, p.value = wilcox_res$p.value)
  } else {
    data.frame(gene = gene_name, genotype = genotype_name, p.value = NA_real_)
  }
}

sex_test_results <- lapply(available_biomarkers, function(g) {
  rbind(
    plot_sex_box(gene_expr_df, g, "WT", expr_WT_dir),
    plot_sex_box(gene_expr_df, g, "KO", expr_KO_dir)
  )
}) %>% bind_rows()

sex_test_results$p.adj <- p.adjust(sex_test_results$p.value, method = "BH")

write.csv(
  sex_test_results,
  file = file.path(outdir, "Female_vs_Male_wilcox_biomarkers_by_genotype.csv"),
  row.names = FALSE
)

group_summary <- gene_expr_df %>%
  group_by(gene, genotype, sex) %>%
  summarise(
    n           = sum(is.finite(expr)),
    mean_expr   = mean(expr, na.rm = TRUE),
    median_expr = median(expr, na.rm = TRUE),
    sd_expr     = sd(expr, na.rm = TRUE),
    .groups     = "drop"
  )

write.csv(
  group_summary,
  file = file.path(outdir, "biomarker_group_summary_by_genotype_sex.csv"),
  row.names = FALSE
)

# =========================================================
# FOLH1 CORRELATIONS STRATIFIED BY GENOTYPE AND SEX
# =========================================================

correlation_root_dir <- file.path(
  outdir,
  "05_Correlation_results_by_genotype_sex"
)

spearman_dir <- file.path(
  correlation_root_dir,
  "01_Spearman"
)

spearman_table_dir <- file.path(
  spearman_dir,
  "tables"
)

spearman_plot_dir <- file.path(
  spearman_dir,
  "barplots"
)

pearson_dir <- file.path(
  correlation_root_dir,
  "02_Pearson"
)

pearson_table_dir <- file.path(
  pearson_dir,
  "tables"
)

pearson_plot_dir <- file.path(
  pearson_dir,
  "barplots"
)

for (directory in c(
  correlation_root_dir,
  spearman_dir,
  spearman_table_dir,
  spearman_plot_dir,
  pearson_dir,
  pearson_table_dir,
  pearson_plot_dir
)) {
  dir.create(
    directory,
    showWarnings = FALSE,
    recursive = TRUE
  )
}


# ==========================================
# DEFINE THE FOUR GROUPS
# ==========================================

correlation_groups <- list(
  Males_KO = list(
    genotype = "KO",
    sex = "Male"
  ),
  
  Males_WT = list(
    genotype = "WT",
    sex = "Male"
  ),
  
  Females_KO = list(
    genotype = "KO",
    sex = "Female"
  ),
  
  Females_WT = list(
    genotype = "WT",
    sex = "Female"
  )
)


# Check that the expected metadata columns exist
required_sample_info_columns <- c(
  "sample",
  "genotype",
  "sex"
)

missing_sample_info_columns <- setdiff(
  required_sample_info_columns,
  colnames(sample_info)
)

if (length(missing_sample_info_columns) > 0) {
  stop(
    paste(
      "The following columns are missing from sample_info:",
      paste(
        missing_sample_info_columns,
        collapse = ", "
      )
    )
  )
}


# Check the number of samples in each group
correlation_group_sizes <- lapply(
  names(correlation_groups),
  function(group_name) {
    
    group_definition <- correlation_groups[[group_name]]
    
    samples_in_group <- rownames(sample_info)[
      sample_info$genotype == group_definition$genotype &
        sample_info$sex == group_definition$sex
    ]
    
    data.frame(
      correlation_group = group_name,
      genotype = group_definition$genotype,
      sex = group_definition$sex,
      n_samples = length(samples_in_group),
      samples = paste(
        samples_in_group,
        collapse = "; "
      ),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()


write.csv(
  correlation_group_sizes,
  file = file.path(
    correlation_root_dir,
    "correlation_group_sample_sizes.csv"
  ),
  row.names = FALSE
)


# =============================
# CORRELATION CALCULATION FUNCTION
# =============================

calculate_group_correlations <- function(
    group_name,
    method = c("spearman", "pearson")
) {
  
  method <- match.arg(method)
  
  group_definition <- correlation_groups[[group_name]]
  
  samples_in_group <- rownames(sample_info)[
    sample_info$genotype == group_definition$genotype &
      sample_info$sex == group_definition$sex
  ]
  
  message(
    "Calculating ",
    method,
    " correlations for ",
    group_name,
    " | n = ",
    length(samples_in_group)
  )
  
  if (length(samples_in_group) == 0) {
    warning(
      "No samples found for correlation group: ",
      group_name
    )
    
    return(
      data.frame(
        gene = character(0),
        correlation_group = character(0),
        genotype = character(0),
        sex = character(0),
        method = character(0),
        cor = numeric(0),
        p.value = numeric(0),
        n = integer(0),
        p.adj = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  folh1_rows <- get_symbol_rows("FOLH1")
  
  if (length(folh1_rows) == 0) {
    stop(
      "FOLH1 was not found in the VST expression matrix."
    )
  }
  
  folh1_expr_group <- as.numeric(
    colMeans(
      vst_mat[
        folh1_rows,
        samples_in_group,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  )
  
  correlation_results <- lapply(
    available_biomarkers[
      available_biomarkers != "FOLH1"
    ],
    function(gene_name) {
      
      gene_rows <- get_symbol_rows(
        gene_name
      )
      
      if (length(gene_rows) == 0) {
        return(
          data.frame(
            gene = gene_name,
            correlation_group = group_name,
            genotype = group_definition$genotype,
            sex = group_definition$sex,
            method = method,
            cor = NA_real_,
            p.value = NA_real_,
            n = 0,
            stringsAsFactors = FALSE
          )
        )
      }
      
      biomarker_expr_group <- as.numeric(
        colMeans(
          vst_mat[
            gene_rows,
            samples_in_group,
            drop = FALSE
          ],
          na.rm = TRUE
        )
      )
      
      valid_values <- is.finite(
        folh1_expr_group
      ) &
        is.finite(
          biomarker_expr_group
        )
      
      n_valid <- sum(
        valid_values
      )
      
      if (n_valid < 3) {
        return(
          data.frame(
            gene = gene_name,
            correlation_group = group_name,
            genotype = group_definition$genotype,
            sex = group_definition$sex,
            method = method,
            cor = NA_real_,
            p.value = NA_real_,
            n = n_valid,
            stringsAsFactors = FALSE
          )
        )
      }
      
      folh1_values <- folh1_expr_group[
        valid_values
      ]
      
      biomarker_values <- biomarker_expr_group[
        valid_values
      ]
      
      # Correlation cannot be calculated if either vector is constant
      if (
        sd(folh1_values) == 0 ||
        sd(biomarker_values) == 0
      ) {
        return(
          data.frame(
            gene = gene_name,
            correlation_group = group_name,
            genotype = group_definition$genotype,
            sex = group_definition$sex,
            method = method,
            cor = NA_real_,
            p.value = NA_real_,
            n = n_valid,
            stringsAsFactors = FALSE
          )
        )
      }
      
      correlation_test <- if (
        method == "spearman"
      ) {
        suppressWarnings(
          cor.test(
            x = biomarker_values,
            y = folh1_values,
            method = "spearman",
            exact = FALSE
          )
        )
      } else {
        suppressWarnings(
          cor.test(
            x = biomarker_values,
            y = folh1_values,
            method = "pearson"
          )
        )
      }
      
      data.frame(
        gene = gene_name,
        correlation_group = group_name,
        genotype = group_definition$genotype,
        sex = group_definition$sex,
        method = method,
        cor = unname(
          correlation_test$estimate
        ),
        p.value = correlation_test$p.value,
        n = n_valid,
        stringsAsFactors = FALSE
      )
    }
  ) %>%
    bind_rows()
  
  correlation_results$p.adj <- p.adjust(
    correlation_results$p.value,
    method = "BH"
  )
  
  correlation_results %>%
    arrange(
      desc(cor)
    )
}


# =============================
# CALCULATE SPEARMAN RESULTS
# =============================

spearman_results_by_group <- lapply(
  names(correlation_groups),
  function(group_name) {
    
    calculate_group_correlations(
      group_name = group_name,
      method = "spearman"
    )
  }
)

names(spearman_results_by_group) <- names(
  correlation_groups
)

spearman_results_all <- bind_rows(
  spearman_results_by_group
)


write.csv(
  spearman_results_all,
  file = file.path(
    spearman_table_dir,
    "FOLH1_biomarker_Spearman_correlations_all_four_groups.csv"
  ),
  row.names = FALSE
)


for (group_name in names(
  spearman_results_by_group
)) {
  
  write.csv(
    spearman_results_by_group[[group_name]],
    file = file.path(
      spearman_table_dir,
      paste0(
        "FOLH1_biomarker_Spearman_",
        group_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


# =============================
# CALCULATE PEARSON RESULTS
# =============================

pearson_results_by_group <- lapply(
  names(correlation_groups),
  function(group_name) {
    
    calculate_group_correlations(
      group_name = group_name,
      method = "pearson"
    )
  }
)

names(pearson_results_by_group) <- names(
  correlation_groups
)

pearson_results_all <- bind_rows(
  pearson_results_by_group
)


write.csv(
  pearson_results_all,
  file = file.path(
    pearson_table_dir,
    "FOLH1_biomarker_Pearson_correlations_all_four_groups.csv"
  ),
  row.names = FALSE
)


for (group_name in names(
  pearson_results_by_group
)) {
  
  write.csv(
    pearson_results_by_group[[group_name]],
    file = file.path(
      pearson_table_dir,
      paste0(
        "FOLH1_biomarker_Pearson_",
        group_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}


# =========================================
# CORRELATION BARPLOT FUNCTION
# =========================================

plot_group_correlation_bar <- function(
    correlation_df,
    group_name,
    method,
    output_dir
) {
  
  plot_df <- correlation_df %>%
    filter(
      is.finite(cor)
    ) %>%
    mutate(
      direction = ifelse(
        cor >= 0,
        "Positive",
        "Negative"
      ),
      gene = factor(
        gene,
        levels = gene[
          order(cor)
        ]
      )
    )
  
  if (nrow(plot_df) == 0) {
    warning(
      "No valid ",
      method,
      " correlations available for ",
      group_name
    )
    
    return(NULL)
  }
  
  correlation_bar_colors <- c(
    "Positive" = "#E64B35FF",
    "Negative" = "#4DBBD5FF"
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = cor,
      y = gene,
      fill = direction
    )
  ) +
    geom_col(
      width = 0.76,
      colour = "white",
      linewidth = 0.5
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 1,
      colour = "black"
    ) +
    scale_fill_manual(
      values = correlation_bar_colors
    ) +
    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(
        -1,
        1,
        by = 0.25
      )
    ) +
    labs(
      title = paste0(
        method,
        " correlation of FOLH1 with biomarkers - ",
        group_name
      ),
      x = paste0(
        method,
        " correlation with FOLH1"
      ),
      y = NULL
    ) +
    theme_classic(
      base_size = 18,
      base_family = "Arial"
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 20
      ),
      axis.title = element_text(
        face = "bold",
        colour = "black",
        size = 18
      ),
      axis.text.x = element_text(
        colour = "black",
        size = 14
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 15
      ),
      axis.line = element_line(
        linewidth = 0.9,
        colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.8,
        colour = "black"
      ),
      legend.position = "none"
    )
  
  output_file <- file.path(
    output_dir,
    paste0(
      "FOLH1_biomarker_",
      method,
      "_barplot_",
      group_name,
      ".png"
    )
  )
  
  ggsave(
    output_file,
    p,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  invisible(p)
}


# =========================================
# GENERATE SPEARMAN BARPLOTS
# =========================================

for (group_name in names(
  spearman_results_by_group
)) {
  
  plot_group_correlation_bar(
    correlation_df = spearman_results_by_group[[group_name]],
    group_name = group_name,
    method = "Spearman",
    output_dir = spearman_plot_dir
  )
}


# =========================================
# GENERATE PEARSON BARPLOTS
# =========================================

for (group_name in names(
  pearson_results_by_group
)) {
  
  plot_group_correlation_bar(
    correlation_df = pearson_results_by_group[[group_name]],
    group_name = group_name,
    method = "Pearson",
    output_dir = pearson_plot_dir
  )
}


# =========================================
# COMBINED CORRELATION TABLE
# =========================================

all_correlation_results <- bind_rows(
  spearman_results_all,
  pearson_results_all
)


write.csv(
  all_correlation_results,
  file = file.path(
    correlation_root_dir,
    "FOLH1_biomarker_Spearman_and_Pearson_all_four_groups.csv"
  ),
  row.names = FALSE
)

# =========================================================================
# 15) COMBINED SCATTERPLOTS FOR EACH BIOMARKER (genotype-coloured)
# =========================================================================

scatter_fun_group <- function(g, save_dir = scatter_root_dir) {
  samples_use <- rownames(sample_info)
  
  folh1_rows <- get_symbol_rows("FOLH1")
  folh1_expr <- as.numeric(colMeans(vst_mat[folh1_rows, samples_use, drop = FALSE], na.rm = TRUE))
  
  rows          <- get_symbol_rows(g)
  biomarker_expr<- as.numeric(colMeans(vst_mat[rows, samples_use, drop = FALSE], na.rm = TRUE))
  
  df <- data.frame(
    sample    = samples_use,
    FOLH1     = folh1_expr,
    biomarker = biomarker_expr,
    genotype  = sample_info[samples_use, "genotype"],
    sex       = sample_info[samples_use, "sex"],
    stringsAsFactors = FALSE
  )
  
  df <- df[is.finite(df$FOLH1) & is.finite(df$biomarker), , drop = FALSE]
  if (nrow(df) < 3) return(NULL)
  
  ct_all <- suppressWarnings(cor.test(df$FOLH1, df$biomarker, method = "spearman", exact = FALSE))
  lab_all <- paste0(
    "All samples: Spearman rho = ", round(unname(ct_all$estimate), 3),
    ", p = ", signif(ct_all$p.value, 3)
  )
  
  p <- ggplot(df, aes(x = FOLH1, y = biomarker, color = genotype)) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_smooth(aes(group = genotype, fill = genotype), method = "lm", se = FALSE, linewidth = 0.85) +
    scale_color_manual(values = group_cols_box) +
    scale_fill_manual(values = group_cols_box) +
    labs(
      title    = paste0("FOLH1 vs ", g, " (", tumor_label, ")"),
      subtitle = lab_all,
      x        = "FOLH1 (VST)",
      y        = paste0(g, " (VST)")
    ) +
    theme_nature(base_size = 15) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle= element_text(size = 12, hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  ggsave(
    file.path(save_dir, paste0("FOLH1_vs_", g, "_scatter_combined_WT_KO_selected_samples.png")),
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

# ========================================================================
# 16) GO BP ORA + GSEA: Female vs Male within WT and within KO
# ========================================================================

run_pathway_for_sex_comparison <- function(deg_obj, genotype_label) {
  
  if (is.null(deg_obj)) return(NULL)
  
  res_df <- deg_obj$res_df
  deg_sig<- deg_obj$deg_sig
  
  # Ranked gene list (log2FC: Female vs Male)
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
  
  # sets up in Female and up in Male
  sig_entrez_up_F <- deg_sig %>%
    filter(log2FoldChange >  0.58, !is.na(ENTREZID)) %>%
    pull(ENTREZID) %>% unique()
  
  sig_entrez_up_M <- deg_sig %>%
    filter(log2FoldChange < -0.58, !is.na(ENTREZID)) %>%
    pull(ENTREZID) %>% unique()
  
  ego_up_F <- if (length(sig_entrez_up_F) > 0) {
    enrichGO(
      gene          = sig_entrez_up_F,
      universe      = universe_entrez,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  } else NULL
  
  ego_up_M <- if (length(sig_entrez_up_M) > 0) {
    enrichGO(
      gene          = sig_entrez_up_M,
      universe      = universe_entrez,
      OrgDb         = org.Mm.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  } else NULL
  
  # GO BP ORA plots
  if (!is.null(ego_up_F) && nrow(as.data.frame(ego_up_F)) > 0) {
    csv_name  <- paste0("GO_BP_ORA_up_in_Female_vs_Male_", genotype_label, ".csv")
    plot_name <- paste0("GO_BP_ORA_up_in_Female_vs_Male_", genotype_label, "_dotplot.png")
    
    write.csv(as.data.frame(ego_up_F), file.path(go_bp_dir, csv_name), row.names = FALSE)
    
    png(file.path(go_bp_dir, plot_name), width = 2400, height = 1800, res = 300)
    print(
      dotplot(ego_up_F, showCategory = 20) +
        ggtitle(paste0("GO BP enriched in Female vs Male (", genotype_label, ")")) +
        theme_enrich_medium()
    )
    dev.off()
  }
  
  if (!is.null(ego_up_M) && nrow(as.data.frame(ego_up_M)) > 0) {
    csv_name  <- paste0("GO_BP_ORA_up_in_Male_vs_Female_", genotype_label, ".csv")
    plot_name <- paste0("GO_BP_ORA_up_in_Male_vs_Female_", genotype_label, "_dotplot.png")
    
    write.csv(as.data.frame(ego_up_M), file.path(go_bp_dir, csv_name), row.names = FALSE)
    
    png(file.path(go_bp_dir, plot_name), width = 2400, height = 1800, res = 300)
    print(
      dotplot(ego_up_M, showCategory = 20) +
        ggtitle(paste0("GO BP enriched in Male vs Female (", genotype_label, ")")) +
        theme_enrich_medium()
    )
    dev.off()
  }
  
  # MSigDB sets (Hallmark + Immune)
  msig_hallmark_mm <- msigdbr(species = "Mus musculus", category = "H")
  hallmark_term2gene <- msig_hallmark_mm %>% select(gs_name, entrez_gene)
  
  msig_immune_mm <- msigdbr(species = "Mus musculus", category = "C7")
  immune_term2gene <- msig_immune_mm %>% select(gs_name, entrez_gene)
  
  gsea_hallmark <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = hallmark_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    verbose       = FALSE
  )
  
  if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
    csv_name  <- paste0("GSEA_Hallmark_ranked_Female_vs_Male_", genotype_label, ".csv")
    plot_name <- paste0("GSEA_Hallmark_dotplot_Female_vs_Male_", genotype_label, ".png")
    
    write.csv(as.data.frame(gsea_hallmark), file.path(hallmark_gsea_dir, csv_name), row.names = FALSE)
    
    png(file.path(hallmark_gsea_dir, plot_name), width = 2400, height = 1800, res = 300)
    print(
      dotplot(gsea_hallmark, showCategory = 20, split = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle(paste0("Hallmark GSEA (Female vs Male, ", genotype_label, ")")) +
        theme_enrich_small()
    )
    dev.off()
  }
  
  gsea_immune <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = immune_term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    verbose       = FALSE
  )
  
  if (!is.null(gsea_immune) && nrow(as.data.frame(gsea_immune)) > 0) {
    csv_name  <- paste0("GSEA_ImmuneSigDB_ranked_Female_vs_Male_", genotype_label, ".csv")
    plot_name <- paste0("GSEA_ImmuneSigDB_dotplot_Female_vs_Male_", genotype_label, ".png")
    
    write.csv(as.data.frame(gsea_immune), file.path(immune_gsea_dir, csv_name), row.names = FALSE)
    
    png(file.path(immune_gsea_dir, plot_name), width = 2400, height = 1800, res = 300)
    print(
      dotplot(gsea_immune, showCategory = 20, split = ".sign") +
        facet_grid(. ~ .sign) +
        ggtitle(paste0("ImmuneSigDB GSEA (Female vs Male, ", genotype_label, ")")) +
        theme_enrich_small()
    )
    dev.off()
  }
  
  list(
    ego_up_F     = ego_up_F,
    ego_up_M     = ego_up_M,
    gsea_hallmark= gsea_hallmark,
    gsea_immune  = gsea_immune
  )
}

pathway_WT <- run_pathway_for_sex_comparison(deg_WT, "WT")
pathway_KO <- run_pathway_for_sex_comparison(deg_KO, "KO")

# =========================================================================
# mMCP-Counter: four genotype-by-sex comparisons
# =========================================================================
#
# Comparisons:
# 1. WT_Female vs WT_Male
# 2. KO_Female vs KO_Male
# 3. WT_Male vs KO_Male
# 4. WT_Female vs KO_Female
#
# Required objects created earlier in the script:
# norm_mat
# row_map
# sample_info
# outdir
#
# Required sample_info columns:
# sample
# genotype
# sex
# =========================================================================


# =========================
# mMCP-Counter DIRECTORIES
# =========================

mmcp_root_dir <- file.path(
  outdir,
  "07_mMCP_counter"
)

mmcp_input_dir <- file.path(
  mmcp_root_dir,
  "01_input"
)

mmcp_scores_dir <- file.path(
  mmcp_root_dir,
  "02_scores"
)

mmcp_group_dir <- file.path(
  mmcp_root_dir,
  "03_group_comparisons"
)

mmcp_plot_dir <- file.path(
  mmcp_root_dir,
  "04_plots"
)

mmcp_comparison_plot_dir <- file.path(
  mmcp_plot_dir,
  "four_comparisons"
)

mmcp_comparison_table_dir <- file.path(
  mmcp_group_dir,
  "four_comparisons"
)

for (directory in c(
  mmcp_root_dir,
  mmcp_input_dir,
  mmcp_scores_dir,
  mmcp_group_dir,
  mmcp_plot_dir,
  mmcp_comparison_plot_dir,
  mmcp_comparison_table_dir
)) {
  dir.create(
    directory,
    showWarnings = FALSE,
    recursive = TRUE
  )
}


# =========================
# PREPARE mMCP INPUT MATRIX
# =========================

mmcp_map <- row_map %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  select(
    ENSEMBL,
    SYMBOL
  ) %>%
  distinct()


mmcp_expr <- data.frame(
  ENSEMBL = rownames(norm_mat),
  norm_mat,
  check.names = FALSE
) %>%
  left_join(
    mmcp_map,
    by = "ENSEMBL"
  ) %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  select(
    -ENSEMBL
  ) %>%
  group_by(SYMBOL) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )


mmcp_mat <- as.data.frame(mmcp_expr)

rownames(mmcp_mat) <- mmcp_mat$SYMBOL

mmcp_mat$SYMBOL <- NULL

mmcp_mat <- as.matrix(mmcp_mat)


write.csv(
  data.frame(
    gene_symbol = rownames(mmcp_mat),
    mmcp_mat,
    check.names = FALSE
  ),
  file = file.path(
    mmcp_input_dir,
    "mMCP_counter_input_symbol_matrix.csv"
  ),
  row.names = FALSE
)


# =========================
# RUN mMCP-COUNTER
# =========================

if (!"mMCPcounter" %in% rownames(installed.packages())) {
  stop(
    paste0(
      "Package 'mMCPcounter' is not installed. ",
      "Run: remotes::install_github(",
      "'cit-bioinfo/mMCP-counter')"
    )
  )
}


mmcp_res <- mMCPcounter.estimate(
  exp = mmcp_mat,
  features = "Gene.Symbol",
  genomeVersion = "GCRm39"
)


mmcp_scores <- as.data.frame(mmcp_res)

mmcp_scores$cell_type <- rownames(mmcp_scores)


# =========================
# LONG-FORMAT mMCP RESULTS
# =========================

mmcp_scores_long <- mmcp_scores %>%
  pivot_longer(
    cols = -cell_type,
    names_to = "sample",
    values_to = "score"
  ) %>%
  left_join(
    sample_info %>%
      select(
        sample,
        genotype,
        sex
      ),
    by = "sample"
  ) %>%
  mutate(
    group4 = paste(
      genotype,
      sex,
      sep = "_"
    )
  )


write.csv(
  mmcp_scores_long,
  file = file.path(
    mmcp_scores_dir,
    "mMCP_counter_scores_long_with_genotype_sex.csv"
  ),
  row.names = FALSE
)


mmcp_wide <- mmcp_scores_long %>%
  select(
    cell_type,
    sample,
    score
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = score
  )


write.csv(
  mmcp_wide,
  file = file.path(
    mmcp_scores_dir,
    "mMCP_counter_scores_wide.csv"
  ),
  row.names = FALSE
)


# =========================
# CHECK GROUP LABELS
# =========================

expected_groups <- c(
  "WT_Female",
  "WT_Male",
  "KO_Female",
  "KO_Male"
)

observed_groups <- sort(
  unique(
    na.omit(mmcp_scores_long$group4)
  )
)

missing_groups <- setdiff(
  expected_groups,
  observed_groups
)

if (length(missing_groups) > 0) {
  warning(
    paste(
      "The following expected groups were not found:",
      paste(missing_groups, collapse = ", ")
    )
  )
}


# ==========================================
# DEFINE THE FOUR REQUESTED COMPARISONS
# ==========================================

comparison_definitions <- list(
  WT_Female_vs_WT_Male = c(
    "WT_Female",
    "WT_Male"
  ),
  
  KO_Female_vs_KO_Male = c(
    "KO_Female",
    "KO_Male"
  ),
  
  WT_Male_vs_KO_Male = c(
    "WT_Male",
    "KO_Male"
  ),
  
  WT_Female_vs_KO_Female = c(
    "WT_Female",
    "KO_Female"
  )
)


comparison_labels <- c(
  WT_Female_vs_WT_Male = "WT female vs WT male",
  KO_Female_vs_KO_Male = "KO female vs KO male",
  WT_Male_vs_KO_Male = "WT male vs KO male",
  WT_Female_vs_KO_Female = "WT female vs KO female"
)


# Colours for the four genotype-by-sex groups
mMCP_group4_cols <- c(
  "WT_Female" = "#E31A1C",
  "WT_Male" = "#1F78B4",
  "KO_Female" = "#FB6A4A",
  "KO_Male" = "#6A3D9A"
)


# ==========================================
# GROUP SUMMARY FOR ALL FOUR GROUPS
# ==========================================

mmcp_group_summary <- mmcp_scores_long %>%
  group_by(
    cell_type,
    group4,
    genotype,
    sex
  ) %>%
  summarise(
    n = sum(is.finite(score)),
    mean_score = mean(
      score,
      na.rm = TRUE
    ),
    median_score = median(
      score,
      na.rm = TRUE
    ),
    sd_score = sd(
      score,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


write.csv(
  mmcp_group_summary,
  file = file.path(
    mmcp_group_dir,
    "mMCP_group_summary_by_genotype_sex.csv"
  ),
  row.names = FALSE
)


# ==========================================
# WILCOXON TEST FUNCTION
# ==========================================

run_mmcp_pairwise_test <- function(
    data,
    comparison_name,
    group_pair
) {
  
  comparison_data <- data %>%
    filter(
      group4 %in% group_pair,
      is.finite(score)
    ) %>%
    mutate(
      group4 = factor(
        group4,
        levels = group_pair
      )
    )
  
  results <- comparison_data %>%
    group_by(cell_type) %>%
    group_modify(
      ~ {
        
        test_data <- .x %>%
          filter(
            is.finite(score)
          )
        
        group_counts <- table(
          test_data$group4
        )
        
        if (
          length(group_counts) < 2 ||
          any(group_counts < 1)
        ) {
          return(
            data.frame(
              n_group_1 = ifelse(
                length(group_counts) >= 1,
                as.numeric(group_counts[1]),
                0
              ),
              n_group_2 = ifelse(
                length(group_counts) >= 2,
                as.numeric(group_counts[2]),
                0
              ),
              statistic = NA_real_,
              p.value = NA_real_
            )
          )
        }
        
        wilcox_result <- tryCatch(
          wilcox.test(
            score ~ group4,
            data = test_data,
            exact = FALSE
          ),
          error = function(e) NULL
        )
        
        if (is.null(wilcox_result)) {
          return(
            data.frame(
              n_group_1 = as.numeric(group_counts[1]),
              n_group_2 = as.numeric(group_counts[2]),
              statistic = NA_real_,
              p.value = NA_real_
            )
          )
        }
        
        data.frame(
          n_group_1 = as.numeric(group_counts[1]),
          n_group_2 = as.numeric(group_counts[2]),
          statistic = unname(
            wilcox_result$statistic
          ),
          p.value = wilcox_result$p.value
        )
      }
    ) %>%
    ungroup() %>%
    mutate(
      comparison = comparison_name,
      comparison_label = comparison_labels[
        comparison_name
      ],
      group_1 = group_pair[1],
      group_2 = group_pair[2],
      p.adj = p.adjust(
        p.value,
        method = "BH"
      )
    ) %>%
    select(
      comparison,
      comparison_label,
      cell_type,
      group_1,
      group_2,
      n_group_1,
      n_group_2,
      statistic,
      p.value,
      p.adj
    )
  
  return(results)
}


# ==========================================
# RUN ALL FOUR COMPARISONS
# ==========================================

mmcp_pairwise_results <- purrr::map2_dfr(
  names(comparison_definitions),
  comparison_definitions,
  ~ run_mmcp_pairwise_test(
    data = mmcp_scores_long,
    comparison_name = .x,
    group_pair = .y
  )
)


write.csv(
  mmcp_pairwise_results,
  file = file.path(
    mmcp_group_dir,
    "mMCP_pairwise_tests_four_comparisons.csv"
  ),
  row.names = FALSE
)


# ==========================================
# mMCP BOXPLOT FUNCTION
# ==========================================

plot_mmcp_comparison <- function(
    cell_type_name,
    comparison_name,
    group_pair,
    data
) {
  
  plot_df <- data %>%
    filter(
      cell_type == cell_type_name,
      group4 %in% group_pair,
      is.finite(score)
    ) %>%
    mutate(
      group4 = factor(
        group4,
        levels = group_pair
      )
    )
  
  if (nrow(plot_df) == 0) {
    return(NULL)
  }
  
  group_counts <- table(
    plot_df$group4
  )
  
  if (
    length(group_counts) < 2 ||
    any(group_counts < 1)
  ) {
    p_label <- "Wilcoxon p = NA"
  } else {
    
    wilcox_result <- tryCatch(
      wilcox.test(
        score ~ group4,
        data = plot_df,
        exact = FALSE
      ),
      error = function(e) NULL
    )
    
    p_label <- if (!is.null(wilcox_result)) {
      paste0(
        "Wilcoxon p = ",
        signif(
          wilcox_result$p.value,
          3
        )
      )
    } else {
      "Wilcoxon p = NA"
    }
  }
  
  y_range <- range(
    plot_df$score,
    na.rm = TRUE
  )
  
  y_pad <- ifelse(
    diff(y_range) == 0,
    0.5,
    0.15 * diff(y_range)
  )
  
  y_pos <- y_range[2] + 0.45 * y_pad
  
  p <- ggplot(
    plot_df,
    aes(
      x = group4,
      y = score,
      fill = group4
    )
  ) +
    geom_boxplot(
      width = 0.62,
      outlier.shape = NA,
      alpha = 0.9,
      linewidth = 0.85
    ) +
    geom_jitter(
      width = 0.12,
      size = 2.0,
      alpha = 0.65,
      colour = "black"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_pos,
      label = p_label,
      size = 4.2,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = mMCP_group4_cols,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "mMCP-counter: ",
        cell_type_name
      ),
      subtitle = paste0(
        comparison_labels[
          comparison_name
        ],
        " (",
        tumor_label,
        ")"
      ),
      x = NULL,
      y = "mMCP-counter score"
    ) +
    coord_cartesian(
      ylim = c(
        y_range[1],
        y_range[2] + y_pad
      )
    ) +
    theme_nature(
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 16
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 12
      ),
      axis.text.x = element_text(
        angle = 30,
        hjust = 1
      ),
      legend.position = "none"
    )
  
  safe_cell_type <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    cell_type_name
  )
  
  output_file <- file.path(
    mmcp_comparison_plot_dir,
    paste0(
      "mMCP_boxplot_",
      safe_cell_type,
      "_",
      comparison_name,
      ".png"
    )
  )
  
  ggsave(
    output_file,
    p,
    width = 7.4,
    height = 5.8,
    dpi = 300
  )
  
  invisible(p)
}


# ==========================================
# GENERATE BOXPLOTS FOR ALL COMPARISONS
# ==========================================

mmcp_cell_types <- sort(
  unique(
    mmcp_scores_long$cell_type
  )
)

for (comparison_name in names(comparison_definitions)) {
  
  group_pair <- comparison_definitions[[comparison_name]]
  
  comparison_results <- mmcp_pairwise_results %>%
    filter(
      comparison == comparison_name
    )
  
  write.csv(
    comparison_results,
    file = file.path(
      mmcp_comparison_table_dir,
      paste0(
        "mMCP_pairwise_",
        comparison_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  for (cell_type_name in mmcp_cell_types) {
    
    plot_mmcp_comparison(
      cell_type_name = cell_type_name,
      comparison_name = comparison_name,
      group_pair = group_pair,
      data = mmcp_scores_long
    )
  }
}


# ==========================================
# HEATMAP ANNOTATED BY GENOTYPE AND SEX
# ==========================================

mmcp_heat <- mmcp_scores_long %>%
  select(
    cell_type,
    sample,
    score
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = score
  ) %>%
  as.data.frame()

rownames(mmcp_heat) <- mmcp_heat$cell_type

mmcp_heat$cell_type <- NULL

mmcp_heat <- as.matrix(mmcp_heat)


mmcp_ann_col <- sample_info[
  colnames(mmcp_heat),
  c(
    "genotype",
    "sex"
  ),
  drop = FALSE
]

mmcp_ann_col$group4 <- paste(
  mmcp_ann_col$genotype,
  mmcp_ann_col$sex,
  sep = "_"
)


pheatmap(
  mmcp_heat,
  scale = "row",
  annotation_col = mmcp_ann_col,
  main = paste0(
    "mMCP-counter scores (",
    tumor_label,
    ")"
  ),
  fontsize_col = 8,
  fontsize_row = 8,
  filename = file.path(
    mmcp_plot_dir,
    "mMCP_heatmap_all_samples_genotype_sex.png"
  ),
  width = 10,
  height = 8
)


# ==========================================
# PCA ANNOTATED BY GENOTYPE AND SEX
# ==========================================

mmcp_pca_obj <- prcomp(
  t(mmcp_heat),
  scale. = TRUE
)

mmcp_pct <- (
  mmcp_pca_obj$sdev^2 /
    sum(mmcp_pca_obj$sdev^2)
) * 100


mmcp_pca <- data.frame(
  sample = rownames(
    mmcp_pca_obj$x
  ),
  PC1 = mmcp_pca_obj$x[, 1],
  PC2 = mmcp_pca_obj$x[, 2],
  genotype = sample_info[
    rownames(mmcp_pca_obj$x),
    "genotype"
  ],
  sex = sample_info[
    rownames(mmcp_pca_obj$x),
    "sex"
  ],
  stringsAsFactors = FALSE
)

mmcp_pca$group4 <- paste(
  mmcp_pca$genotype,
  mmcp_pca$sex,
  sep = "_"
)


p_mmcp_pca <- ggplot(
  mmcp_pca,
  aes(
    PC1,
    PC2,
    color = group4
  )
) +
  geom_point(
    size = 4
  ) +
  scale_color_manual(
    values = mMCP_group4_cols,
    drop = FALSE
  ) +
  labs(
    title = paste0(
      "PCA of mMCP-counter scores (",
      tumor_label,
      ")"
    ),
    x = paste0(
      "PC1 (",
      round(mmcp_pct[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(mmcp_pct[2], 1),
      "%)"
    ),
    color = "Genotype × sex"
  ) +
  theme_nature(
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right"
  )


ggsave(
  file.path(
    mmcp_plot_dir,
    "mMCP_PCA_genotype_sex.png"
  ),
  p_mmcp_pca,
  width = 7,
  height = 6,
  dpi = 300
)
# ===================
# 18) EXCEL WORKBOOK
# ===================

wb <- createWorkbook()

addWorksheet(wb, "Sample_info")
writeData(wb, "Sample_info", sample_info)

addWorksheet(wb, "PCA_QC_all")
writeData(wb, "PCA_QC_all", pca_all$pca_df)

if (!is.null(deg_WT)) {
  addWorksheet(wb, "DEG_summary_WT")
  writeData(wb, "DEG_summary_WT", deg_WT$deg_summary)
  
  addWorksheet(wb, "DESeq2_WT_all_results")
  writeData(wb, "DESeq2_WT_all_results", deg_WT$res_df)
  
  addWorksheet(wb, "DESeq2_WT_sig_DEGs")
  writeData(wb, "DESeq2_WT_sig_DEGs", deg_WT$deg_sig)
  
  addWorksheet(wb, "DESeq2_WT_top50")
  writeData(wb, "DESeq2_WT_top50", deg_WT$deg_top50)
  
  addWorksheet(wb, "DESeq2_WT_top100")
  writeData(wb, "DESeq2_WT_top100", deg_WT$deg_top100)
}

if (!is.null(deg_KO)) {
  addWorksheet(wb, "DEG_summary_KO")
  writeData(wb, "DEG_summary_KO", deg_KO$deg_summary)
  
  addWorksheet(wb, "DESeq2_KO_all_results")
  writeData(wb, "DESeq2_KO_all_results", deg_KO$res_df)
  
  addWorksheet(wb, "DESeq2_KO_sig_DEGs")
  writeData(wb, "DESeq2_KO_sig_DEGs", deg_KO$deg_sig)
  
  addWorksheet(wb, "DESeq2_KO_top50")
  writeData(wb, "DESeq2_KO_top50", deg_KO$deg_top50)
  
  addWorksheet(wb, "DESeq2_KO_top100")
  writeData(wb, "DESeq2_KO_top100", deg_KO$deg_top100)
}

addWorksheet(wb, "Biomarker_presence")
writeData(wb, "Biomarker_presence", biomarker_presence)

addWorksheet(wb, "Biomarker_group_summary")
writeData(wb, "Biomarker_group_summary", group_summary)

addWorksheet(wb, "Sex_wilcox_results")
writeData(wb, "Sex_wilcox_results", sex_test_results)

addWorksheet(wb, "Spearman_FOLH1")
writeData(wb, "Spearman_FOLH1", cor_results_all)

addWorksheet(wb, "Pearson_FOLH1")
writeData(wb, "Pearson_FOLH1", pearson_results_all)

addWorksheet(wb, "Expression_long")
writeData(wb, "Expression_long", gene_expr_df)

if (!is.null(pathway_WT)) {
  if (!is.null(pathway_WT$ego_up_F)) {
    addWorksheet(wb, "GO_BP_up_FvM_WT")
    writeData(wb, "GO_BP_up_FvM_WT", as.data.frame(pathway_WT$ego_up_F))
  }
  if (!is.null(pathway_WT$ego_up_M)) {
    addWorksheet(wb, "GO_BP_up_MvF_WT")
    writeData(wb, "GO_BP_up_MvF_WT", as.data.frame(pathway_WT$ego_up_M))
  }
  if (!is.null(pathway_WT$gsea_hallmark)) {
    addWorksheet(wb, "GSEA_Hallmark_WT")
    writeData(wb, "GSEA_Hallmark_WT", as.data.frame(pathway_WT$gsea_hallmark))
  }
  if (!is.null(pathway_WT$gsea_immune)) {
    addWorksheet(wb, "GSEA_ImmuneSigDB_WT")
    writeData(wb, "GSEA_ImmuneSigDB_WT", as.data.frame(pathway_WT$gsea_immune))
  }
}

if (!is.null(pathway_KO)) {
  if (!is.null(pathway_KO$ego_up_F)) {
    addWorksheet(wb, "GO_BP_up_FvM_KO")
    writeData(wb, "GO_BP_up_FvM_KO", as.data.frame(pathway_KO$ego_up_F))
  }
  if (!is.null(pathway_KO$ego_up_M)) {
    addWorksheet(wb, "GO_BP_up_MvF_KO")
    writeData(wb, "GO_BP_up_MvF_KO", as.data.frame(pathway_KO$ego_up_M))
  }
  if (!is.null(pathway_KO$gsea_hallmark)) {
    addWorksheet(wb, "GSEA_Hallmark_KO")
    writeData(wb, "GSEA_Hallmark_KO", as.data.frame(pathway_KO$gsea_hallmark))
  }
  if (!is.null(pathway_KO$gsea_immune)) {
    addWorksheet(wb, "GSEA_ImmuneSigDB_KO")
    writeData(wb, "GSEA_ImmuneSigDB_KO", as.data.frame(pathway_KO$gsea_immune))
  }
}

addWorksheet(wb, "mMCP_scores_long")
writeData(wb, "mMCP_scores_long", mmcp_scores_long)

addWorksheet(wb, "mMCP_scores_wide")
writeData(wb, "mMCP_scores_wide", mmcp_wide)

addWorksheet(wb, "mMCP_group_summary")
writeData(wb, "mMCP_group_summary", mmcp_group_summary)


saveWorkbook(
  wb,
  file = file.path(outdir, "Overall_mouse_M_vs_F_within_genotype_results.xlsx"),
  overwrite = TRUE
)

message("Analysis complete. Results written to: ", outdir)