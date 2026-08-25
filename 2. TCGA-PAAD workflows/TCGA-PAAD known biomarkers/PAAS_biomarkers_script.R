# =====================================================================
# Ravindu: Comparison of FOLH1 with Known Prostate Cancer Biomarkers in PAAD
# =====================================================================

library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(survival)
library(survminer)
library(pheatmap)

outdir <- "output_PAAD_FOLH1_biomarkers_stage_DFS"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

km_dir <- file.path(outdir, "KM plots")
dir.create(km_dir, showWarnings = FALSE, recursive = TRUE)

stage_dir <- file.path(outdir, "Stage boxplots")
dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)

theme_nature <- function(base_size = 14, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks = element_line(linewidth = 0.6, colour = "black"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(face = "bold")
    )
}
theme_set(theme_nature())

# -------------------
# STEP 1: LOAD DATA
# -------------------

data <- readRDS("TCGA_panCancer_FOLH1.rds")
meta <- as.data.frame(colData(data))

# Keep only PAAD
paad_idx <- as.character(meta$project) == "TCGA-PAAD"
data_paad <- data[, paad_idx]
meta_paad <- as.data.frame(colData(data_paad))

# Store original sample IDs
original_sample_ids <- colnames(data_paad)

# -----------------------------
# IMPORT DFS CLINICAL TXT
# -----------------------------

cat("\n=== READING CLINICAL FILE ===\n")
clinical_df <- read.delim(
  "PAAD_clinical.txt.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Show all column names
cat("Clinical file columns:\n")
print(colnames(clinical_df))

# Find patient identifier column
patient_cols <- grep("Patient|patient|PatientID|patient_id|Sample|sample", colnames(clinical_df), value = TRUE)
cat("\nPatient-related columns found:\n")
print(patient_cols)

patient_col <- patient_cols[1]
cat("Using patient column:", patient_col, "\n")

# Extract DFS data
dfs_time_col <- grep("Disease Free \\(Months\\)|DFS \\(Months\\)", colnames(clinical_df), value = TRUE)[1]
dfs_status_col <- grep("Disease Free Status|DFS Status", colnames(clinical_df), value = TRUE)[1]

cat("DFS time column:", dfs_time_col, "\n")
cat("DFS status column:", dfs_status_col, "\n")

# Extract the three columns we need
dfs_extract <- clinical_df[, c(patient_col, dfs_time_col, dfs_status_col), drop = FALSE]

# Create clean column names
colnames(dfs_extract) <- c("patient", "DFS_months", "DFS_status")

# Convert to character
dfs_extract$patient <- as.character(dfs_extract$patient)
dfs_extract$DFS_months <- as.character(dfs_extract$DFS_months)
dfs_extract$DFS_status <- as.character(dfs_extract$DFS_status)

# Remove duplicates - keep first occurrence per patient
dfs_extract <- dfs_extract[!duplicated(dfs_extract$patient), ]

# Show what we extracted
cat("\nExtracted DFS data (first 5 rows):\n")
print(head(dfs_extract, 5))
cat("\nDFS_status values:\n")
print(table(dfs_extract$DFS_status, useNA = "ifany"))
cat("Number of unique patients in DFS data:", nrow(dfs_extract), "\n")

# Now merge with meta_paad
# First, ensure meta_paad has a patient column
meta_paad$patient <- as.character(meta_paad$patient)

cat("\nmeta_paad patient column (first 5):\n")
print(head(meta_paad$patient, 5))
cat("Number of samples in meta_paad:", nrow(meta_paad), "\n")

# Merge - this ADDS DFS columns to meta_paad
meta_paad <- merge(
  meta_paad,
  dfs_extract,
  by = "patient",
  all.x = TRUE,
  sort = FALSE
)

# Verify merge worked
cat("\n=== AFTER MERGE ===\n")
cat("DFS_months column exists:", "DFS_months" %in% colnames(meta_paad), "\n")
cat("DFS_status column exists:", "DFS_status" %in% colnames(meta_paad), "\n")
cat("DFS_months missing:", sum(is.na(meta_paad$DFS_months)), "\n")
cat("DFS_status missing:", sum(is.na(meta_paad$DFS_status)), "\n")
cat("DFS_status values in merged data:\n")
print(table(meta_paad$DFS_status, useNA = "ifany"))

# DON'T set rownames to patient - keep as sample IDs
# Instead, just ensure the data is properly aligned
rownames(meta_paad) <- original_sample_ids

# Reorder to match expression matrix
meta_paad <- meta_paad[original_sample_ids, , drop = FALSE]

count_mat <- assay(data_paad)
gene_annot <- as.data.frame(rowData(data_paad))

# ---------------------------
# STEP 2: NORMALIZATION + VST
# ---------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta_paad,
  design = ~ 1
)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
norm_mat <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

# Reorder metadata to match expression matrix columns
meta_paad <- meta_paad[colnames(expr_mat), , drop = FALSE]
gene_annot <- gene_annot[rownames(dds), , drop = FALSE]

gene_symbols <- if ("gene_name" %in% colnames(gene_annot)) {
  gene_annot$gene_name
} else if ("symbol" %in% colnames(gene_annot)) {
  gene_annot$symbol
} else {
  rownames(gene_annot)
}

rownames(expr_mat) <- rownames(dds)


# ----------------------------------------
# TARGET + EXPANDED BIOMARKER DEFINITIONS
# ----------------------------------------

target_gene <- "FOLH1"

biomarkers <- c(
  "KLK3", "AR", "NKX3-1", "TMPRSS2", "AMACR", "ERG", "MKI67", "CDH1", "VIM",
  "AURKA", "BRCA1", "BRCA2", "MYC", "TP53", "PTEN", "RB1",
  "TTF1", "INSM1", "NKX2-1", "ACP3", "CHGA", "CHGB", "TFRC", "NCAM1", "SCG2", "SYP",
  "FOLH1"
)

biomarkers <- unique(biomarkers)

# ----------------------------------------
# STEP 2B: CHECK BIOMARKER ANNOTATION
# ----------------------------------------

biomarker_presence <- data.frame(
  gene = biomarkers,
  found_in_annotation = biomarkers %in% gene_symbols,
  stringsAsFactors = FALSE
)

print(biomarker_presence)

cat(
  "Found",
  sum(biomarker_presence$found_in_annotation),
  "of",
  nrow(biomarker_presence),
  "biomarkers in gene annotation.\n"
)

if (!all(biomarker_presence$found_in_annotation)) {
  cat(
    "Missing biomarkers:",
    paste(biomarker_presence$gene[!biomarker_presence$found_in_annotation], collapse = ", "),
    "\n"
  )
}

write.csv(
  biomarker_presence,
  file = file.path(outdir, "biomarker_presence_check.csv"),
  row.names = FALSE
)

gene_to_row <- function(gene) {
  idx <- which(gene_symbols == gene)
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

target_idx <- gene_to_row(target_gene)
stopifnot(!is.na(target_idx))

bio_idx <- sapply(biomarkers, gene_to_row)
avail_mask <- !is.na(bio_idx)
biomarkers_avail <- biomarkers[avail_mask]
bio_idx <- bio_idx[avail_mask]

missing_genes <- biomarkers[!avail_mask]
if (length(missing_genes) > 0) {
  warning(
    "These biomarkers were not found and were skipped: ",
    paste(missing_genes, collapse = ", ")
  )
}

stopifnot(length(biomarkers_avail) > 0)

# PRIMARY TUMOR SAMPLES ONLY
cat("\n=== IDENTIFYING PRIMARY TUMOR SAMPLES ===\n")

# Use TCGA barcode sample type codes (01 = Primary) - most reliable
sample_type_codes <- substr(colnames(data_paad), 14, 15)
tumor_samples <- colnames(data_paad)[sample_type_codes == "01"]
tumor_samples <- intersect(tumor_samples, colnames(expr_mat))
cat("Number of primary tumor samples:", length(tumor_samples), "\n")
stopifnot(length(tumor_samples) > 2)

# ----------------------------
# STEP 3: SPEARMAN CORRELATION
# ----------------------------

cor_results <- function(mat, target_idx, gene_idx_vec, sample_ids, method = "spearman") {
  sample_ids <- intersect(sample_ids, colnames(mat))
  
  res <- lapply(seq_along(gene_idx_vec), function(i) {
    gidx <- gene_idx_vec[i]
    gene_name <- gene_symbols[gidx]
    x <- as.numeric(mat[gidx, sample_ids])
    y <- as.numeric(mat[target_idx, sample_ids])
    ok <- is.finite(x) & is.finite(y)
    
    if (sum(ok) < 3) {
      return(data.frame(
        gene = gene_name,
        method = method,
        cor = NA_real_,
        p.value = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
    
    data.frame(
      gene = gene_name,
      method = method,
      cor = unname(ct$estimate),
      p.value = ct$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  out <- bind_rows(res)
  out$p.adj <- p.adjust(out$p.value, method = "BH")
  out
}

tumor_cor_spear <- cor_results(
  mat = expr_mat,
  target_idx = target_idx,
  gene_idx_vec = bio_idx,
  sample_ids = tumor_samples,
  method = "spearman"
)

tumor_cor_spear <- tumor_cor_spear %>%
  arrange(desc(cor))

print(tumor_cor_spear)

write.csv(
  tumor_cor_spear,
  file = file.path(outdir, "primary_tumor_spearman_FOLH1_vs_biomarkers.csv"),
  row.names = FALSE
)

# ---------------------------
# STEP 4: BAR CHARTS
# ---------------------------

plot_df <- tumor_cor_spear %>%
  filter(gene != "FOLH1") %>%  # <-- ADD THIS LINE TO REMOVE FOLH1
  mutate(
    gene = factor(gene, levels = rev(gene)),
    direction = ifelse(cor >= 0, "Positive", "Negative")
  )

p_bar <- ggplot(plot_df, aes(x = cor, y = gene, fill = direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.7, colour = "black") +
  scale_fill_manual(values = c("Positive" = "#E64B35FF", "Negative" = "#4DBBD5FF")) +
  labs(
    title = "Correlation of FOLH1 with known biomarkers (TCGA-PAAD)",
    x = "Spearman correlation with FOLH1",
    y = NULL
  ) +
  theme_classic(base_size = 16, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title = element_text(face = "bold", colour = "black", size = 15),
    axis.text = element_text(colour = "black", size = 13),
    axis.line = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.7, colour = "black"),
    legend.position = "none"
  )

ggsave(
  filename = file.path(outdir, "primary_tumor_spearman_barplot_FOLH1_biomarkers.png"),
  plot = p_bar,
  width = 9,
  height = 7,
  dpi = 300
)

# -------------------------
# STEP 5: SUMMARY TABLE
# -------------------------

summary_df <- lapply(c(target_gene, biomarkers_avail), function(g) {
  ridx <- gene_to_row(g)
  x <- as.numeric(expr_mat[ridx, tumor_samples])
  
  data.frame(
    gene = g,
    n = sum(is.finite(x)),
    mean_expr = mean(x, na.rm = TRUE),
    median_expr = median(x, na.rm = TRUE),
    sd_expr = sd(x, na.rm = TRUE),
    min_expr = min(x, na.rm = TRUE),
    max_expr = max(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(
  summary_df,
  file = file.path(outdir, "primary_tumor_expression_summary.csv"),
  row.names = FALSE
)

# ---------------------------------------------
# STEP 6: SCATTERPLOTS
# ---------------------------------------------

plot_scatter_single <- function(gene, ridx) {
  df <- data.frame(
    FOLH1 = as.numeric(expr_mat[target_idx, tumor_samples]),
    biomarker = as.numeric(expr_mat[ridx, tumor_samples]),
    sample = tumor_samples,
    stringsAsFactors = FALSE
  )
  
  ct <- suppressWarnings(cor.test(df$FOLH1, df$biomarker, method = "spearman", exact = FALSE))
  lab <- paste0(
    "Spearman rho = ", round(unname(ct$estimate), 3),
    ", p = ", signif(ct$p.value, 3)
  )
  
  p <- ggplot(df, aes(x = FOLH1, y = biomarker)) +
    geom_point(size = 1.4, alpha = 0.7, color = "#E64B35FF") +
    geom_smooth(method = "lm", se = FALSE, color = "#3C5488FF", linewidth = 0.8) +  # <-- Changed se = FALSE
    labs(
      x = "FOLH1 (VST)",
      y = paste0(gene, " (VST)"),
      title = paste0("FOLH1 vs ", gene, " in primary TCGA-PAAD tumors"),
      subtitle = lab
    ) +
    theme_nature()
  
  ggsave(
    filename = file.path(outdir, paste0("FOLH1_vs_", gene, "_primary_tumor_scatter.png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

for (i in seq_along(biomarkers_avail)) {
  plot_scatter_single(biomarkers_avail[i], bio_idx[i])
}

# -----------------------------
# STEP 7: COMBINED SCATTER GRID
# -----------------------------

plot_scatter_grid_single <- function(gene, ridx) {
  df <- data.frame(
    FOLH1 = as.numeric(expr_mat[target_idx, tumor_samples]),
    biomarker = as.numeric(expr_mat[ridx, tumor_samples]),
    sample = tumor_samples,
    stringsAsFactors = FALSE
  )
  
  ct <- suppressWarnings(cor.test(df$FOLH1, df$biomarker, method = "spearman", exact = FALSE))
  lab <- paste0(
    "Spearman rho = ", round(unname(ct$estimate), 3),
    ", p = ", signif(ct$p.value, 3)
  )
  
  ggplot(df, aes(x = FOLH1, y = biomarker)) +
    geom_point(size = 1.4, alpha = 0.7, color = "#E64B35FF") +
    geom_smooth(method = "lm", se = TRUE, color = "#3C5488FF", linewidth = 0.8) +
    labs(
      x = "FOLH1 (VST)",
      y = paste0(gene, " (VST)"),
      title = gene,
      subtitle = lab
    ) +
    theme_nature() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      plot.subtitle = element_text(size = 8)
    )
}

scatter_plots <- lapply(seq_along(biomarkers_avail), function(i) {
  plot_scatter_grid_single(biomarkers_avail[i], bio_idx[i])
})

scatter_grid <- ggpubr::ggarrange(
  plotlist = scatter_plots,
  ncol = 2,
  nrow = ceiling(length(scatter_plots) / 2),
  align = "hv"
)

ggsave(
  filename = file.path(outdir, "primary_tumor_FOLH1_biomarker_scatter_grid.png"),
  plot = scatter_grid,
  width = 14,
  height = 16,
  dpi = 300
)

# -----------------------------------
# STEP 8: HEATMAP
# -----------------------------------

heat_idx <- c(target_idx, bio_idx)
heat_mat <- expr_mat[heat_idx, tumor_samples, drop = FALSE]
rownames(heat_mat) <- c(target_gene, biomarkers_avail)

annotation_col <- data.frame(
  Tissue = factor(rep("Primary Tumor", length(tumor_samples)))
)
rownames(annotation_col) <- tumor_samples

pheatmap::pheatmap(
  heat_mat,
  scale = "row",
  annotation_col = annotation_col,
  show_colnames = FALSE,
  main = "FOLH1 and known biomarkers in primary TCGA-PAAD tumors",
  filename = file.path(outdir, "primary_tumor_FOLH1_biomarker_heatmap.png"),
  width = 10,
  height = 7
)

# ----------------------------------------------------
# STEP 9: DFS KAPLAN-MEIER PLOTS
# ----------------------------------------------------

cat("\n=== DFS KM ANALYSIS ===\n")
cat("DFS_status column exists:", "DFS_status" %in% colnames(meta_paad), "\n")
if ("DFS_status" %in% colnames(meta_paad)) {
  print(table(meta_paad$DFS_status, useNA = "ifany"))
}

if ("DFS_status" %in% colnames(meta_paad) && sum(!is.na(meta_paad$DFS_status)) > 0) {
  
  meta_paad$DFS_status_clean <- trimws(as.character(meta_paad$DFS_status))
  
  meta_paad$DFS_event <- dplyr::case_when(
    meta_paad$DFS_status_clean == "1:Recurred/Progressed" ~ 1,
    meta_paad$DFS_status_clean == "0:DiseaseFree" ~ 0,
    TRUE ~ NA_real_
  )
  
  meta_paad$DFS_months <- suppressWarnings(as.numeric(meta_paad$DFS_months))
  
  cat("DFS_event missing:", sum(is.na(meta_paad$DFS_event)), "\n")
  cat("DFS_months missing:", sum(is.na(meta_paad$DFS_months)), "\n")
  
  plot_km_single <- function(gene, ridx) {
    expr_vec <- as.numeric(expr_mat[ridx, tumor_samples])
    
    df_km <- data.frame(
      sample = tumor_samples,
      patient = meta_paad[tumor_samples, "patient"],
      expr = expr_vec,
      DFS_months = suppressWarnings(as.numeric(meta_paad[tumor_samples, "DFS_months"])),
      DFS_event = suppressWarnings(as.numeric(meta_paad[tumor_samples, "DFS_event"])),
      stringsAsFactors = FALSE
    ) %>%
      filter(is.finite(expr), is.finite(DFS_months), is.finite(DFS_event))
    
    if (nrow(df_km) < 10) return(NULL)
    
    med_cut <- median(df_km$expr, na.rm = TRUE)
    df_km$group <- factor(ifelse(df_km$expr >= med_cut, "High", "Low"), levels = c("Low", "High"))
    
    if (length(unique(df_km$group)) < 2) return(NULL)
    
    fit <- do.call(
      survfit,
      list(Surv(DFS_months, DFS_event) ~ group, data = df_km)
    )
    
    sdiff <- survdiff(Surv(DFS_months, DFS_event) ~ group, data = df_km)
    pval <- 1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1)
    
    cox_fit <- coxph(Surv(DFS_months, DFS_event) ~ group, data = df_km)
    cox_sum <- summary(cox_fit)
    
    hr <- unname(cox_sum$coefficients[1, "exp(coef)"])
    hr_low <- unname(cox_sum$conf.int[1, "lower .95"])
    hr_high <- unname(cox_sum$conf.int[1, "upper .95"])
    
    km_plot <- ggsurvplot(
      fit,
      data = df_km,
      pval = TRUE,
      risk.table = TRUE,
      conf.int = FALSE,
      palette = c("Low" = "#1F4E79", "High" = "#8B0000"),
      legend.title = NULL,
      legend.labs = c("Low", "High"),
      xlab = "Disease-free survival (months)",
      ylab = "Disease-free survival probability",
      title = paste0(gene, " DFS in primary TCGA-PAAD"),
      risk.table.height = 0.22,
      ggtheme = theme_nature()
    )
    
    km_plot$plot <- km_plot$plot +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        plot.subtitle = element_text(size = 10),
        legend.position = c(0.85, 0.85)
      )
    
    ggsave(
      filename = file.path(km_dir, paste0(gene, "_DFS_KM.png")),
      plot = km_plot$plot,
      width = 8,
      height = 8,
      dpi = 300
    )
    
    data.frame(
      gene = gene,
      n = nrow(df_km),
      cutoff_median = med_cut,
      n_low = sum(df_km$group == "Low"),
      n_high = sum(df_km$group == "High"),
      logrank_p = pval,
      HR_high_vs_low = hr,
      CI95_lower = hr_low,
      CI95_upper = hr_high,
      stringsAsFactors = FALSE
    )
  }
  
  km_results_list <- lapply(seq_along(biomarkers_avail), function(i) {
    plot_km_single(biomarkers_avail[i], bio_idx[i])
  })
  
  km_results_list <- Filter(Negate(is.null), km_results_list)
  km_results_df <- bind_rows(km_results_list)
  
  if (nrow(km_results_df) > 0) {
    km_results_df <- km_results_df %>% arrange(logrank_p)
  }
  
  if (nrow(km_results_df) > 0) {
    km_results_df <- km_results_df %>%
      arrange(logrank_p) %>%
      mutate(
        p_value = logrank_p,
        p_value_fmt = signif(logrank_p, 3)
      )
  }
  
  write.csv(
    km_results_df,
    file = file.path(outdir, "DFS_KM_results_all_biomarkers.csv"),
    row.names = FALSE
  )
  
} else {
  cat("\nWARNING: DFS_status column missing or empty. Skipping KM analysis.\n")
  km_results_df <- data.frame()
}

# ----------------------------------------------
# STEP 10: STAGE-BASED EXPRESSION ANALYSIS
# ----------------------------------------------

cat("\n=== STAGE-BASED EXPRESSION ANALYSIS ===\n")

# Check available stage columns
cat("Available columns with 'stage' or 'tumor':\n")
print(colnames(meta_paad)[grepl("stage|Stage|STAGE|tumor|Tumor|metast|Mets|METS|sample|Sample", colnames(meta_paad), ignore.case = TRUE)])

# Use sample_type_code from TCGA barcodes (more reliable)
meta_paad$sample_type_code <- substr(colnames(data_paad), 14, 15)
meta_paad$sample_type <- dplyr::case_when(
  meta_paad$sample_type_code == "01" ~ "Primary",
  meta_paad$sample_type_code == "06" ~ "Metastatic",
  meta_paad$sample_type_code == "05" ~ "Additional Primary",
  meta_paad$sample_type_code == "07" ~ "Additional Metastatic",
  TRUE ~ "Other"
)

cat("\nSample type distribution:\n")
print(table(meta_paad$sample_type, useNA = "ifany"))

# Filter samples by type
primary_samples_type <- colnames(data_paad)[meta_paad$sample_type == "Primary" & colnames(data_paad) %in% colnames(expr_mat)]
metastatic_samples_type <- colnames(data_paad)[meta_paad$sample_type %in% c("Metastatic", "Additional Metastatic") & colnames(data_paad) %in% colnames(expr_mat)]

cat("Primary samples:", length(primary_samples_type), "\n")
cat("Metastatic samples:", length(metastatic_samples_type), "\n")

if (length(primary_samples_type) >= 3 && length(metastatic_samples_type) >= 3) {
  
  stage_comparison_df <- data.frame(
    gene = character(),
    group = character(),
    expression = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(biomarkers_avail)) {
    gene <- biomarkers_avail[i]
    ridx <- bio_idx[i]
    
    expr_primary <- as.numeric(expr_mat[ridx, primary_samples_type])
    expr_mets <- as.numeric(expr_mat[ridx, metastatic_samples_type])
    
    stage_comparison_df <- rbind(
      stage_comparison_df,
      data.frame(
        gene = gene,
        group = c(rep("Primary", length(expr_primary)), rep("Metastatic", length(expr_mets))),
        expression = c(expr_primary, expr_mets),
        stringsAsFactors = FALSE
      )
    )
  }
  
  stage_comparison_df$gene <- factor(stage_comparison_df$gene, levels = biomarkers_avail)
  stage_comparison_df$group <- factor(stage_comparison_df$group, levels = c("Primary", "Metastatic"))
  
  p_stage_all <- ggplot(stage_comparison_df, aes(x = gene, y = expression, fill = group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "black") +
    scale_fill_manual(values = c("Primary" = "#4DBBD5FF", "Metastatic" = "#E64B35FF")) +
    labs(
      title = "FOLH1 and biomarkers: Primary vs Metastatic (TCGA-PAAD)",
      x = "Gene",
      y = "Expression (VST)"
    ) +
    theme_nature() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
      plot.title = element_text(size = 14)
    )
  
  ggsave(
    filename = file.path(stage_dir, "all_biomarkers_primary_vs_metastatic_boxplot.png"),
    plot = p_stage_all,
    width = 12,
    height = 7,
    dpi = 300
  )
  
  genes_to_plot <- c(target_gene, head(biomarkers_avail[!biomarkers_avail %in% target_gene], 10))
  
  for (gene in genes_to_plot) {
    ridx <- gene_to_row(gene)
    
    df_gene <- data.frame(
      gene = gene,
      group = c(rep("Primary", length(as.numeric(expr_mat[ridx, primary_samples_type]))), 
                rep("Metastatic", length(as.numeric(expr_mat[ridx, metastatic_samples_type])))),
      expression = c(as.numeric(expr_mat[ridx, primary_samples_type]), 
                     as.numeric(expr_mat[ridx, metastatic_samples_type])),
      stringsAsFactors = FALSE
    )
    
    if (nrow(df_gene) < 6) next
    
    expr_p <- df_gene$expression[df_gene$group == "Primary"]
    expr_m <- df_gene$expression[df_gene$group == "Metastatic"]
    
    if (length(expr_p) >= 3 && length(expr_m) >= 3) {
      wt <- wilcox.test(expr_p, expr_m)
      p_label <- paste0("Wilcoxon p = ", signif(wt$p.value, 3))
      
      y_max <- max(df_gene$expression, na.rm = TRUE)
      y_pos <- y_max + 0.1 * diff(range(df_gene$expression, na.rm = TRUE))
      
      p_gene <- ggplot(df_gene, aes(x = group, y = expression, fill = group)) +
        geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
        geom_jitter(width = 0.15, size = 1.4, alpha = 0.7, color = "black") +
        scale_fill_manual(values = c("Primary" = "#4DBBD5FF", "Metastatic" = "#E64B35FF")) +
        annotate("text", x = 1.5, y = y_pos, label = p_label, size = 5, fontface = "bold") +
        labs(
          title = paste0(gene, " expression by tumor type"),
          x = "Tumor type",
          y = "Expression (VST)"
        ) +
        theme_nature()
      
      ggsave(
        filename = file.path(stage_dir, paste0(gene, "_primary_vs_metastatic_boxplot.png")),
        plot = p_gene,
        width = 7,
        height = 6,
        dpi = 300
      )
    }
  }
  
  stage_stats_df <- lapply(biomarkers_avail, function(gene) {
    ridx <- gene_to_row(gene)
    
    expr_p <- as.numeric(expr_mat[ridx, primary_samples_type])
    expr_m <- as.numeric(expr_mat[ridx, metastatic_samples_type])
    
    if (length(expr_p) < 3 || length(expr_m) < 3) {
      return(data.frame(
        gene = gene,
        n_primary = length(expr_p),
        n_metastatic = length(expr_m),
        mean_primary = NA_real_,
        mean_metastatic = NA_real_,
        median_primary = NA_real_,
        median_metastatic = NA_real_,
        p_wilcox = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    wt <- wilcox.test(expr_p, expr_m)
    
    data.frame(
      gene = gene,
      n_primary = length(expr_p),
      n_metastatic = length(expr_m),
      mean_primary = mean(expr_p, na.rm = TRUE),
      mean_metastatic = mean(expr_m, na.rm = TRUE),
      median_primary = median(expr_p, na.rm = TRUE),
      median_metastatic = median(expr_m, na.rm = TRUE),
      p_wilcox = wt$p.value,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  stage_stats_df <- stage_stats_df %>%
    mutate(
      p_adj = p.adjust(p_wilcox, method = "BH"),
      log2FC = log2((mean_metastatic + 1) / (mean_primary + 1))
    ) %>%
    arrange(p_adj)
  
  write.csv(
    stage_stats_df,
    file = file.path(stage_dir, "primary_vs_metastatic_expression_statistics.csv"),
    row.names = FALSE
  )
} else {
  cat("\nWARNING: Insufficient metastatic samples for comparison.\n")
  cat("Primary:", length(primary_samples_type), ", Metastatic:", length(metastatic_samples_type), "\n")
}

# ---------------------------
# STEP 10B: PATHOLOGICAL STAGE
# ---------------------------

if ("pathological_stage" %in% colnames(meta_paad)) {
  
  stage_col <- "pathological_stage"
  
  stage_df <- meta_paad %>%
    filter(sample_type == "Primary") %>%
    filter(!is.na(.[[stage_col]])) %>%
    mutate(
      stage_simple = case_when(
        grepl("Stage I", .[[stage_col]]) ~ "Stage I",
        grepl("Stage II", .[[stage_col]]) ~ "Stage II",
        grepl("Stage III", .[[stage_col]]) ~ "Stage III",
        grepl("Stage IV", .[[stage_col]]) ~ "Stage IV",
        TRUE ~ as.character(.[[stage_col]])
      )
    )
  
  cat("\nStage distribution (simplified):\n")
  print(table(stage_df$stage_simple, useNA = "ifany"))
  
  if (length(unique(stage_df$stage_simple)) >= 2) {
    
    stage_expr_list <- lapply(biomarkers_avail, function(gene) {
      ridx <- gene_to_row(gene)
      
      df_stage <- data.frame(
        gene = gene,
        stage = stage_df$stage_simple,
        expression = as.numeric(expr_mat[ridx, stage_df$patient]),
        stringsAsFactors = FALSE
      )
      
      df_stage <- df_stage %>% filter(is.finite(expression))
      
      if (nrow(df_stage) < 5) return(NULL)
      
      df_stage
    })
    
    stage_expr_list <- Filter(Negate(is.null), stage_expr_list)
    stage_expr_df <- bind_rows(stage_expr_list)
    
    if (nrow(stage_expr_df) > 0) {
      
      stage_expr_df$gene <- factor(stage_expr_df$gene, levels = biomarkers_avail)
      
      p_stage_overall <- ggplot(stage_expr_df, aes(x = stage, y = expression, fill = stage)) +
        geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
        geom_jitter(width = 0.12, size = 0.6, alpha = 0.4, color = "black") +
        scale_fill_brewer(palette = "Set2") +
        facet_wrap(~gene, scales = "free_y", ncol = 4) +
        labs(
          title = "Biomarker expression by pathological stage (TCGA-PAAD)",
          x = "Pathological Stage",
          y = "Expression (VST)"
        ) +
        theme_nature() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text = element_text(size = 8),
          plot.title = element_text(size = 12)
        )
      
      ggsave(
        filename = file.path(stage_dir, "all_biomarkers_by_pathological_stage_facet.png"),
        plot = p_stage_overall,
        width = 14,
        height = 12,
        dpi = 300
      )
      
      for (gene in c(target_gene, biomarkers_avail)) {
        df_gene <- stage_expr_df %>% filter(gene == !!gene)
        
        if (nrow(df_gene) < 5 || length(unique(df_gene$stage)) < 2) next
        
        kw_test <- kruskal.test(expression ~ stage, data = df_gene)
        p_label <- paste0("Kruskal-Wallis p = ", signif(kw_test$p.value, 3))
        
        p_gene_stage <- ggplot(df_gene, aes(x = stage, y = expression, fill = stage)) +
          geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
          geom_jitter(width = 0.15, size = 1.2, alpha = 0.6, color = "black") +
          scale_fill_brewer(palette = "Set2") +
          annotate("text", x = 1.5, y = max(df_gene$expression, na.rm = TRUE) * 1.05, 
                   label = p_label, size = 4, fontface = "bold") +
          labs(
            title = paste0(gene, " by pathological stage"),
            x = "Pathological Stage",
            y = "Expression (VST)"
          ) +
          theme_nature() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
          )
        
        ggsave(
          filename = file.path(stage_dir, paste0(gene, "_by_pathological_stage_boxplot.png")),
          plot = p_gene_stage,
          width = 7,
          height = 6,
          dpi = 300
        )
      }
      
      stage_stats_path_df <- lapply(biomarkers_avail, function(gene) {
        df_gene <- stage_expr_df %>% filter(gene == !!gene)
        
        if (nrow(df_gene) < 5 || length(unique(df_gene$stage)) < 2) {
          return(data.frame(
            gene = gene,
            n = nrow(df_gene),
            n_stages = length(unique(df_gene$stage)),
            p_kruskal = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
        
        kw_test <- kruskal.test(expression ~ stage, data = df_gene)
        
        data.frame(
          gene = gene,
          n = nrow(df_gene),
          n_stages = length(unique(df_gene$stage)),
          p_kruskal = kw_test$p.value,
          stringsAsFactors = FALSE
        )
      }) %>% bind_rows()
      
      stage_stats_path_df <- stage_stats_path_df %>%
        mutate(
          p_adj = p.adjust(p_kruskal, method = "BH")
        ) %>%
        arrange(p_adj)
      
      write.csv(
        stage_stats_path_df,
        file = file.path(stage_dir, "expression_by_pathological_stage_statistics.csv"),
        row.names = FALSE
      )
    }
  }
} else {
  cat("\nWARNING: pathological_stage column not found in metadata.\n")
}

# --------------------------
# STEP 11: EXCEL OUTPUT
# --------------------------

wb <- createWorkbook()

addWorksheet(wb, "Biomarker_presence")
writeData(wb, "Biomarker_presence", biomarker_presence)

addWorksheet(wb, "Spearman_correlations")
writeData(wb, "Spearman_correlations", tumor_cor_spear)

addWorksheet(wb, "Expression_summary")
writeData(wb, "Expression_summary", summary_df)

addWorksheet(wb, "DFS_KM_results")
writeData(wb, "DFS_KM_results", km_results_df)

if (exists("stage_stats_df")) {
  addWorksheet(wb, "Primary_vs_Mets")
  writeData(wb, "Primary_vs_Mets", stage_stats_df)
}

if (exists("stage_stats_path_df")) {
  addWorksheet(wb, "By_Path_Stage")
  writeData(wb, "By_Path_Stage", stage_stats_path_df)
}

saveWorkbook(
  wb,
  file = file.path(outdir, "PAAD_FOLH1_biomarker_stage_results_DFS.xlsx"),
  overwrite = TRUE
)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Output directory:", outdir, "\n")
cat("KM plots directory:", km_dir, "\n")
cat("Stage plots directory:", stage_dir, "\n")