# ==================================================================
# Ravindu: Comparison of FOLH1 with Known Prostate Cancer Biomarkers
# Expanded biomarker panel + annotation presence check + DFS KM plots
# ==================================================================

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

outdir <- "output_new_biomarkers_included"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

km_dir <- file.path(outdir, "KM plots")
dir.create(km_dir, showWarnings = FALSE, recursive = TRUE)

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

# Keep only PRAD
prad_idx <- as.character(meta$project) == "TCGA-PRAD"
data_prad <- data[, prad_idx]
meta_prad <- as.data.frame(colData(data_prad))

# -----------------------------
# IMPORT DFS CLINICAL TXT
# -----------------------------

clinical_df <- read.delim(
  "PRAD_clinical.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dfs_df <- clinical_df[, c("Patient Identifier", "Disease Free (Months)", "Disease Free Status"), drop = FALSE]

dfs_df$`Patient Identifier` <- as.character(dfs_df$`Patient Identifier`)
meta_prad$patient <- as.character(meta_prad$patient)

meta_prad <- merge(
  meta_prad,
  dfs_df,
  by.x = "patient",
  by.y = "Patient Identifier",
  all.x = TRUE,
  sort = FALSE
)

names(meta_prad)[names(meta_prad) == "Disease Free (Months)"] <- "DFS_months"
names(meta_prad)[names(meta_prad) == "Disease Free Status"] <- "DFS_status"

print(c("DFS_months", "DFS_status") %in% colnames(meta_prad))
cat("DFS_months missing:", sum(is.na(meta_prad$DFS_months)), "\n")
cat("DFS_status missing:", sum(is.na(meta_prad$DFS_status)), "\n")

# Keep ONLY PRIMARY tumors
stopifnot("tumor_descriptor" %in% colnames(meta_prad))
keep_idx <- as.character(meta_prad$tumor_descriptor) %in% "Primary"
data_prad_keep <- data_prad[, keep_idx]
meta_prad_keep <- as.data.frame(colData(data_prad_keep))
rownames(meta_prad_keep) <- colnames(data_prad_keep)

# Bring DFS columns into primary-only metadata by patient
meta_prad_keep$patient <- as.character(meta_prad_keep$patient)
dfs_map <- meta_prad[, c("patient", "DFS_months", "DFS_status"), drop = FALSE] %>%
  distinct(patient, .keep_all = TRUE)

meta_prad_keep <- meta_prad_keep %>%
  left_join(dfs_map, by = "patient")

rownames(meta_prad_keep) <- colnames(data_prad_keep)

count_mat <- assay(data_prad_keep)
gene_annot <- as.data.frame(rowData(data_prad_keep))

# ---------------------------
# STEP 2: NORMALIZATION + VST
# ---------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta_prad_keep,
  design = ~ 1
)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)
norm_mat <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

meta_prad_keep <- meta_prad_keep[colnames(expr_mat), , drop = FALSE]
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
  "AURKA", "BRCA1", "BRCA2", "MYC", "TP53", "PTEN", "RB1", "TTF1", "INSM1",
  "NKX2-1", "ACP3", "CHGA", "CHGB", "TFRC", "NCAM1", "SCG2", "SYP", "FOLH1"
)

# ----------------------------------------
# STEP 2B: CHECK BIOMARKER ANNOTATION HIT
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

if (all(biomarker_presence$found_in_annotation)) {
  cat("All biomarkers were found in the gene annotation.\n")
} else {
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
tumor_samples <- rownames(meta_prad_keep)[as.character(meta_prad_keep$tumor_descriptor) == "Primary"]
tumor_samples <- intersect(tumor_samples, colnames(expr_mat))
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
# STEP 4: PLOTTING BAR CHARTS
# ---------------------------

plot_df <- tumor_cor_spear %>%
  mutate(
    gene = factor(gene, levels = rev(gene)),
    direction = ifelse(cor >= 0, "Positive", "Negative")
  )

p_bar <- ggplot(plot_df, aes(x = cor, y = gene, fill = direction)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linetype = 1, linewidth = 0.7, colour = "black") +
  scale_fill_manual(values = c("Positive" = "#E64B35FF", "Negative" = "#4DBBD5FF")) +
  labs(
    title = "Correlation of FOLH1 with known biomarkers",
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
# STEP 6: SCATTERPLOTS: FOLH1 VS EACH BIOMARKER
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
    geom_smooth(method = "lm", se = TRUE, color = "#3C5488FF", linewidth = 0.8) +
    labs(
      x = "FOLH1 (VST)",
      y = paste0(gene, " (VST)"),
      title = paste0("FOLH1 vs ", gene, " in primary TCGA-PRAD tumors"),
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

# --------------------------
# STEP 8: ERG FUSION BOXPLOT
# --------------------------

erg_df <- data.frame(
  sample = rownames(meta_prad_keep),
  ERG_status = as.character(meta_prad_keep$paper_ERG_status),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(ERG_status)) %>%
  filter(ERG_status %in% c("fusion", "none")) %>%
  mutate(ERG_status = factor(ERG_status, levels = c("none", "fusion")))

if (nrow(erg_df) > 2) {
  erg_df$FOLH1 <- as.numeric(expr_mat[target_idx, erg_df$sample])
  erg_df <- erg_df %>% filter(is.finite(FOLH1))
  
  erg_test <- wilcox.test(FOLH1 ~ ERG_status, data = erg_df)
  
  p_label <- paste0(
    "Wilcoxon p = ",
    signif(erg_test$p.value, 3)
  )
  
  y_pos <- max(erg_df$FOLH1, na.rm = TRUE) + 0.4 * diff(range(erg_df$FOLH1, na.rm = TRUE))
  
  p_erg <- ggplot(erg_df, aes(x = ERG_status, y = FOLH1, fill = ERG_status)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.15, size = 1.4, alpha = 0.7, color = "black") +
    scale_fill_manual(values = c("none" = "#4DBBD5FF", "fusion" = "#E64B35FF")) +
    annotate("text", x = 1.5, y = y_pos, label = p_label, size = 5, fontface = "bold") +
    labs(
      title = "FOLH1 expression by ERG fusion status",
      x = "ERG status",
      y = "FOLH1 (VST)"
    ) +
    theme_nature()
  
  ggsave(
    filename = file.path(outdir, "FOLH1_by_ERG_status_boxplot.png"),
    plot = p_erg,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  writeLines(
    capture.output(erg_test),
    file.path(outdir, "FOLH1_by_ERG_status_wilcox.txt")
  )
}

# -----------------------------------
# STEP 9: HEATMAP: FOLH1 + BIOMARKERS
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
  main = "FOLH1 and known biomarkers in primary TCGA-PRAD tumors",
  filename = file.path(outdir, "primary_tumor_FOLH1_biomarker_heatmap.png"),
  width = 10,
  height = 7
)

# ----------------------------------------------------
# STEP 10: DFS KAPLAN-MEIER PLOTS FOR EACH BIOMARKER
# ----------------------------------------------------

print(table(meta_prad_keep$DFS_status, useNA = "ifany"))

meta_prad_keep$DFS_status_clean <- trimws(as.character(meta_prad_keep$DFS_status))

meta_prad_keep$DFS_event <- dplyr::case_when(
  meta_prad_keep$DFS_status_clean == "1:Recurred/Progressed" ~ 1,
  meta_prad_keep$DFS_status_clean == "0:DiseaseFree" ~ 0,
  TRUE ~ NA_real_
)

meta_prad_keep$DFS_months <- suppressWarnings(as.numeric(meta_prad_keep$DFS_months))

cat("DFS_event missing:", sum(is.na(meta_prad_keep$DFS_event)), "\n")
cat("DFS_months missing:", sum(is.na(meta_prad_keep$DFS_months)), "\n")

plot_km_single <- function(gene, ridx) {
  expr_vec <- as.numeric(expr_mat[ridx, tumor_samples])
  
  df_km <- data.frame(
    sample = tumor_samples,
    patient = meta_prad_keep[tumor_samples, "patient"],
    expr = expr_vec,
    DFS_months = suppressWarnings(as.numeric(meta_prad_keep[tumor_samples, "DFS_months"])),
    DFS_event = suppressWarnings(as.numeric(meta_prad_keep[tumor_samples, "DFS_event"])),
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
    title = paste0(gene, " DFS in primary TCGA-PRAD"),
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

if (exists("erg_df")) {
  addWorksheet(wb, "ERG_FOLH1")
  writeData(wb, "ERG_FOLH1", erg_df)
}

saveWorkbook(
  wb,
  file = file.path(outdir, "primary_tumor_FOLH1_biomarker_results.xlsx"),
  overwrite = TRUE
)