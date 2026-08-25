# ==============================================================
# Ravindu: Comparison of FOLH1 Expression Across TCGA Cancer Types
# ===============================================================

library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)

# =====================================
# STEP 1: Load TCGA pan-cancer object
# =====================================

tcga_obj <- readRDS("TCGA_panCancer_FOLH1.rds")

expr_mat <- assay(tcga_obj, "fpkm_uq_unstrand")
meta <- as.data.frame(colData(tcga_obj))

# Inspect expression matrix
dim(expr_mat)
head(rownames(expr_mat))
colnames(meta)

# Create gene annotation table from rowData()
gene_annot <- as.data.frame(rowData(tcga_obj))

print("Annotation table column names:")
print(colnames(gene_annot))

gene_name_col <- colnames(gene_annot)[
  colnames(gene_annot) %in% c("gene_name")
][1]

print(paste("Gene name column found:", gene_name_col))

gene_annot$ENSEMBL <- rownames(gene_annot)
ensembl_col <- "ENSEMBL"

# =====================
# STEP 2: Find FOLH1 
# =====================

target_gene <- "FOLH1"
match_idx <- match(target_gene, gene_annot[[gene_name_col]])

folh1_ensembl <- gene_annot[[ensembl_col]][match_idx]
print(paste("FOLH1 ENSEMBL ID:", folh1_ensembl))

clean_ensembls <- sub("\\..*", "", rownames(expr_mat))
clean_folh1 <- sub("\\..*", "", folh1_ensembl)

row_idx <- which(clean_ensembls == clean_folh1)

actual_rowname <- rownames(expr_mat)[row_idx[1]]

folh1_expr <- expr_mat[actual_rowname, ]

# =======================================
# STEP 3: Build expression data frame
# =======================================

sample_ids <- if ("sample_id" %in% colnames(meta)) {
  meta$sample_id
} else {
  colnames(expr_mat)
}

expression_df <- data.frame(
  sample_id = sample_ids,
  project_id = if ("project_id" %in% colnames(meta)) meta$project_id else NA,
  folh1_fpkm = as.numeric(folh1_expr),
  stringsAsFactors = FALSE
)

# ============================================
# STEP 4: Keep PRIMARY TUMOR samples only
# ============================================

expression_df$sample_type_code <- meta$tumor_descriptor

expression_df_primary <- expression_df %>%
  filter(sample_type_code == "Primary")

print("Sample counts before and after primary tumor filtering:")
print(paste("All samples:", nrow(expression_df)))
print(paste("Primary tumor only:", nrow(expression_df_primary)))

# ================================================
# STEP 5: Map TCGA project codes to cancer types
# ================================================

expression_df_primary$cancer_code <- sub("^TCGA-", "", expression_df_primary$project_id)

tcga_cancer_map <- c(
  "LGG"  = "Brain Lower Grade Glioma",
  "KICH" = "Kidney Chromophobe",
  "TGCT" = "Testicular Germ Cell Tumors",
  "UCS"  = "Uterine Carcinosarcoma",
  "PRAD" = "Prostate Adenocarcinoma",
  "BRCA" = "Breast Invasive Carcinoma",
  "LUAD" = "Lung Adenocarcinoma",
  "LUSC" = "Lung Squamous Cell Carcinoma",
  "KIRC" = "Kidney Renal Clear Cell Carcinoma",
  "KIRP" = "Kidney Renal Papillary Cell Carcinoma",
  "COAD" = "Colon Adenocarcinoma",
  "READ" = "Rectal Adenocarcinoma",
  "GBM"  = "Glioblastoma Multiforme",
  "HNSC" = "Head and Neck Squamous Cell Carcinoma",
  "THCA" = "Thyroid Carcinoma",
  "OV"   = "Ovarian Serous Cystadenocarcinoma",
  "UCEC" = "Uterine Corpus Endometrial Carcinoma",
  "CESC" = "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma",
  "LIHC" = "Liver Hepatocellular Carcinoma",
  "PAAD" = "Pancreatic Adenocarcinoma",
  "ESCA" = "Esophageal Carcinoma",
  "STAD" = "Stomach Adenocarcinoma",
  "BLCA" = "Bladder Urothelial Carcinoma",
  "SKCM" = "Skin Cutaneous Melanoma",
  "PCPG" = "Pheochromocytoma and Paraganglioma",
  "SARC" = "Sarcoma",
  "LAML" = "Acute Myeloid Leukemia",
  "ACC"  = "Adrenocortical Carcinoma",
  "MESO" = "Mesothelioma",
  "CHOL" = "Cholangiocarcinoma",
  "UVM"  = "Uveal Melanoma",
  "DLBC" = "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
  "THYM" = "Thymoma"
)

expression_df_primary$cancer_type <- tcga_cancer_map[expression_df_primary$cancer_code]

# =====================================
# STEP 6: Clean dataset for analysis
# =====================================

expression_df_clean <- expression_df_primary %>%
  filter(!is.na(cancer_type), cancer_type != "", !is.na(folh1_fpkm))

print("Samples retained for final analysis:")
print(nrow(expression_df_clean))

# Log transform for analysis and plotting
expression_df_clean <- expression_df_clean %>%
  mutate(folh1_log2_fpkm = log2(folh1_fpkm + 1))

# ====================================================
# STEP 7: Summarize FOLH1 expression by cancer type
# ====================================================

cancer_expression <- expression_df_clean %>%
  group_by(cancer_type) %>%
  summarise(
    median_fpkm = median(folh1_log2_fpkm, na.rm = TRUE),
    sd_fpkm = sd(folh1_log2_fpkm, na.rm = TRUE),
    median_log2_fpkm = median(folh1_log2_fpkm, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_fpkm))

print("=== TOP 10 CANCERS WITH HIGHEST FOLH1 EXPRESSION (PRIMARY TUMORS, FPKM_unstrand) ===")
print(
  cancer_expression %>%
    head(10) %>%
    dplyr::select(cancer_type, median_fpkm, n_samples)
)

# ===============================================================
# STEP 8: Visualization - Bar plot (TOP 15 by median log2 FPKM)
# ===============================================================

top_n <- 15

bar_plot <- cancer_expression %>%
  head(top_n) %>%
  ggplot(aes(x = reorder(cancer_type, median_fpkm), y = median_fpkm)) +
  geom_col(fill = "#4472C4", color = "black") +
  coord_flip() +
  labs(
    title = "FOLH1 (PSMA) Expression Across TCGA Pan-Cancer Types",
    subtitle = paste("Primary tumors only; Top", top_n, "cancers by median log2(FPKM_unstrand + 1)"),
    x = "Cancer Type",
    y = "Median log2(FPKM_unstrand + 1)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11),
    panel.grid = element_blank()
  )

ggsave("folh1_pan_cancer_barplot_fpkm_primary.png", bar_plot, width = 10, height = 8, dpi = 300)
print("Bar plot saved: folh1_pan_cancer_barplot_fpkm_primary.png")

# =======================================================
# STEP 9: Visualization - Box plot using log2(FPKM + 1)
# =======================================================

top_cancers <- cancer_expression %>%
  head(top_n) %>%
  pull(cancer_type)

box_plot <- expression_df_clean %>%
  filter(cancer_type %in% top_cancers) %>%
  ggplot(aes(x = cancer_type, y = folh1_log2_fpkm, fill = cancer_type)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  labs(
    title = "FOLH1 (PSMA) Expression Distribution Across Cancer Types",
    subtitle = paste("Primary tumors only; Top", top_n, "cancers by median log2(FPKM_unstrand + 1)"),
    x = "Cancer Type",
    y = "log2(FPKM_unstrand + 1)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none"
  )

ggsave("folh1_pan_cancer_boxplot_fpkm_primary.png", box_plot, width = 12, height = 7, dpi = 300)
print("Box plot saved: folh1_pan_cancer_boxplot_fpkm_primary.png")

# ======================================================
# STEP 10: Visualization - Dot plot with sample sizes
# ======================================================

top_df <- cancer_expression %>%
  head(top_n)

dot_plot <- top_df %>%
  ggplot(aes(
    x = reorder(cancer_type, median_fpkm),
    y = median_fpkm,
    size = n_samples,
    color = median_fpkm
  )) +
  geom_point(alpha = 0.8) +
  coord_flip() +
  scale_color_gradient(low = "#87CEEB", high = "#4472C4") +
  labs(
    title = "FOLH1 (PSMA) Median Expression by Cancer Type",
    subtitle = paste(
      "Primary tumors only; point size = sample size (n =",
      min(top_df$n_samples, na.rm = TRUE), "-",
      max(top_df$n_samples, na.rm = TRUE), ")"
    ),
    x = "Cancer Type",
    y = "Median log2(FPKM_unstrand + 1)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

ggsave("folh1_pan_cancer_dotplot_fpkm_primary.png", dot_plot, width = 10, height = 8, dpi = 300)
print("Dot plot saved: folh1_pan_cancer_dotplot_fpkm_primary.png")

# =========================
# STEP 13: Export results
# ==========================

write.csv(cancer_expression, "folh1_cancer_expression_summary_fpkm_primary.csv", row.names = FALSE)
write.csv(expression_df_clean, "folh1_full_expression_data_fpkm_primary.csv", row.names = FALSE)

