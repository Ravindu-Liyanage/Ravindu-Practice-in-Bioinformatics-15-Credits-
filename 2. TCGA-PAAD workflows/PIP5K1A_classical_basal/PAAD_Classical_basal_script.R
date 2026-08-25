# ========================================================================
# Ravindu: PIP5K1A Survival Analysis based on Classical and Basal Sub-types
# =========================================================================

library(SummarizedExperiment)

outdir <- "output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load data
data <- readRDS("TCGA_panCancer_FOLH1.rds")
count_mat <- assay(data)
meta <- as.data.frame(colData(data))
gene_annot <- as.data.frame(rowData(data))

# Extract PAAD samples
meta_paad <- meta[meta$project == "TCGA-PAAD", ]
sub_col <- "paper_mRNA.Moffitt.clusters..All.150.Samples..1basal..2classical"

# Keep only PAAD samples with non-missing sub-type labels
meta_paad2 <- meta_paad[!is.na(meta_paad[[sub_col]]), ]

# Match count matrix columns to the filtered metadata rows
count_paad <- count_mat[, rownames(meta_paad2)]

# -------------------------------------
# 1) Normalizing and variance stabilizing 
# -------------------------------------

library(DESeq2)
# keep as integer matrix for DESeq2
count_paad <- round(count_paad)

# build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_paad,
  colData = meta_paad2,
  design = ~ 1
)

# pre-filter very low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# variance-stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# normalized expression matrix
vst_mat <- assay(vsd)

# ---------------------
# 2) Split sub-type labels
# ---------------------

basal_samples <- rownames(meta_paad2)[meta_paad2[[sub_col]] == 1]
classical_samples <- rownames(meta_paad2)[meta_paad2[[sub_col]] == 2]

# Separate matrices
count_basal <- vst_mat[, basal_samples, drop = FALSE]
count_classical <- vst_mat[, classical_samples, drop = FALSE]

dim(count_basal)
dim(count_classical)

table(meta_paad2[[sub_col]])

# --------------------------------
# 3) Extracting the PIP5K1A data
# --------------------------------

# make sure gene_annot rownames are Ensembl IDs with versions, matching count_mat
gene_map <- data.frame(
  ensembl_id = rownames(gene_annot),
  gene_name = gene_annot$gene_name,
  stringsAsFactors = FALSE
)

pip5k1a_rows <- gene_map$ensembl_id[gene_map$gene_name == "PIP5K1A"]

pip5k1a_rows

# Extracting expression
pip5k1a_expr_basal <- count_basal[pip5k1a_rows, , drop = FALSE]
pip5k1a_expr_classical <- count_classical[pip5k1a_rows, , drop = FALSE]
rownames(pip5k1a_expr_basal) <- "PIP5K1A"
rownames(pip5k1a_expr_classical) <- "PIP5K1A"

# ----------------------------------------------------
# Setting high and low gene expression based on Median
# ----------------------------------------------------

# convert to data frames for easier handling
pip5k1a_basal_df <- as.data.frame(t(pip5k1a_expr_basal))
pip5k1a_basal_df$sample <- rownames(pip5k1a_basal_df)

pip5k1a_classical_df <- as.data.frame(t(pip5k1a_expr_classical))
pip5k1a_classical_df$sample <- rownames(pip5k1a_classical_df)

# add subtype labels
pip5k1a_basal_df$subtype <- "Basal-like"
pip5k1a_classical_df$subtype <- "Classical-like"

# median cutoffs within each subtype
basal_median <- median(pip5k1a_basal_df$PIP5K1A, na.rm = TRUE)
classical_median <- median(pip5k1a_classical_df$PIP5K1A, na.rm = TRUE)

# add high/low groups
pip5k1a_basal_df$expr_group <- ifelse(pip5k1a_basal_df$PIP5K1A >= basal_median, "High", "Low")
pip5k1a_classical_df$expr_group <- ifelse(pip5k1a_classical_df$PIP5K1A >= classical_median, "High", "Low")

# convert to factor
pip5k1a_basal_df$expr_group <- factor(pip5k1a_basal_df$expr_group, levels = c("Low", "High"))
pip5k1a_classical_df$expr_group <- factor(pip5k1a_classical_df$expr_group, levels = c("Low", "High"))

pip5k1a_all <- rbind(pip5k1a_basal_df, pip5k1a_classical_df)
write.csv(pip5k1a_all, file.path(outdir, "PIP5K1A_basal_classical_high_low.csv"), row.names = FALSE)

# --------------------------------------------
# Uploading new clinical data that include DFS
# ---------------------------------------------

# shorten expression matrix sample IDs to patient-level TCGA barcodes
pip5k1a_basal_df$sample <- substr(pip5k1a_basal_df$sample, 1, 12)
pip5k1a_classical_df$sample <- substr(pip5k1a_classical_df$sample, 1, 12)

clin <- read.delim("paad_tcga_gdc_clinical_data.tsv", check.names = FALSE, stringsAsFactors = FALSE)
clin_df <- clin[!is.na(clin$`Disease Free (Months)`) & !is.na(clin$`Disease Free Status`), ]

pip5k1a_all <- rbind(pip5k1a_basal_df, pip5k1a_classical_df)

merged_df <- merge(
  pip5k1a_all,
  clin_df,
  by.x = "sample",
  by.y = "Patient ID"
)

#============================
# DFS analysis and KM plots
#============================

library(survival)
library(survminer)

# clean column names if needed
merged_df$DFS_time <- as.numeric(merged_df$`Disease Free (Months)`)
merged_df$DFS_event <- ifelse(grepl("^1", merged_df$`Disease Free Status`), 1, 0)

# Basal-like subset
basal_df <- merged_df[merged_df$subtype == "Basal-like", ]
basal_df <- basal_df[!is.na(basal_df$DFS_time) & !is.na(basal_df$DFS_event) & !is.na(basal_df$expr_group), ]

fit_basal <- survfit(Surv(DFS_time, DFS_event) ~ expr_group, data = basal_df)

basal_plot <- ggsurvplot(
  fit_basal,
  data = basal_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  title = "PIP5K1A DFS in Basal-like subtype",
  legend.title = "PIP5K1A",
  legend.labs = c("Low", "High"),
  xlab = "Disease-free months",
  ylab = "Disease-free survival probability"
)

png(file.path(outdir, "PIP5K1A_DFS_Basal_like_KM.png"), width = 2200, height = 1800, res = 300)
print(basal_plot)
dev.off()

# Classical-like subset
classical_df <- merged_df[merged_df$subtype == "Classical-like", ]
classical_df <- classical_df[!is.na(classical_df$DFS_time) & !is.na(classical_df$DFS_event) & !is.na(classical_df$expr_group), ]

fit_classical <- survfit(Surv(DFS_time, DFS_event) ~ expr_group, data = classical_df)


classical_plot <- ggsurvplot(
  fit_classical,
  data = classical_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  title = "PIP5K1A DFS in Classical-like subtype",
  legend.title = "PIP5K1A",
  legend.labs = c("Low", "High"),
  xlab = "Disease-free months",
  ylab = "Disease-free survival probability"
)

png(file.path(outdir, "PIP5K1A_DFS_Classical_like_KM.png"), width = 2200, height = 1800, res = 300)
print(classical_plot)
dev.off()

#===============
# Log Rank Test
#===============

# Basal-like log-rank test
survdiff(Surv(DFS_time, DFS_event) ~ expr_group, data = basal_df)

# Classical-like log-rank test
survdiff(Surv(DFS_time, DFS_event) ~ expr_group, data = classical_df)

# Basal-like
sd_basal <- survdiff(Surv(DFS_time, DFS_event) ~ expr_group, data = basal_df)
p_basal <- 1 - pchisq(sd_basal$chisq, df = length(sd_basal$n) - 1)

# Classical-like
sd_classical <- survdiff(Surv(DFS_time, DFS_event) ~ expr_group, data = classical_df)
p_classical <- 1 - pchisq(sd_classical$chisq, df = length(sd_classical$n) - 1)

logrank_results <- data.frame(
  subtype = c("Basal-like", "Classical-like"),
  chisq = c(sd_basal$chisq, sd_classical$chisq),
  p_value = c(p_basal, p_classical)
)

write.csv(logrank_results, file.path(outdir, "PIP5K1A_DFS_logrank_results.csv"), row.names = FALSE)

#================
# Cox regression
#================

# Cox regression - Basal-like
cox_basal <- coxph(Surv(DFS_time, DFS_event) ~ expr_group, data = basal_df)
sum_basal <- summary(cox_basal)

# Cox regression - Classical-like
cox_classical <- coxph(Surv(DFS_time, DFS_event) ~ expr_group, data = classical_df)
sum_classical <- summary(cox_classical)

# Extract HR table
cox_results <- data.frame(
  subtype = c("Basal-like", "Classical-like"),
  HR = c(
    sum_basal$conf.int[1, "exp(coef)"],
    sum_classical$conf.int[1, "exp(coef)"]
  ),
  lower_95_CI = c(
    sum_basal$conf.int[1, "lower .95"],
    sum_classical$conf.int[1, "lower .95"]
  ),
  upper_95_CI = c(
    sum_basal$conf.int[1, "upper .95"],
    sum_classical$conf.int[1, "upper .95"]
  ),
  p_value = c(
    sum_basal$coefficients[1, "Pr(>|z|)"],
    sum_classical$coefficients[1, "Pr(>|z|)"]
  )
)

write.csv(cox_results, file.path(outdir, "PIP5K1A_DFS_cox_results.csv"), row.names = FALSE)

#==================
# Interaction Test
#==================

library(survival)

merged_df$subtype <- factor(merged_df$subtype, levels = c("Basal-like", "Classical-like"))
merged_df$PIP5K1A_z <- scale(merged_df$PIP5K1A)

# interaction model
cox_int <- coxph(Surv(DFS_time, DFS_event) ~ PIP5K1A_z * subtype, data = merged_df)
sum_int <- summary(cox_int)

# extract interaction p-value
int_p <- sum_int$coefficients[grep("PIP5K1A_z:subtype", rownames(sum_int$coefficients)), "Pr(>|z|)"]

# save full summary and p-value
interaction_results <- data.frame(
  term = rownames(sum_int$coefficients),
  coef = sum_int$coefficients[, "coef"],
  HR = sum_int$coefficients[, "exp(coef)"],
  lower_95_CI = sum_int$conf.int[, "lower .95"],
  upper_95_CI = sum_int$conf.int[, "upper .95"],
  p_value = sum_int$coefficients[, "Pr(>|z|)"]
)

interaction_results$signif <- cut(
  interaction_results$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns"),
  right = FALSE
)

write.csv(interaction_results, file.path(outdir, "PIP5K1A_DFS_interaction_results.csv"), row.names = FALSE)

#===========================
# Uploading OS clinical data
#===========================

# Check OS fields in clinical metadata
names(meta_paad)

# Make sure sample IDs are patient-level barcodes
pip5k1a_basal_df$sample <- substr(pip5k1a_basal_df$sample, 1, 12)
pip5k1a_classical_df$sample <- substr(pip5k1a_classical_df$sample, 1, 12)

# Keep only rows with OS information
clin_os <- clin[!is.na(clin$`Overall Survival (Months)`) & !is.na(clin$`Overall Survival Status`), ]

pip5k1a_all <- rbind(pip5k1a_basal_df, pip5k1a_classical_df)

merged_os_df <- merge(
  pip5k1a_all,
  clin_os,
  by.x = "sample",
  by.y = "Patient ID"
)

#==========================
# OS analysis and KM plots
#==========================

library(survival)
library(survminer)

merged_os_df$OS_time <- as.numeric(merged_os_df$`Overall Survival (Months)`)
merged_os_df$OS_event <- ifelse(grepl("^1", merged_os_df$`Overall Survival Status`), 1, 0)

# Basal-like subset
basal_os <- merged_os_df[merged_os_df$subtype == "Basal-like", ]
basal_os <- basal_os[!is.na(basal_os$OS_time) & !is.na(basal_os$OS_event) & !is.na(basal_os$expr_group), ]

fit_basal_os <- survfit(Surv(OS_time, OS_event) ~ expr_group, data = basal_os)

basal_os_plot <- ggsurvplot(
  fit_basal_os,
  data = basal_os,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  title = "PIP5K1A OS in Basal-like subtype",
  legend.title = "PIP5K1A",
  legend.labs = c("Low", "High"),
  xlab = "Overall survival months",
  ylab = "Overall survival probability"
)

png(file.path(outdir, "PIP5K1A_OS_Basal_like_KM.png"), width = 2200, height = 1800, res = 300)
print(basal_os_plot)
dev.off()

# Classical-like subset
classical_os <- merged_os_df[merged_os_df$subtype == "Classical-like", ]
classical_os <- classical_os[!is.na(classical_os$OS_time) & !is.na(classical_os$OS_event) & !is.na(classical_os$expr_group), ]

fit_classical_os <- survfit(Surv(OS_time, OS_event) ~ expr_group, data = classical_os)

classical_os_plot <- ggsurvplot(
  fit_classical_os,
  data = classical_os,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  title = "PIP5K1A OS in Classical-like subtype",
  legend.title = "PIP5K1A",
  legend.labs = c("Low", "High"),
  xlab = "Overall survival months",
  ylab = "Overall survival probability"
)

png(file.path(outdir, "PIP5K1A_OS_Classical_like_KM.png"), width = 2200, height = 1800, res = 300)
print(classical_os_plot)
dev.off()

#================
# Log-rank test
#================

sd_basal_os <- survdiff(Surv(OS_time, OS_event) ~ expr_group, data = basal_os)
p_basal_os <- 1 - pchisq(sd_basal_os$chisq, df = length(sd_basal_os$n) - 1)

sd_classical_os <- survdiff(Surv(OS_time, OS_event) ~ expr_group, data = classical_os)
p_classical_os <- 1 - pchisq(sd_classical_os$chisq, df = length(sd_classical_os$n) - 1)

logrank_os_results <- data.frame(
  subtype = c("Basal-like", "Classical-like"),
  chisq = c(sd_basal_os$chisq, sd_classical_os$chisq),
  p_value = c(p_basal_os, p_classical_os)
)

write.csv(logrank_os_results, file.path(outdir, "PIP5K1A_OS_logrank_results.csv"), row.names = FALSE)

#=================
# Cox regression
#=================

cox_basal_os <- coxph(Surv(OS_time, OS_event) ~ expr_group, data = basal_os)
sum_basal_os <- summary(cox_basal_os)

cox_classical_os <- coxph(Surv(OS_time, OS_event) ~ expr_group, data = classical_os)
sum_classical_os <- summary(cox_classical_os)

cox_os_results <- data.frame(
  subtype = c("Basal-like", "Classical-like"),
  HR = c(
    sum_basal_os$conf.int[1, "exp(coef)"],
    sum_classical_os$conf.int[1, "exp(coef)"]
  ),
  lower_95_CI = c(
    sum_basal_os$conf.int[1, "lower .95"],
    sum_classical_os$conf.int[1, "lower .95"]
  ),
  upper_95_CI = c(
    sum_basal_os$conf.int[1, "upper .95"],
    sum_classical_os$conf.int[1, "upper .95"]
  ),
  p_value = c(
    sum_basal_os$coefficients[1, "Pr(>|z|)"],
    sum_classical_os$coefficients[1, "Pr(>|z|)"]
  )
)

write.csv(cox_os_results, file.path(outdir, "PIP5K1A_OS_cox_results.csv"), row.names = FALSE)

#====================
# Interaction test
#====================

merged_os_df$subtype <- factor(merged_os_df$subtype, levels = c("Basal-like", "Classical-like"))
merged_os_df$PIP5K1A_z <- scale(merged_os_df$PIP5K1A)

cox_int_os <- coxph(Surv(OS_time, OS_event) ~ PIP5K1A_z * subtype, data = merged_os_df)
sum_int_os <- summary(cox_int_os)

interaction_os_results <- data.frame(
  term = rownames(sum_int_os$coefficients),
  coef = sum_int_os$coefficients[, "coef"],
  HR = sum_int_os$coefficients[, "exp(coef)"],
  lower_95_CI = sum_int_os$conf.int[, "lower .95"],
  upper_95_CI = sum_int_os$conf.int[, "upper .95"],
  p_value = sum_int_os$coefficients[, "Pr(>|z|)"]
)

interaction_os_results$signif <- cut(
  interaction_os_results$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns"),
  right = FALSE
)

write.csv(interaction_os_results, file.path(outdir, "PIP5K1A_OS_interaction_results.csv"), row.names = FALSE)


# ===================
# KRAS MUTANT ONLY
# ==================
=
# Set KRAS mutation column name here
kras_col <- "paper_KRAS.Mutated..1.or.0."   # <-- CHANGE THIS

# Check KRAS column values first
table(meta_paad2[[kras_col]], useNA = "ifany")

# Keep only KRAS-mutant PAAD cases with subtype labels
meta_paad_kras <- meta_paad2[meta_paad2[[kras_col]] %in% c("Mutant", "Mutated", "Yes", 1), ]

# Match expression matrix columns to KRAS-mutant metadata
vst_paad_kras <- vst_mat[, rownames(meta_paad_kras), drop = FALSE]

# Split KRAS-mutant samples by subtype
basal_samples_kras <- rownames(meta_paad_kras)[meta_paad_kras[[sub_col]] == 1]
classical_samples_kras <- rownames(meta_paad_kras)[meta_paad_kras[[sub_col]] == 2]

count_basal_kras <- vst_paad_kras[, basal_samples_kras, drop = FALSE]
count_classical_kras <- vst_paad_kras[, classical_samples_kras, drop = FALSE]

dim(count_basal_kras)
dim(count_classical_kras)
table(meta_paad_kras[[sub_col]])

#=========================================
# Extract PIP5K1A expression in KRAS cases
#=========================================

gene_map <- data.frame(
  ensembl_id = rownames(gene_annot),
  gene_name = gene_annot$gene_name,
  stringsAsFactors = FALSE
)

pip5k1a_rows <- gene_map$ensembl_id[gene_map$gene_name == "PIP5K1A"]

pip5k1a_expr_basal_kras <- count_basal_kras[pip5k1a_rows, , drop = FALSE]
pip5k1a_expr_classical_kras <- count_classical_kras[pip5k1a_rows, , drop = FALSE]

rownames(pip5k1a_expr_basal_kras) <- "PIP5K1A"
rownames(pip5k1a_expr_classical_kras) <- "PIP5K1A"

#=========================================
# High/Low grouping within KRAS cases
#=========================================

pip5k1a_basal_kras_df <- as.data.frame(t(pip5k1a_expr_basal_kras))
pip5k1a_basal_kras_df$sample <- rownames(pip5k1a_basal_kras_df)

pip5k1a_classical_kras_df <- as.data.frame(t(pip5k1a_expr_classical_kras))
pip5k1a_classical_kras_df$sample <- rownames(pip5k1a_classical_kras_df)

pip5k1a_basal_kras_df$subtype <- "Basal-like"
pip5k1a_classical_kras_df$subtype <- "Classical-like"

basal_median_kras <- median(pip5k1a_basal_kras_df$PIP5K1A, na.rm = TRUE)
classical_median_kras <- median(pip5k1a_classical_kras_df$PIP5K1A, na.rm = TRUE)

pip5k1a_basal_kras_df$expr_group <- ifelse(pip5k1a_basal_kras_df$PIP5K1A >= basal_median_kras, "High", "Low")
pip5k1a_classical_kras_df$expr_group <- ifelse(pip5k1a_classical_kras_df$PIP5K1A >= classical_median_kras, "High", "Low")

pip5k1a_basal_kras_df$expr_group <- factor(pip5k1a_basal_kras_df$expr_group, levels = c("Low", "High"))
pip5k1a_classical_kras_df$expr_group <- factor(pip5k1a_classical_kras_df$expr_group, levels = c("Low", "High"))

# shorten to patient-level TCGA barcode
pip5k1a_basal_kras_df$sample <- substr(pip5k1a_basal_kras_df$sample, 1, 12)
pip5k1a_classical_kras_df$sample <- substr(pip5k1a_classical_kras_df$sample, 1, 12)

pip5k1a_kras_all <- rbind(pip5k1a_basal_kras_df, pip5k1a_classical_kras_df)
write.csv(pip5k1a_kras_all, file.path(outdir, "PIP5K1A_KRAS_basal_classical_high_low.csv"), row.names = FALSE)

#===========================
# KRAS-mutant DFS analysis
#============================

clin_dfs_kras <- clin[!is.na(clin$`Disease Free (Months)`) & !is.na(clin$`Disease Free Status`), ]

merged_dfs_kras <- merge(
  pip5k1a_kras_all,
  clin_dfs_kras,
  by.x = "sample",
  by.y = "Patient ID"
)

merged_dfs_kras$DFS_time <- as.numeric(merged_dfs_kras$`Disease Free (Months)`)
merged_dfs_kras$DFS_event <- ifelse(grepl("^1", merged_dfs_kras$`Disease Free Status`), 1, 0)

# Basal-like DFS
basal_dfs_kras <- merged_dfs_kras[merged_dfs_kras$subtype == "Basal-like", ]
basal_dfs_kras <- basal_dfs_kras[!is.na(basal_dfs_kras$DFS_time) & !is.na(basal_dfs_kras$DFS_event) & !is.na(basal_dfs_kras$expr_group), ]

fit_basal_dfs_kras <- survfit(Surv(DFS_time, DFS_event) ~ expr_group, data = basal_dfs_kras)

basal_dfs_kras_plot <- ggsurvplot(
  fit_basal_dfs_kras,
  data = basal_dfs_kras,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  title = "PIP5K1A DFS in KRAS-mutant Basal-like subtype",
  legend.title = "PIP5K1A",
  legend.labs = c("Low", "High"),
  xlab = "Disease-free months",
  ylab = "Disease-free survival probability"
)

png(file.path(outdir, "PIP5K1A_DFS_KRAS_Basal_like_KM.png"), width = 2200, height = 1800, res = 300)
print(basal_dfs_kras_plot)
dev.off()

# Classical-like DFS
classical_dfs_kras <- merged_dfs_kras[merged_dfs_kras$subtype == "Classical-like", ]
classical_dfs_kras <- classical_dfs_kras[!is.na(classical_dfs_kras$DFS_time) & !is.na(classical_dfs_kras$DFS_event) & !is.na(classical_dfs_kras$expr_group), ]

fit_classical_dfs_kras <- survfit(Surv(DFS_time, DFS_event) ~ expr_group, data = classical_dfs_kras)

classical_dfs_kras_plot <- ggsurvplot(
  fit_classical_dfs_kras,
  data = classical_dfs_kras,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  title = "PIP5K1A DFS in KRAS-mutant Classical-like subtype",
  legend.title = "PIP5K1A",
  legend.labs = c("Low", "High"),
  xlab = "Disease-free months",
  ylab = "Disease-free survival probability"
)

png(file.path(outdir, "PIP5K1A_DFS_KRAS_Classical_like_KM.png"), width = 2200, height = 1800, res = 300)
print(classical_dfs_kras_plot)
dev.off()

# Log-rank DFS
sd_basal_dfs_kras <- survdiff(Surv(DFS_time, DFS_event) ~ expr_group, data = basal_dfs_kras)
p_basal_dfs_kras <- 1 - pchisq(sd_basal_dfs_kras$chisq, df = length(sd_basal_dfs_kras$n) - 1)

sd_classical_dfs_kras <- survdiff(Surv(DFS_time, DFS_event) ~ expr_group, data = classical_dfs_kras)
p_classical_dfs_kras <- 1 - pchisq(sd_classical_dfs_kras$chisq, df = length(sd_classical_dfs_kras$n) - 1)

logrank_dfs_kras_results <- data.frame(
  subtype = c("Basal-like", "Classical-like"),
  chisq = c(sd_basal_dfs_kras$chisq, sd_classical_dfs_kras$chisq),
  p_value = c(p_basal_dfs_kras, p_classical_dfs_kras)
)

write.csv(logrank_dfs_kras_results, file.path(outdir, "PIP5K1A_DFS_KRAS_logrank_results.csv"), row.names = FALSE)

# Cox DFS
cox_basal_dfs_kras <- coxph(Surv(DFS_time, DFS_event) ~ expr_group, data = basal_dfs_kras)
sum_basal_dfs_kras <- summary(cox_basal_dfs_kras)

cox_classical_dfs_kras <- coxph(Surv(DFS_time, DFS_event) ~ expr_group, data = classical_dfs_kras)
sum_classical_dfs_kras <- summary(cox_classical_dfs_kras)

cox_dfs_kras_results <- data.frame(
  subtype = c("Basal-like", "Classical-like"),
  HR = c(sum_basal_dfs_kras$conf.int[1, "exp(coef)"], sum_classical_dfs_kras$conf.int[1, "exp(coef)"]),
  lower_95_CI = c(sum_basal_dfs_kras$conf.int[1, "lower .95"], sum_classical_dfs_kras$conf.int[1, "lower .95"]),
  upper_95_CI = c(sum_basal_dfs_kras$conf.int[1, "upper .95"], sum_classical_dfs_kras$conf.int[1, "upper .95"]),
  p_value = c(sum_basal_dfs_kras$coefficients[1, "Pr(>|z|)"], sum_classical_dfs_kras$coefficients[1, "Pr(>|z|)"])
)

write.csv(cox_dfs_kras_results, file.path(outdir, "PIP5K1A_DFS_KRAS_cox_results.csv"), row.names = FALSE)

# Interaction DFS
merged_dfs_kras$subtype <- factor(merged_dfs_kras$subtype, levels = c("Basal-like", "Classical-like"))
merged_dfs_kras$PIP5K1A_z <- scale(merged_dfs_kras$PIP5K1A)

cox_int_dfs_kras <- coxph(Surv(DFS_time, DFS_event) ~ PIP5K1A_z * subtype, data = merged_dfs_kras)
sum_int_dfs_kras <- summary(cox_int_dfs_kras)

interaction_dfs_kras_results <- data.frame(
  term = rownames(sum_int_dfs_kras$coefficients),
  coef = sum_int_dfs_kras$coefficients[, "coef"],
  HR = sum_int_dfs_kras$coefficients[, "exp(coef)"],
  lower_95_CI = sum_int_dfs_kras$conf.int[, "lower .95"],
  upper_95_CI = sum_int_dfs_kras$conf.int[, "upper .95"],
  p_value = sum_int_dfs_kras$coefficients[, "Pr(>|z|)"]
)

interaction_dfs_kras_results$signif <- cut(
  interaction_dfs_kras_results$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns"),
  right = FALSE
)

write.csv(interaction_dfs_kras_results, file.path(outdir, "PIP5K1A_DFS_KRAS_interaction_results.csv"), row.names = FALSE)

#==========================
# KRAS-mutant OS analysis
#==========================

clin_os_kras <- clin[!is.na(clin$`Overall Survival (Months)`) & !is.na(clin$`Overall Survival Status`), ]

merged_os_kras <- merge(
  pip5k1a_kras_all,
  clin_os_kras,
  by.x = "sample",
  by.y = "Patient ID"
)

merged_os_kras$OS_time <- as.numeric(merged_os_kras$`Overall Survival (Months)`)
merged_os_kras$OS_event <- ifelse(grepl("^1", merged_os_kras$`Overall Survival Status`), 1, 0)

# Basal-like OS
basal_os_kras <- merged_os_kras[merged_os_kras$subtype == "Basal-like", ]
basal_os_kras <- basal_os_kras[!is.na(basal_os_kras$OS_time) & !is.na(basal_os_kras$OS_event) & !is.na(basal_os_kras$expr_group), ]

fit_basal_os_kras <- survfit(Surv(OS_time, OS_event) ~ expr_group, data = basal_os_kras)

basal_os_kras_plot <- ggsurvplot(
  fit_basal_os_kras,
  data = basal_os_kras,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  title = "PIP5K1A OS in KRAS-mutant Basal-like subtype",
  legend.title = "PIP5K1A",
  legend.labs = c("Low", "High"),
  xlab = "Overall survival months",
  ylab = "Overall survival probability"
)

png(file.path(outdir, "PIP5K1A_OS_KRAS_Basal_like_KM.png"), width = 2200, height = 1800, res = 300)
print(basal_os_kras_plot)
dev.off()

# Classical-like OS
classical_os_kras <- merged_os_kras[merged_os_kras$subtype == "Classical-like", ]
classical_os_kras <- classical_os_kras[!is.na(classical_os_kras$OS_time) & !is.na(classical_os_kras$OS_event) & !is.na(classical_os_kras$expr_group), ]

fit_classical_os_kras <- survfit(Surv(OS_time, OS_event) ~ expr_group, data = classical_os_kras)

classical_os_kras_plot <- ggsurvplot(
  fit_classical_os_kras,
  data = classical_os_kras,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  title = "PIP5K1A OS in KRAS-mutant Classical-like subtype",
  legend.title = "PIP5K1A",
  legend.labs = c("Low", "High"),
  xlab = "Overall survival months",
  ylab = "Overall survival probability"
)

png(file.path(outdir, "PIP5K1A_OS_KRAS_Classical_like_KM.png"), width = 2200, height = 1800, res = 300)
print(classical_os_kras_plot)
dev.off()

# Log-rank OS
sd_basal_os_kras <- survdiff(Surv(OS_time, OS_event) ~ expr_group, data = basal_os_kras)
p_basal_os_kras <- 1 - pchisq(sd_basal_os_kras$chisq, df = length(sd_basal_os_kras$n) - 1)

sd_classical_os_kras <- survdiff(Surv(OS_time, OS_event) ~ expr_group, data = classical_os_kras)
p_classical_os_kras <- 1 - pchisq(sd_classical_os_kras$chisq, df = length(sd_classical_os_kras$n) - 1)

logrank_os_kras_results <- data.frame(
  subtype = c("Basal-like", "Classical-like"),
  chisq = c(sd_basal_os_kras$chisq, sd_classical_os_kras$chisq),
  p_value = c(p_basal_os_kras, p_classical_os_kras)
)

write.csv(logrank_os_kras_results, file.path(outdir, "PIP5K1A_OS_KRAS_logrank_results.csv"), row.names = FALSE)

# Cox OS
cox_basal_os_kras <- coxph(Surv(OS_time, OS_event) ~ expr_group, data = basal_os_kras)
sum_basal_os_kras <- summary(cox_basal_os_kras)

cox_classical_os_kras <- coxph(Surv(OS_time, OS_event) ~ expr_group, data = classical_os_kras)
sum_classical_os_kras <- summary(cox_classical_os_kras)

cox_os_kras_results <- data.frame(
  subtype = c("Basal-like", "Classical-like"),
  HR = c(sum_basal_os_kras$conf.int[1, "exp(coef)"], sum_classical_os_kras$conf.int[1, "exp(coef)"]),
  lower_95_CI = c(sum_basal_os_kras$conf.int[1, "lower .95"], sum_classical_os_kras$conf.int[1, "lower .95"]),
  upper_95_CI = c(sum_basal_os_kras$conf.int[1, "upper .95"], sum_classical_os_kras$conf.int[1, "upper .95"]),
  p_value = c(sum_basal_os_kras$coefficients[1, "Pr(>|z|)"], sum_classical_os_kras$coefficients[1, "Pr(>|z|)"])
)

write.csv(cox_os_kras_results, file.path(outdir, "PIP5K1A_OS_KRAS_cox_results.csv"), row.names = FALSE)

# Interaction OS
merged_os_kras$subtype <- factor(merged_os_kras$subtype, levels = c("Basal-like", "Classical-like"))
merged_os_kras$PIP5K1A_z <- scale(merged_os_kras$PIP5K1A)

cox_int_os_kras <- coxph(Surv(OS_time, OS_event) ~ PIP5K1A_z * subtype, data = merged_os_kras)
sum_int_os_kras <- summary(cox_int_os_kras)

interaction_os_kras_results <- data.frame(
  term = rownames(sum_int_os_kras$coefficients),
  coef = sum_int_os_kras$coefficients[, "coef"],
  HR = sum_int_os_kras$coefficients[, "exp(coef)"],
  lower_95_CI = sum_int_os_kras$conf.int[, "lower .95"],
  upper_95_CI = sum_int_os_kras$conf.int[, "upper .95"],
  p_value = sum_int_os_kras$coefficients[, "Pr(>|z|)"]
)

interaction_os_kras_results$signif <- cut(
  interaction_os_kras_results$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns"),
  right = FALSE
)

write.csv(interaction_os_kras_results, file.path(outdir, "PIP5K1A_OS_KRAS_interaction_results.csv"), row.names = FALSE)