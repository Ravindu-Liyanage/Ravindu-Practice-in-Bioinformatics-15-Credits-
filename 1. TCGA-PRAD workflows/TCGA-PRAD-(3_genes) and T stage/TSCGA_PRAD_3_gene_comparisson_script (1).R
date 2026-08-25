# ==========================================================================
# Ravindu: 3 GENE COMPARISSON & FOLH1 Expression Across Pathological T Stage
# ==========================================================================

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
library(RColorBrewer)
library(grid)

outdir <- "new_output"
sepdir <- file.path(outdir, "separate_plots")
kdir <- file.path(outdir, "KM_plots")
hmdir <- file.path(outdir, "heatmaps")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(sepdir, showWarnings = FALSE, recursive = TRUE)
dir.create(kdir, showWarnings = FALSE, recursive = TRUE)
dir.create(hmdir, showWarnings = FALSE, recursive = TRUE)

# collect Excel sheets here, then write all workbooks at the very end
gene_sheets <- list()
tstage_sheets <- list()

# THEME
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

nature_cols <- c(
  "#4DBBD5FF",
  "#E64B35FF",
  "#00A087FF",
  "#3C5488FF",
  "#F39B7FFF",
  "#8491B4FF",
  "#91D1C2FF",
  "#DC0000FF"
)

# ===================
# STEP 1: LOAD DATA
# ===================

data <- readRDS("TCGA_panCancer_FOLH1.rds")

meta_all <- as.data.frame(colData(data))
meta_all$barcode <- rownames(meta_all)

# SUBSET TO TCGA-PRAD EARLY
prad_keep <- as.character(meta_all$project_id) == "TCGA-PRAD"

data_prad <- data[, prad_keep]
meta_prad <- as.data.frame(colData(data_prad))
meta_prad$barcode <- rownames(meta_prad)

count_mat <- assay(data_prad)
gene_annot <- as.data.frame(rowData(data_prad))
target_genes <- c("FCGR3A", "FOLH1", "PTGS2")

dim(count_mat)
dim(meta_prad)
dim(gene_annot)

# ===========================================================
# STEP 2: NORMALIZATION + VARIANCE STABILIZATION ON PRAD ONLY
# ===========================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta_prad,
  design = ~ 1
)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- estimateSizeFactors(dds)

norm_mat <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

meta_prad <- meta_prad[colnames(expr_mat), , drop = FALSE]

# GENE ANNOTATION
if ("gene_name" %in% colnames(gene_annot)) {
  gene_name_col <- "gene_name"
} else {
  stop("gene_name column not found in gene annotation.")
}

idx <- match(target_genes, gene_annot[[gene_name_col]])
missing_genes <- target_genes[is.na(idx)]

expr_genes <- expr_mat[idx[!is.na(idx)], , drop = FALSE]
rownames(expr_genes) <- target_genes[!is.na(idx)]

folh1_idx <- which(rownames(expr_mat) == rownames(expr_mat)[idx[which(target_genes == "FOLH1")]][1])
folh1_idx <- folh1_idx[1]

# DEFINE TCGA-PRAD SAMPLE GROUPS
meta_prad$group <- case_when(
  substr(meta_prad$barcode, 14, 15) == "01" ~ "Primary Tumor",
  substr(meta_prad$barcode, 14, 15) == "11" ~ "Normal",
  substr(meta_prad$barcode, 14, 15) == "06" ~ "Metastatic",
  TRUE ~ NA_character_
)

print(table(meta_prad$group, useNA = "ifany"))

group_tbl <- as.data.frame(table(meta_prad$group))
colnames(group_tbl) <- c("Group", "Count")
gene_sheets[["PRAD_sample_groups"]] <- group_tbl

common_samples <- intersect(meta_prad$barcode, colnames(expr_genes))
meta_prad <- meta_prad[match(common_samples, meta_prad$barcode), , drop = FALSE]
expr_prad <- expr_genes[, common_samples, drop = FALSE]

# =================================
# STEP 3: PRIMARY TUMOR VS NORMAL
# =================================

meta_tn <- meta_prad %>%
  filter(group %in% c("Primary Tumor", "Normal"))

expr_tn <- expr_prad[, meta_tn$barcode, drop = FALSE]

df_long <- as.data.frame(expr_tn)
df_long$gene <- rownames(df_long)

df_long <- df_long %>%
  pivot_longer(cols = -gene, names_to = "barcode", values_to = "expression") %>%
  left_join(meta_tn[, c("barcode", "group")], by = "barcode") %>%
  filter(!is.na(group)) %>%
  mutate(expression = as.numeric(expression))

stats_tbl <- df_long %>%
  group_by(gene) %>%
  summarise(
    p.value = t.test(expression ~ group)$p.value,
    t.statistic = unname(t.test(expression ~ group)$statistic),
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

gene_sheets[["PRAD_gene_stats"]] <- stats_tbl

for (g in rownames(expr_prad)) {
  p <- df_long %>%
    filter(gene == g) %>%
    ggplot(aes(x = group, y = expression, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6) +
    stat_compare_means(method = "t.test") +
    labs(
      title = g,
      y = "VST expression",
      x = NULL
    ) +
    theme_nature()
  
  ggsave(file.path(sepdir, paste0(g, "_PRAD.png")), p, width = 6, height = 5)
}

p_all <- ggplot(df_long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_wrap(~gene) +
  stat_compare_means(method = "t.test") +
  labs(y = "VST expression", x = NULL) +
  theme_nature()

ggsave(file.path(outdir, "PRAD_3genes_combined_new.png"), p_all, width = 10, height = 5)

# =========================
# STEP 4: FOLH1 BY GLEASON SCORE
# =========================

gleason_df <- meta_prad

if (!("gleason_score" %in% colnames(gleason_df))) {
  stop("gleason_score column not found in PRAD metadata.")
}

# Convert possible character/factor values to numeric.
# This also handles values such as "6", "7", "8", "9", and "10".
gleason_df <- gleason_df %>%
  mutate(
    Gleason_numeric = suppressWarnings(
      as.numeric(as.character(gleason_score))
    )
  ) %>%
  filter(
    group == "Primary Tumor",
    !is.na(Gleason_numeric)
  ) %>%
  mutate(
    Gleason_group = case_when(
      Gleason_numeric <= 6 ~ "6",
      Gleason_numeric == 7 ~ "7",
      Gleason_numeric >= 8 ~ "8-10",
      TRUE ~ NA_character_
    ),
    Gleason_group = factor(
      Gleason_group,
      levels = c("6", "7", "8-10"),
      ordered = TRUE
    )
  ) %>%
  filter(!is.na(Gleason_group))

# Match metadata and expression samples by barcode.
common_gleason <- intersect(
  gleason_df$barcode,
  colnames(expr_prad)
)

gleason_df <- gleason_df[
  match(common_gleason, gleason_df$barcode),
  ,
  drop = FALSE
]

expr_gleason <- expr_prad[
  "FOLH1",
  common_gleason,
  drop = FALSE
]

# Long-format FOLH1 data.
folh1_gleason <- data.frame(
  barcode = colnames(expr_gleason),
  FOLH1 = as.numeric(expr_gleason[1, ]),
  stringsAsFactors = FALSE
) %>%
  left_join(
    gleason_df[, c("barcode", "Gleason_numeric", "Gleason_group")],
    by = "barcode"
  ) %>%
  filter(
    !is.na(FOLH1),
    !is.na(Gleason_group)
  )

# Check group sizes before testing.
print(table(folh1_gleason$Gleason_group, useNA = "ifany"))

# Summary statistics.
folh1_gleason_summary <- folh1_gleason %>%
  group_by(Gleason_group) %>%
  summarise(
    n_samples = n(),
    mean_expr = mean(FOLH1, na.rm = TRUE),
    median_expr = median(FOLH1, na.rm = TRUE),
    sd_expr = sd(FOLH1, na.rm = TRUE),
    sem_expr = sd_expr / sqrt(n_samples),
    min_expr = min(FOLH1, na.rm = TRUE),
    max_expr = max(FOLH1, na.rm = TRUE),
    .groups = "drop"
  )

gene_sheets[["FOLH1_Gleason_summary"]] <-
  folh1_gleason_summary

# Omnibus non-parametric test.
kw_gleason <- kruskal.test(
  FOLH1 ~ Gleason_group,
  data = folh1_gleason
)

kw_gleason_table <- data.frame(
  Test = "Kruskal-Wallis",
  Statistic = unname(kw_gleason$statistic),
  df = unname(kw_gleason$parameter),
  p.value = kw_gleason$p.value
)

gene_sheets[["FOLH1_Gleason_KW"]] <-
  kw_gleason_table

# Pairwise Wilcoxon tests with Benjamini-Hochberg correction.
pairwise_gleason <- pairwise.wilcox.test(
  x = folh1_gleason$FOLH1,
  g = folh1_gleason$Gleason_group,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_gleason_table <- as.data.frame(
  as.table(pairwise_gleason$p.value)
)

colnames(pairwise_gleason_table) <- c(
  "Gleason_group_1",
  "Gleason_group_2",
  "p_adj_BH"
)

pairwise_gleason_table <- pairwise_gleason_table %>%
  filter(!is.na(p_adj_BH)) %>%
  mutate(
    significance = case_when(
      p_adj_BH <= 0.001 ~ "***",
      p_adj_BH <= 0.01 ~ "**",
      p_adj_BH <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj_BH)

gene_sheets[["FOLH1_Gleason_pairwise"]] <-
  pairwise_gleason_table

# Optional ordered trend test.
# This tests whether expression changes monotonically with Gleason score.
trend_test <- cor.test(
  folh1_gleason$FOLH1,
  folh1_gleason$Gleason_numeric,
  method = "spearman",
  exact = FALSE
)

trend_test_table <- data.frame(
  Test = "Spearman trend test",
  rho = unname(trend_test$estimate),
  p.value = trend_test$p.value
)

gene_sheets[["FOLH1_Gleason_trend"]] <-
  trend_test_table

cat(
  "FOLH1 Gleason Kruskal-Wallis p-value:",
  format.pval(kw_gleason$p.value),
  "\n"
)

# Plot specifically for FOLH1.
p_folh1_gleason <- ggplot(
  folh1_gleason,
  aes(
    x = Gleason_group,
    y = FOLH1,
    fill = Gleason_group
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.45,
    size = 1.3
  ) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    label.y = max(folh1_gleason$FOLH1, na.rm = TRUE) * 1.08
  ) +
  scale_fill_manual(
    values = c(
      "6" = "#4DBBD5FF",
      "7" = "#00A087FF",
      "8-10" = "#E64B35FF"
    )
  ) +
  labs(
    title = "FOLH1 expression across Gleason groups in TCGA-PRAD",
    subtitle = paste0(
      "Kruskal-Wallis p = ",
      format.pval(kw_gleason$p.value, digits = 3)
    ),
    x = "Gleason group",
    y = "FOLH1 VST expression"
  ) +
  theme_nature() +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(
      hjust = 0.5,
      colour = "black"
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "FOLH1_Gleason_expression.png"
  ),
  plot = p_folh1_gleason,
  width = 7,
  height = 5.5,
  dpi = 300
)

# =========================
# STEP 5: SURVIVAL ANALYSIS
# =========================
clin_prad <- meta_prad
clin_prad$patient_id <- substr(clin_prad$barcode, 1, 12)

clin_surv <- clin_prad %>%
  group_by(patient_id) %>%
  summarise(
    time = suppressWarnings(max(as.numeric(c(days_to_death, days_to_last_follow_up)), na.rm = TRUE)),
    event = ifelse(any(vital_status == "Dead", na.rm = TRUE), 1, 0),
    .groups = "drop"
  ) %>%
  mutate(time = ifelse(is.infinite(time), NA, time)) %>%
  filter(!is.na(time) & time > 0)

expr_df_surv <- as.data.frame(t(expr_prad))
expr_df_surv$barcode <- rownames(expr_df_surv)
expr_df_surv$patient_id <- substr(expr_df_surv$barcode, 1, 12)

km_df <- expr_df_surv %>%
  left_join(clin_surv, by = "patient_id") %>%
  filter(!is.na(time), !is.na(event))

genes <- rownames(expr_prad)

logrank_tbl <- lapply(genes, function(gene) {
  tmp <- km_df %>%
    select(time, event, all_of(gene)) %>%
    rename(expression = all_of(gene)) %>%
    filter(!is.na(expression)) %>%
    mutate(group = ifelse(expression >= median(expression, na.rm = TRUE), "High", "Low"))
  
  sd <- survdiff(Surv(time, event) ~ group, data = tmp)
  p <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  
  sig <- ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns")))
  
  data.frame(
    Gene = gene,
    p.value = p,
    Significance = sig
  )
})

logrank_tbl <- do.call(rbind, logrank_tbl)
gene_sheets[["logrank_table"]] <- logrank_tbl

plot_list <- list()

for (gene in genes) {
  tmp <- km_df %>%
    select(time, event, all_of(gene)) %>%
    rename(expression = all_of(gene)) %>%
    mutate(group = ifelse(expression >= median(expression, na.rm = TRUE), "High", "Low")) %>%
    filter(!is.na(time), !is.na(event), !is.na(expression))
  
  tmp$group <- factor(tmp$group, levels = c("Low", "High"))
  fit <- survfit(Surv(time, event) ~ group, data = tmp)
  
  p <- ggsurvplot(
    fit,
    data = tmp,
    pval = TRUE,
    risk.table = TRUE,
    title = gene,
    ggtheme = theme_nature()
  )
  
  plot_list[[gene]] <- p$plot
  
  ggsave(
    filename = file.path(kdir, paste0(gene, "_KM.png")),
    plot = p$plot,
    width = 6,
    height = 5
  )
}

combined <- ggarrange(
  plotlist = plot_list,
  ncol = 3,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  filename = file.path(outdir, "KM_three_genes_combined.png"),
  plot = combined,
  width = 14,
  height = 5
)

# ===============================
# STEP 6: CO-EXPRESSION HEATMAPS
# ===============================
gene_map <- gene_annot %>%
  select(gene_id, all_of(gene_name_col)) %>%
  mutate(gene_id = sub("\\.[0-9]+$", "", gene_id)) %>%
  distinct()

colnames(gene_map)[2] <- "gene_name"

expr_df_heat <- as.data.frame(expr_mat)
expr_df_heat$gene_id <- sub("\\.[0-9]+$", "", rownames(expr_df_heat))

expr_df_heat <- expr_df_heat %>%
  left_join(gene_map, by = "gene_id") %>%
  filter(!is.na(gene_name)) %>%
  select(gene_name, everything(), -gene_id) %>%
  distinct(gene_name, .keep_all = TRUE)

expr_sym <- as.matrix(expr_df_heat[, -1])
rownames(expr_sym) <- expr_df_heat$gene_name
storage.mode(expr_sym) <- "numeric"

make_coexpression_heatmap <- function(target_gene, expr_sym, top_n = 15, out_file = NULL) {
  if (!target_gene %in% rownames(expr_sym)) return(NULL)
  
  target_vec <- expr_sym[target_gene, ]
  
  cors <- apply(
    expr_sym, 1,
    function(x) cor(x, target_vec, method = "pearson", use = "pairwise.complete.obs")
  )
  
  cors <- cors[!is.na(cors)]
  cors <- sort(cors, decreasing = TRUE)
  
  pos_genes <- names(head(cors[cors > 0], top_n))
  neg_genes <- names(head(sort(cors[cors < 0], decreasing = FALSE), top_n))
  
  genes_keep <- unique(c(target_gene, pos_genes, neg_genes))
  genes_keep <- genes_keep[genes_keep %in% rownames(expr_sym)]
  
  sub_mat <- expr_sym[genes_keep, , drop = FALSE]
  cor_mat <- cor(t(sub_mat), use = "pairwise.complete.obs")
  
  ann_row <- data.frame(
    Type = ifelse(
      rownames(cor_mat) == target_gene, "Target",
      ifelse(rownames(cor_mat) %in% pos_genes, "Positive", "Negative")
    )
  )
  rownames(ann_row) <- rownames(cor_mat)
  
  ann_col <- ann_row
  rownames(ann_col) <- colnames(cor_mat)
  
  ann_colors <- list(
    Type = c(
      Target = "#000000",
      Positive = "#D55E00",
      Negative = "#0072B2"
    )
  )
  
  hm_cols <- colorRampPalette(c("#313695", "#F7F7F7", "#A50026"))(120)
  
  pheatmap(
    cor_mat,
    color = hm_cols,
    breaks = seq(-1, 1, length.out = 121),
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    treeheight_row = 20,
    treeheight_col = 20,
    angle_col = 45,
    fontsize_row = 9,
    fontsize_col = 9,
    fontsize = 10,
    main = paste0(target_gene, " co-expression"),
    legend = TRUE,
    annotation_row = ann_row,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_rownames = TRUE,
    show_colnames = TRUE,
    na_col = "#EEEEEE",
    silent = TRUE
  ) -> ph
  
  if (!is.null(out_file)) {
    png(out_file, width = 2200, height = 2200, res = 300, type = "cairo-png")
    grid::grid.newpage()
    grid::grid.draw(ph$gtable)
    dev.off()
  }
  
  invisible(ph)
}

for (g in genes) {
  if (g %in% rownames(expr_sym)) {
    make_coexpression_heatmap(
      target_gene = g,
      expr_sym = expr_sym,
      top_n = 15,
      out_file = file.path(hmdir, paste0(g, "_heatmap.png"))
    )
  }
}

# ==========================
# STEP 7: T-STAGE ANALYSIS
# ==========================

stage_candidates <- c("ajcc_pathologic_t", "pathologic_T", "t_stage", "T_stage", "tumor_stage")
stage_col <- stage_candidates[stage_candidates %in% colnames(meta_prad)][1]

if (is.na(stage_col)) {
  stop("No T-stage column found in meta_prad. Check colnames(meta_prad).")
}

cat("Using T-stage column:", stage_col, "\n")
print(table(meta_prad[[stage_col]], useNA = "ifany"))

meta_tstage <- meta_prad %>%
  mutate(
    T_stage_raw = as.character(.data[[stage_col]]),
    T_stage = case_when(
      grepl("^T1", T_stage_raw, ignore.case = TRUE) ~ "T1",
      grepl("^T2", T_stage_raw, ignore.case = TRUE) ~ "T2",
      grepl("^T3", T_stage_raw, ignore.case = TRUE) ~ "T3",
      grepl("^T4", T_stage_raw, ignore.case = TRUE) ~ "T4",
      TRUE ~ NA_character_
    ),
    T_stage = factor(T_stage, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE)
  ) %>%
  filter(!is.na(T_stage))

tstage_samples <- intersect(meta_tstage$barcode, colnames(expr_mat))
meta_tstage <- meta_tstage[match(tstage_samples, meta_tstage$barcode), , drop = FALSE]
expr_tstage <- expr_mat[, tstage_samples, drop = FALSE]

folh1_expr_tstage <- data.frame(
  barcode = colnames(expr_tstage),
  FOLH1 = as.numeric(expr_tstage[folh1_idx, ]),
  stringsAsFactors = FALSE
) %>%
  left_join(
    meta_tstage[, c("barcode", "T_stage")],
    by = "barcode"
  ) %>%
  filter(!is.na(FOLH1), !is.na(T_stage))

tstage_summary <- folh1_expr_tstage %>%
  group_by(T_stage) %>%
  summarise(
    n_samples = n(),
    mean_expr = mean(FOLH1, na.rm = TRUE),
    median_expr = median(FOLH1, na.rm = TRUE),
    sd_expr = sd(FOLH1, na.rm = TRUE),
    .groups = "drop"
  )

tstage_sheets[["FOLH1_Tstage_summary"]] <- tstage_summary

kw_tstage <- kruskal.test(FOLH1 ~ T_stage, data = folh1_expr_tstage)

pairwise_tstage <- pairwise.wilcox.test(
  x = folh1_expr_tstage$FOLH1,
  g = folh1_expr_tstage$T_stage,
  p.adjust.method = "BH",
  exact = FALSE
)

pairwise_tstage_df <- as.data.frame(as.table(pairwise_tstage$p.value))
colnames(pairwise_tstage_df) <- c("T_stage_1", "T_stage_2", "p_adj_bh")
pairwise_tstage_df <- pairwise_tstage_df %>%
  filter(!is.na(p_adj_bh)) %>%
  arrange(p_adj_bh)

tstage_sheets[["FOLH1_Tstage_pairwise_wilcox"]] <- pairwise_tstage_df

cat("Kruskal-Wallis p-value for T-stage:", format(kw_tstage$p.value, scientific = TRUE), "\n")

p_tstage <- ggplot(folh1_expr_tstage, aes(x = T_stage, y = FOLH1, fill = T_stage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.5, size = 1) +
  labs(
    title = "FOLH1 expression across T stage in TCGA-PRAD",
    x = "T stage",
    y = "FOLH1 VST expression"
  ) +
  theme_nature() +
  theme(legend.position = "none")

ggsave(
  file.path(outdir, "FOLH1_Tstage_boxplot.png"),
  p_tstage, width = 7, height = 5, dpi = 300
)

# =========================================
# FINAL STEP: WRITE ONLY 2 EXCEL WORKBOOKS
# =========================================

safe_sheet_name <- function(x, existing = character()) {
  x <- gsub("[\\\\/:*?\\[\\]]", "_", x)
  x <- substr(x, 1, 31)
  if (!(x %in% existing)) return(x)
  
  base <- substr(x, 1, 28)
  i <- 1
  new_x <- paste0(base, "_", i)
  while (new_x %in% existing) {
    i <- i + 1
    new_x <- paste0(base, "_", i)
  }
  substr(new_x, 1, 31)
}

# Workbook 1: 3-gene comparison related outputs
wb_gene <- createWorkbook()

existing_names <- character()
for (nm in names(gene_sheets)) {
  sn <- safe_sheet_name(nm, existing_names)
  addWorksheet(wb_gene, sn)
  writeData(wb_gene, sheet = sn, x = gene_sheets[[nm]])
  existing_names <- c(existing_names, sn)
}

saveWorkbook(
  wb_gene,
  file = file.path(outdir, "three_gene_comparison.xlsx"),
  overwrite = TRUE
)

# Workbook 2: T-stage outputs
wb_tstage <- createWorkbook()

existing_names <- character()
for (nm in names(tstage_sheets)) {
  sn <- safe_sheet_name(nm, existing_names)
  addWorksheet(wb_tstage, sn)
  writeData(wb_tstage, sheet = sn, x = tstage_sheets[[nm]])
  existing_names <- c(existing_names, sn)
}

saveWorkbook(
  wb_tstage,
  file = file.path(outdir, "tstage_comparison.xlsx"),
  overwrite = TRUE
)

cat("Excel workbooks created successfully:\n")
cat(" -", file.path(outdir, "three_gene_comparison.xlsx"), "\n")
cat(" -", file.path(outdir, "tstage_comparison.xlsx"), "\n")