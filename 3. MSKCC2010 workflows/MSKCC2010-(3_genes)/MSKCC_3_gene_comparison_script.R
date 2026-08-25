library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)

#=======================
# NAture themes
#========================

library(ggplot2)

theme_nature <- function(base_size = 14, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = base_size + 2
      ),
      axis.title = element_text(
        face = "bold",
        colour = "black"
      ),
      axis.text = element_text(
        colour = "black"
      ),
      axis.line = element_line(
        linewidth = 0.7,
        colour = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.6,
        colour = "black"
      ),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      strip.background = element_rect(
        fill = "white",
        colour = "black"
      ),
      strip.text = element_text(face = "bold")
    )
}

theme_set(theme_nature())

nature_cols <- c(
  "#4DBBD5FF",  # Blue
  "#E64B35FF",  # Red
  "#00A087FF",  # Green
  "#3C5488FF",  # Navy
  "#F39B7FFF",  # Orange
  "#8491B4FF",  # Purple
  "#91D1C2FF",  # Teal
  "#DC0000FF"   # Dark red
)



# =========================================================
# OUTPUT DIRECTORIES
# =========================================================

outdir <- "output_MSKCC2010"
sepdir <- file.path(outdir, "separate_plots")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(sepdir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# LOAD DATA
# =========================================================

data <- readRDS("MSKCC2010.rds")

rna <- experiments(data)[["mrna_agilent_microarray"]]

expr_mat <- assay(rna)
rownames(expr_mat) <- as.character(rownames(expr_mat))

# =========================================================
# METADATA (IMPORTANT: ONLY USE rna-level metadata)
# =========================================================

meta <- as.data.frame(colData(data))   # NOT colData(rna)
meta$sample_id <- rownames(meta)

# =========================================================
# GENE EXTRACTION
# =========================================================

target_entrez <- c(
  FCGR3A = "2214",
  FOLH1  = "2346",
  PTGS2  = "5743"
)

available_genes <- target_entrez[target_entrez %in% rownames(expr_mat)]

expr_genes <- expr_mat[available_genes, , drop = FALSE]

rownames(expr_genes) <- names(target_entrez)[
  names(target_entrez) %in% names(available_genes)
]

# =========================================================
# GROUP ASSIGNMENT (DO NOT FILTER YET)
# =========================================================

meta$group <- case_when(
  grepl("primary", meta$SAMPLE_TYPE, ignore.case = TRUE) ~ "Primary Tumor",
  grepl("met", meta$SAMPLE_TYPE, ignore.case = TRUE) ~ "Metastatic Tumor",
  TRUE ~ NA_character_
)

# =========================================================
# ALIGN SAMPLES (CRITICAL STEP)
# =========================================================

common_samples <- intersect(meta$sample_id, colnames(expr_genes))

stopifnot(length(common_samples) > 0)

expr <- expr_genes[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$sample_id), ]

# remove unmatched metadata rows safely
meta <- meta %>% filter(!is.na(group))

# =========================================================
# LONG FORMAT
# =========================================================

df_long <- as.data.frame(expr)
df_long$gene <- rownames(df_long)

df_long <- df_long %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  left_join(meta[, c("sample_id", "group")], by = "sample_id") %>%
  filter(!is.na(group)) %>%
  mutate(expression = as.numeric(expression))

# =========================================================
# STATISTICS
# =========================================================

stats_tbl <- df_long %>%
  group_by(gene) %>%
  summarise(
    p.value = t.test(expression ~ group)$p.value,
    statistic = t.test(expression ~ group)$statistic,
    .groups = "drop"
  )

write.xlsx(stats_tbl, file.path(outdir, "MSKCC_gene_stats.xlsx"))

# =========================================================
# INDIVIDUAL PLOTS
# =========================================================

for (g in unique(df_long$gene)) {
  
  p <- df_long %>%
    filter(gene == g) %>%
    ggplot(aes(x = group, y = expression, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6) +
    stat_compare_means(method = "t.test") +
    labs(title = g, x = NULL, y = "Expression")
  
  ggsave(
    file.path(sepdir, paste0(g, "_MSKCC.png")),
    p,
    width = 6,
    height = 5
  )
}

# =========================================================
# COMBINED PLOT
# =========================================================

p_all <- ggplot(df_long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_wrap(~gene) +
  stat_compare_means(method = "t.test")

ggsave(
  file.path(outdir, "MSKCC_3genes_combined.png"),
  p_all,
  width = 10,
  height = 5
)

# =========================================================
# GROUP COUNTS
# =========================================================

group_tbl <- as.data.frame(table(meta$group))
colnames(group_tbl) <- c("Group", "Count")

write.xlsx(group_tbl,
           file.path(outdir, "MSKCC_sample_groups.xlsx"))

# =========================================================
# GLEASON ANALYSIS (PRIMARY TUMORS ONLY)
# =========================================================

# ---------------------------------------------------------
# 1. KEEP PRIMARY TUMORS ONLY
# ---------------------------------------------------------

meta_gleason <- meta %>%
  filter(group == "Primary Tumor")

expr_gleason <- expr[, meta_gleason$sample_id, drop = FALSE]

# ---------------------------------------------------------
# 2. CREATE GLEASON GROUPS
# ---------------------------------------------------------
# IMPORTANT: adjust column name if needed

meta_gleason$gleason_group <- case_when(
  meta_gleason$GLEASON_SCORE %in% c("6", "3+3") ~ "Low",
  meta_gleason$GLEASON_SCORE %in% c("7", "3+4", "4+3") ~ "Intermediate",
  meta_gleason$GLEASON_SCORE %in% c("8", "9", "10", "4+4", "4+5", "5+4", "5+5") ~ "High",
  TRUE ~ NA_character_
)

# remove missing Gleason
meta_gleason <- meta_gleason %>%
  filter(!is.na(gleason_group))

# re-align expression AFTER filtering
expr_gleason <- expr_gleason[, meta_gleason$sample_id, drop = FALSE]

# ---------------------------------------------------------
# 3. LONG FORMAT
# ---------------------------------------------------------

df_gleason <- as.data.frame(expr_gleason)
df_gleason$gene <- rownames(df_gleason)

df_gleason <- df_gleason %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression"
  ) %>%
  left_join(meta_gleason[, c("sample_id", "gleason_group")],
            by = "sample_id") %>%
  filter(!is.na(gleason_group)) %>%
  mutate(expression = as.numeric(expression))

df_gleason$gleason_group <- factor(
  df_gleason$gleason_group,
  levels = c("Low", "Intermediate", "High")
)

# ---------------------------------------------------------
# 4. ANOVA PER GENE
# ---------------------------------------------------------

library(dplyr)
library(openxlsx)

anova_tbl <- df_gleason %>%
  group_by(gene) %>%
  group_modify(~{
    fit <- aov(expression ~ gleason_group, data = .x)
    s <- summary(fit)[[1]]
    
    ss_between <- s[["Sum Sq"]][1]
    ss_within  <- s[["Sum Sq"]][2]
    ss_total   <- ss_between + ss_within
    
    df_between <- s[["Df"]][1]
    df_within  <- s[["Df"]][2]
    df_total   <- df_between + df_within
    
    ms_between <- s[["Mean Sq"]][1]
    ms_within  <- s[["Mean Sq"]][2]
    
    f_val <- s[["F value"]][1]
    p_val <- s[["Pr(>F)"]][1]
    
    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*", "ns")))
    
    tibble(
      Source = c("Between", "Within", "Total"),
      Sum_of_Squares = c(ss_between, ss_within, ss_total),
      df = c(df_between, df_within, df_total),
      Mean_Square = c(ms_between, ms_within, NA_real_),
      F = c(f_val, NA_real_, NA_real_),
      Sig = c(sig, "", "")
    )
  }) %>%
  ungroup() %>%
  rename(Gene = gene)

write.xlsx(anova_tbl, file.path(outdir, "MSKCC_gleason_anova.xlsx"))

# ---------------------------------------------------------
# 5. PLOT
# ---------------------------------------------------------

p_gleason <- ggplot(df_gleason,
                    aes(x = gleason_group,
                        y = expression,
                        fill = gleason_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  facet_wrap(~gene) +
  stat_compare_means(method = "anova") +
  labs(title = "MSKCC2010: Gleason Score vs Gene Expression")

ggsave(file.path(outdir, "MSKCC_gleason_plot.png"),
       p_gleason,
       width = 10,
       height = 5)

# ======================================
# SURVIAL ANALYSIS WITH KM PLOTS
# ======================================

# Build DFS survival object
library(survival)
library(survminer)

table(meta$DFS_STATUS, useNA = "ifany")

dfs_df <- meta %>%
  mutate(sample_id = sample_id) %>%
  select(sample_id, DFS_MONTHS, DFS_STATUS) %>%
  mutate(
    time = as.numeric(DFS_MONTHS),
    event = case_when(
      DFS_STATUS %in% c("1:Recurred") ~ 1,
      DFS_STATUS %in% c("0:DiseaseFree") ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(time), !is.na(event))

# Merge with expression 

expr_df <- as.data.frame(t(expr))
expr_df$sample_id <- rownames(expr_df)

km_df <- expr_df %>%
  left_join(dfs_df, by = "sample_id") %>%
  filter(!is.na(time), !is.na(event))

# KM + Log Rank per gene

genes <- rownames(expr)

logrank_tbl <- lapply(genes, function(gene) {
  
  tmp <- km_df %>%
    select(time, event, all_of(gene)) %>%
    rename(expression = all_of(gene)) %>%
    mutate(expression = as.numeric(expression)) %>%
    filter(!is.na(time), !is.na(event), !is.na(expression))
  
  # ❗ SAFETY CHECK 1: enough samples
  if (nrow(tmp) < 10) {
    return(data.frame(Gene = gene, p.value = NA))
  }
  
  # median split
  med <- median(tmp$expression, na.rm = TRUE)
  
  tmp <- tmp %>%
    mutate(group = ifelse(expression >= med, "High", "Low"))
  
  # ❗ SAFETY CHECK 2: must have both groups
  if (length(unique(tmp$group)) < 2) {
    return(data.frame(Gene = gene, p.value = NA))
  }
  
  sd <- survdiff(Surv(time, event) ~ group, data = tmp)
  p <- 1 - pchisq(sd$chisq, df = length(sd$n) - 1)
  
  data.frame(Gene = gene, p.value = p)
})

# KM plots
dir.create(file.path(outdir, "DFS_KM_plots"), showWarnings = FALSE)

plot_list <- list()

for (gene in genes) {
  
  tmp <- km_df %>%
    select(time, event, all_of(gene)) %>%
    rename(expression = all_of(gene)) %>%
    mutate(
      group = ifelse(expression >= median(expression, na.rm = TRUE),
                     "High", "Low")
    ) %>%
    filter(!is.na(time), !is.na(event), !is.na(expression))
  
  tmp$group <- factor(tmp$group, levels = c("Low", "High"))
  
  fit <- survfit(Surv(time, event) ~ group, data = tmp)
  
  p <- ggsurvplot(
    fit,
    data = tmp,
    pval = TRUE,
    risk.table = TRUE,
    title = paste0(gene, " (DFS)")
  )
  
  plot_list[[gene]] <- p$plot
  
  ggsave(
    file.path(outdir, "DFS_KM_plots", paste0(gene, "_DFS_KM.png")),
    p$plot,
    width = 6,
    height = 5
  )
}

# Combined KM Plots 
library(ggpubr)

combined <- ggarrange(
  plotlist = plot_list,
  ncol = 3,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom"
)

ggsave(
  file.path(outdir, "MSKCC_KM_DFS_combined.png"),
  combined,
  width = 14,
  height = 5
)

# =========================================================
# HEATMAPS: CO-EXPRESSION STYLE
# =========================================================

library(pheatmap)

heatdir <- file.path(outdir, "heatmaps")
dir.create(heatdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# 1. Ensure full expression matrix is numeric
# ------------------------------------------------------------------
expr_full <- as.matrix(expr_mat)
mode(expr_full) <- "numeric"

# ------------------------------------------------------------------
# 2. Sample annotation
# ------------------------------------------------------------------
meta_heat <- meta %>%
  select(sample_id, group) %>%
  filter(!is.na(sample_id), !is.na(group)) %>%
  group_by(sample_id) %>%
  summarise(group = first(group), .groups = "drop") %>%
  as.data.frame()

meta_heat <- meta_heat[match(colnames(expr_full), meta_heat$sample_id), ]
meta_heat <- meta_heat[!is.na(meta_heat$sample_id), ]
meta_heat <- meta_heat[!duplicated(meta_heat$sample_id), ]
rownames(meta_heat) <- meta_heat$sample_id
meta_heat <- meta_heat[, "group", drop = FALSE]
meta_heat$group <- factor(meta_heat$group,
                          levels = c("Primary Tumor", "Metastatic Tumor"))

# ------------------------------------------------------------------
# 3. Target gene ID -> symbol mapping
# ------------------------------------------------------------------
target_entrez <- c(
  FCGR3A = "2214",
  FOLH1  = "2346",
  PTGS2  = "5743"
)

target_entrez <- target_entrez[target_entrez %in% rownames(expr_full)]
stopifnot(length(target_entrez) > 0)

# Build a symbol-named copy of the expression matrix
expr_sym <- expr_full
rownames(expr_sym) <- rownames(expr_full)

for (sym in names(target_entrez)) {
  entrez_id <- target_entrez[[sym]]
  row_idx <- which(rownames(expr_sym) == entrez_id)
  if (length(row_idx) == 1) {
    rownames(expr_sym)[row_idx] <- sym
  }
}

# Remove duplicate rownames if any were created
expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), , drop = FALSE]

# ------------------------------------------------------------------
# 4. Co-expression heatmap function
# ------------------------------------------------------------------
make_coexpression_heatmap <- function(target_gene, expr_mat, top_n = 15) {
  if (!target_gene %in% rownames(expr_mat)) {
    warning(paste("Target gene not found:", target_gene))
    return(NULL)
  }
  
  target_vec <- as.numeric(expr_mat[target_gene, ])
  
  cors <- apply(expr_mat, 1, function(x) {
    cor(as.numeric(x), target_vec,
        method = "pearson",
        use = "pairwise.complete.obs")
  })
  
  cors <- cors[!is.na(cors)]
  cors <- sort(cors, decreasing = TRUE)
  
  pos_genes <- names(head(cors[cors > 0], top_n))
  neg_genes <- names(head(sort(cors[cors < 0]), top_n))
  
  genes_keep <- unique(c(target_gene, pos_genes, neg_genes))
  genes_keep <- genes_keep[genes_keep %in% rownames(expr_mat)]
  
  sub_mat <- expr_mat[genes_keep, , drop = FALSE]
  sub_mat <- as.matrix(sub_mat)
  mode(sub_mat) <- "numeric"
  
  cor_mat <- cor(t(sub_mat), use = "pairwise.complete.obs")
  
  pheatmap(
    cor_mat,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = paste0(target_gene, " co-expression"),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    border_color = NA
  )
}

# ------------------------------------------------------------------
# 5. Generate one heatmap per target gene
# ------------------------------------------------------------------
genes <- intersect(names(target_entrez), rownames(expr_sym))

for (g in genes) {
  png(
    file.path(heatdir, paste0(g, "_heatmap.png")),
    width = 1800, height = 1800, res = 250
  )
  make_coexpression_heatmap(g, expr_sym, top_n = 15)
  dev.off()
}
      
      
      