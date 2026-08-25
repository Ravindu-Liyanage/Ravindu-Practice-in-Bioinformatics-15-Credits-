library(SummarizedExperiment)
library(MultiAssayExperiment)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(limma)
library(survival)
library(survminer)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

outdir <- "output_MSKCC2010_mets"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 1) LOAD DATA
# =========================================================
data <- readRDS("MSKCC2010.rds")

rna <- experiments(data)[["mrna_agilent_microarray"]]

expr_mat <- assay(rna)
expr_mat <- as.matrix(expr_mat)
rownames(expr_mat) <- as.character(rownames(expr_mat))

meta <- as.data.frame(colData(data))
meta$sample_id <- rownames(meta)

# =========================================================
# 2) MATCH SAMPLES BETWEEN EXPRESSION AND METADATA
# =========================================================
common_samples <- intersect(colnames(expr_mat), meta$sample_id)
expr_mat <- expr_mat[, common_samples, drop = FALSE]
meta <- meta[match(common_samples, meta$sample_id), , drop = FALSE]

stopifnot(all(colnames(expr_mat) == meta$sample_id))

# =========================================================
# 3) FIND THE SAMPLE-TYPE COLUMN
# =========================================================
meta_cols <- colnames(meta)
group_col <- NA_character_

candidate_cols <- c(
  "SAMPLE_TYPE", "sample_type", "sampleType",
  "sample_group", "primary_metastatic", "condition",
  "group", "status"
)

for (cc in candidate_cols) {
  if (cc %in% meta_cols) {
    group_col <- cc
    break
  }
}

if (is.na(group_col)) {
  stop("No obvious sample-type column found. Inspect colnames(meta) and set group_col manually.")
}

print(table(meta[[group_col]], useNA = "ifany"))

# =========================================================
# 4) KEEP ONLY METASTATIC SAMPLES
# =========================================================
metastatic_levels <- c("Metastatic", "metastasis", "Metastasis", "METS", "METASTATIC")

meta_mets <- meta[as.character(meta[[group_col]]) %in% metastatic_levels, , drop = FALSE]
expr_mets <- expr_mat[, meta_mets$sample_id, drop = FALSE]

if (nrow(meta_mets) < 2) {
  stop("Too few metastatic samples found after filtering.")
}

# =========================================================
# 5) CHECK EXPRESSION SCALE
# =========================================================
expr_range <- range(expr_mets, na.rm = TRUE)
print(expr_range)

# If your values look raw/non-log (very large), uncomment:
# expr_mets <- log2(expr_mets + 1)

# =========================================================
# 6) CLEAN MATRIX
# =========================================================
expr_mets <- expr_mets[rowSums(is.na(expr_mets)) < ncol(expr_mets), , drop = FALSE]
expr_mets <- expr_mets[!duplicated(rownames(expr_mets)), , drop = FALSE]

# Optional: remove very low-variance genes to stabilize analysis
gene_var <- apply(expr_mets, 1, var, na.rm = TRUE)
expr_mets <- expr_mets[gene_var > 0, , drop = FALSE]

# =========================================================
# 7) DEFINE IMMUNE GENE SETS IN SYMBOLS
# =========================================================
immune_symbol_sets <- list(
  T_cells = c("CD3D","CD3E","CD2","TRAC","LCK","IL7R"),
  CD8_T = c("CD8A","CD8B","NKG7","GZMB","PRF1"),
  Cytotoxic = c("NKG7","GZMB","PRF1","GNLY","KLRD1"),
  B_cells = c("MS4A1","CD79A","CD79B","CD19","CD74"),
  NK_cells = c("NKG7","GNLY","KLRD1","FCGR3A"),
  Monocytes_Macrophages = c("LYZ","FCN1","CD68","CSF1R","TYROBP"),
  Treg = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  Checkpoint = c("PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2"),
  Antigen_Presentation = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2")
)

# =========================================================
# 8) CONVERT IMMUNE SIGNATURES: SYMBOL -> ENTREZID
# =========================================================
convert_symbols_to_entrez <- function(symbols) {
  map_df <- suppressMessages(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  unique(as.character(map_df$ENTREZID))
}

immune_entrez_sets <- lapply(immune_symbol_sets, convert_symbols_to_entrez)
immune_entrez_sets <- lapply(immune_entrez_sets, function(gs) intersect(gs, rownames(expr_mets)))
immune_entrez_sets <- immune_entrez_sets[lengths(immune_entrez_sets) >= 2]

if (length(immune_entrez_sets) == 0) {
  stop("No immune gene sets remained after SYMBOL->ENTREZ conversion and intersection with expr_mets.")
}

print(sapply(immune_entrez_sets, length))

# =========================================================
# 9) RUN GSVA USING ENTREZ-BASED GENE SETS
# =========================================================
gsva_param <- gsvaParam(as.matrix(expr_mets), immune_entrez_sets)
gsva_scores <- gsva(gsva_param)
gsva_scores <- as.data.frame(gsva_scores)

write.csv(gsva_scores, file.path(outdir, "MSKCC2010_GSVA_immune_scores_entrez.csv"))

# =========================================================
# 10) CREATE A GLOBAL IMMUNE SCORE
# =========================================================
immune_score_vec <- colMeans(gsva_scores, na.rm = TRUE)
meta_mets$immune_score <- as.numeric(immune_score_vec[meta_mets$sample_id])

meta_mets$immune_group <- ifelse(
  meta_mets$immune_score >= median(meta_mets$immune_score, na.rm = TRUE),
  "High",
  "Low"
)
meta_mets$immune_group <- factor(meta_mets$immune_group, levels = c("Low", "High"))

write.csv(meta_mets, file.path(outdir, "MSKCC2010_metadata_with_immune_score.csv"), row.names = FALSE)

# =========================================================
# 11) PLOT GSVA SCORES BY IMMUNE GROUP
# =========================================================
score_df <- meta_mets %>%
  dplyr::select(sample_id, immune_group) %>%
  left_join(
    as.data.frame(t(gsva_scores)) %>% tibble::rownames_to_column("sample_id"),
    by = "sample_id"
  )

sig_names <- rownames(gsva_scores)

score_long <- tidyr::pivot_longer(
  score_df,
  cols = all_of(sig_names),
  names_to = "signature",
  values_to = "score"
)

pvals <- score_long %>%
  group_by(signature) %>%
  summarise(
    p.value = wilcox.test(score ~ immune_group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("p = ", signif(p.value, 3)))

label_df <- score_long %>%
  group_by(signature) %>%
  summarise(
    x = 1,
    y = max(score, na.rm = TRUE) * 0.85,
    .groups = "drop"
  ) %>%
  left_join(pvals, by = "signature")

p1 <- ggplot(score_long, aes(x = immune_group, y = score, fill = immune_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.7) +
  facet_wrap(~signature, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Immune score group", y = "GSVA score") +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3
  )

ggsave(file.path(outdir, "MSKCC2010_GSVA_scores_by_group.png"), p1, width = 12, height = 8, dpi = 300)

png(file.path(outdir, "MSKCC2010_GSVA_heatmap.png"), width = 9, height = 6, units = "in", res = 300)
pheatmap(as.matrix(gsva_scores), scale = "row", show_colnames = FALSE,
         main = "GSVA immune signatures (Entrez-based)")
dev.off()

# =========================================================
# 12) DIFFERENTIAL EXPRESSION: HIGH vs LOW IMMUNE SCORE
# =========================================================
design <- model.matrix(~ immune_group, data = meta_mets)
fit <- lmFit(expr_mets, design)
fit <- eBayes(fit)

deg <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
deg$ENTREZID <- rownames(deg)

# Map DEG Entrez IDs to gene symbols for readability
deg_symbol_map <- suppressMessages(
  bitr(
    deg$ENTREZID,
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
  )
)

deg <- deg %>%
  left_join(deg_symbol_map, by = "ENTREZID")

write.csv(deg, file.path(outdir, "MSKCC2010_DEG_High_vs_Low_ImmuneScore.csv"), row.names = FALSE)

sig_deg <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

write.csv(sig_deg, file.path(outdir, "MSKCC2010_DEG_significant.csv"), row.names = FALSE)

# =========================================================
# 13) ENRICHMENT ANALYSIS
# =========================================================
if (nrow(sig_deg) > 0) {
  entrez_ids <- unique(sig_deg$ENTREZID)
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  ego <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  ekegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pAdjustMethod = "BH"
  )
  
  write.csv(as.data.frame(ego), file.path(outdir, "MSKCC2010_GO_enrichment.csv"), row.names = FALSE)
  write.csv(as.data.frame(ekegg), file.path(outdir, "MSKCC2010_KEGG_enrichment.csv"), row.names = FALSE)
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    p_go <- dotplot(ego, showCategory = 10, title = "Top GO Biological Processes")
    ggsave(file.path(outdir, "MSKCC2010_GO_dotplot.png"), p_go, width = 8, height = 6, dpi = 300)
  }
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    p_kegg <- dotplot(ekegg, showCategory = 10, title = "Top KEGG Pathways")
    ggsave(file.path(outdir, "MSKCC2010_KEGG_dotplot.png"), p_kegg, width = 8, height = 6, dpi = 300)
  }
}

# =========================================================
# 14) CORRELATE FOLH1 WITH IMMUNE SCORE
# =========================================================
folh1_map <- suppressMessages(
  bitr(
    "FOLH1",
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

folh1_entrez <- unique(as.character(folh1_map$ENTREZID))[1]

if (!is.na(folh1_entrez) && folh1_entrez %in% rownames(expr_mets)) {
  cor_test <- cor.test(
    as.numeric(expr_mets[folh1_entrez, meta_mets$sample_id]),
    meta_mets$immune_score,
    method = "spearman"
  )
  
  rho_lab <- sprintf("Spearman rho = %.3f", unname(cor_test$estimate))
  p_lab <- sprintf("p-value = %.3g", cor_test$p.value)
  
  plot_df <- data.frame(
    gene_expr = as.numeric(expr_mets[folh1_entrez, meta_mets$sample_id]),
    immune_score = meta_mets$immune_score
  )
  
  p2 <- ggplot(plot_df, aes(gene_expr, immune_score)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    labs(
      title = "FOLH1 Expression vs Immune Score",
      x = paste0("FOLH1 expression (Entrez ", folh1_entrez, ")"),
      y = "GSVA immune score"
    ) +
    annotate(
      "text",
      x = -Inf, y = Inf,
      label = paste(rho_lab, p_lab, sep = "\n"),
      hjust = -0.1, vjust = 1.1,
      size = 4
    )
  
  ggsave(
    file.path(outdir, "FOLH1_vs_MSKCC2010_immune.png"),
    p2, width = 5, height = 5, dpi = 300
  )
}

# =========================================================
# 15) SURVIVAL ANALYSIS
# =========================================================
if (all(c("DFS_MONTHS", "DFS_STATUS") %in% colnames(meta_mets))) {
  meta_mets$time <- suppressWarnings(as.numeric(as.character(meta_mets$DFS_MONTHS)))
  meta_mets$event <- ifelse(as.character(meta_mets$DFS_STATUS) == "1:Recurred", 1, 0)
  
  surv_df <- meta_mets[!is.na(meta_mets$time) & !is.na(meta_mets$event) & !is.na(meta_mets$immune_group), , drop = FALSE]
  
  if (nrow(surv_df) >= 6 && length(unique(surv_df$immune_group)) == 2) {
    fit_surv <- survfit(Surv(time, event) ~ immune_group, data = surv_df)
    g <- ggsurvplot(
      fit_surv,
      data = surv_df,
      pval = TRUE,
      risk.table = TRUE,
      title = "Kaplan-Meier Survival Curves by Immune Score Group",
      legend.title = "Immune group",
      legend.labs = c("Low", "High")
    )
    ggsave(file.path(outdir, "aMSKCC2010_KM_immune_group.png"), plot = g$plot, width = 7, height = 6, dpi = 300)
  }
}
# =========================================================
# 16) SAVE FINAL TABLES
# =========================================================
write.csv(meta_mets, file.path(outdir, "MSKCC2010_final_sample_table.csv"), row.names = FALSE)

immune_signature_sizes <- data.frame(
  signature = names(immune_entrez_sets),
  n_genes_used = sapply(immune_entrez_sets, length),
  stringsAsFactors = FALSE
)

write.csv(immune_signature_sizes, file.path(outdir, "MSKCC2010_immune_signature_sizes.csv"), row.names = FALSE)