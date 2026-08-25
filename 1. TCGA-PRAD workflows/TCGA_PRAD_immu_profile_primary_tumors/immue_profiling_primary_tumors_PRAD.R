# ========================================================
# Ravindu: TCGA-PRAD Immune profiling of Primary Tumors
# ========================================================
 
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(edgeR)
library(limma)
library(survival)
library(survminer)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)

outdir <- "output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --------------
# 1) Load data
# --------------
data <- readRDS("TCGA_panCancer_FOLH1.rds")

count_mat <- assay(data)
meta <- as.data.frame(colData(data))
gene_annot <- as.data.frame(rowData(data))

# Subset PRAD only
meta_prad <- meta[meta$project == "TCGA-PRAD", , drop = FALSE]
prad_counts <- count_mat[, rownames(meta_prad), drop = FALSE]

# keep only primary tumor samples
tumor_keep <- meta_prad$sample_type == "Primary Tumor"
meta_prad <- meta_prad[tumor_keep, , drop = FALSE]
prad_counts <- prad_counts[, rownames(meta_prad), drop = FALSE]

# Clean up gene annotation
if ("gene_name" %in% colnames(gene_annot)) {
  rownames(prad_counts) <- make.unique(as.character(gene_annot$gene_name))
} else if ("external_gene_name" %in% colnames(gene_annot)) {
  rownames(prad_counts) <- make.unique(as.character(gene_annot$external_gene_name))
} else if ("symbol" %in% colnames(gene_annot)) {
  rownames(prad_counts) <- make.unique(as.character(gene_annot$symbol))
}

prad_counts <- prad_counts[rowSums(prad_counts, na.rm = TRUE) > 0, ]

# --------------------------
# 2) Normalize expression
# --------------------------
dge <- DGEList(counts = prad_counts)
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

write.csv(as.data.frame(logcpm), file.path(outdir, "PRAD_logCPM.csv"))

# -----------------------------------
# 3) Define immune gene sets and run GSVA
# -----------------------------------
gene_sets <- list(
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

gene_sets <- lapply(gene_sets, function(gs) intersect(gs, rownames(logcpm)))
gene_sets <- gene_sets[lengths(gene_sets) >= 2]


library(GSVA)

gsva_param <- gsvaParam(
  as.matrix(logcpm),
  gene_sets
)
gsva_scores <- gsva(gsva_param)
gsva_scores <- as.data.frame(gsva_scores)

write.csv(gsva_scores, file.path(outdir, "GSVA_immune_scores.csv"))

# -------------------------------------------
# 4) Create an immune-high / immune-low score
# --------------------------------------------

immune_pan_score <- colMeans(gsva_scores, na.rm = TRUE)
meta_prad$immune_score <- as.numeric(immune_pan_score[rownames(meta_prad)])

meta_prad$immune_group <- ifelse(meta_prad$immune_score >= median(meta_prad$immune_score, na.rm = TRUE),
                                 "High", "Low")
meta_prad$immune_group <- factor(meta_prad$immune_group, levels = c("Low", "High"))
meta_prad_export <- meta_prad %>%
  dplyr::mutate(dplyr::across(where(is.list), ~ vapply(.x, toString, character(1))))

write.csv(meta_prad_export, file.path(outdir, "PRAD_metadata_with_GSVA.csv"), row.names = FALSE)

# --------------------
# 5) Plot GSVA scores
# --------------------

library(dplyr)
library(tidyr)
library(ggplot2)

score_df <- meta_prad %>%
  dplyr::select(immune_group) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(
    as.data.frame(t(gsva_scores)) %>%
      tibble::rownames_to_column("sample_id"),
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

ggsave(file.path(outdir, "GSVA_scores_by_group.png"), p1, width = 12, height = 8, dpi = 300)

png(file.path(outdir, "GSVA_heatmap.png"), width = 9, height = 6, units = "in", res = 300)
pheatmap(as.matrix(gsva_scores), scale = "row", show_colnames = FALSE,
         main = "GSVA immune signatures")
dev.off()

# ----------------------------------------------------
# 6) Differential expression: High vs Low immune score
# ----------------------------------------------------

design <- model.matrix(~ immune_group, data = meta_prad)
v <- voom(prad_counts[, rownames(meta_prad)], design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

deg <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
deg$gene <- rownames(deg)
write.csv(deg, file.path(outdir, "DEG_High_vs_Low_GSVAimmune.csv"), row.names = FALSE)

sig_deg <- deg %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(sig_deg, file.path(outdir, "DEG_significant.csv"), row.names = FALSE)

# --------------------------------
# 7) Enrichment analysis for DEGs
# --------------------------------
if (nrow(sig_deg) > 0) {
  entrez_df <- clusterProfiler::bitr(
    sig_deg$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  entrez_ids <- unique(entrez_df$ENTREZID)
  
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
  
  write.csv(as.data.frame(ego), file.path(outdir, "GO_enrichment.csv"), row.names = FALSE)
  write.csv(as.data.frame(ekegg), file.path(outdir, "KEGG_enrichment.csv"), row.names = FALSE)
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    p_go <- dotplot(ego, showCategory = 10, title = "Top GO Biological Processes")
    ggsave(file.path(outdir, "GO_dotplot.png"), p_go, width = 8, height = 6, dpi = 300)
  }
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    p_kegg <- dotplot(ekegg, showCategory = 10, title = "Top KEGG Pathways")
    ggsave(file.path(outdir, "KEGG_dotplot.png"), p_kegg, width = 8, height = 6, dpi = 300)
  }
}

# -----------------------------------------------
# 8) Correlate a candidate gene with immune score
# -----------------------------------------------

candidate <- "FOLH1"

if (candidate %in% rownames(logcpm)) {
  cor_test <- cor.test(
    as.numeric(logcpm[candidate, rownames(meta_prad)]),
    meta_prad$immune_score,
    method = "spearman"
  )
  
  rho_lab <- sprintf("Spearman rho = %.3f", unname(cor_test$estimate))
  p_lab <- sprintf("p-value = %.3g", cor_test$p.value)
  
  plot_df <- data.frame(
    gene_expr = as.numeric(logcpm[candidate, rownames(meta_prad)]),
    immune_score = meta_prad$immune_score
  )
  
  p2 <- ggplot(plot_df, aes(gene_expr, immune_score)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    theme_classic() +
    labs(
      title = "FOLH1 Expression vs Immune Score",
      x = paste0(candidate, " logCPM"),
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
    file.path(outdir, paste0(candidate, "_vs_GSVAimmune.png")),
    p2, width = 5, height = 5, dpi = 300
  )
}

# ------------------------
# 12) Survival analysis
# ------------------------

time_col <- intersect(c("days_to_death", "days_to_last_follow_up"), colnames(meta_prad))
event_col <- intersect(c("vital_status"), colnames(meta_prad))

if (length(time_col) > 0 && length(event_col) > 0) {
  meta_prad$time <- ifelse(!is.na(meta_prad$days_to_death) & meta_prad$days_to_death > 0,
                           as.numeric(meta_prad$days_to_death),
                           as.numeric(meta_prad$days_to_last_follow_up))
  meta_prad$event <- ifelse(meta_prad$vital_status == "Dead", 1, 0)
  
  fit_surv <- survfit(Surv(time, event) ~ immune_group, data = meta_prad)
  g <- ggsurvplot(
    fit_surv,
    data = meta_prad,
    pval = TRUE,
    risk.table = TRUE,
    title = "Kaplan–Meier Survival Curves by Immune Score Group",
    legend.title = "Immune group",
    legend.labs = c("Low", "High")
  )
  
  g$plot <- g$plot + theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 10))
  )
  
  ggsave(
    file.path(outdir, "KM_GSVAimmune_group.png"),
    plot = g$plot,
    width = 7,
    height = 6,
    dpi = 300
  )
}

# -----------------------------
# 13) Save final sample table
# -----------------------------
write.csv(
  meta_prad %>% dplyr::select(where(~ !is.list(.x))),
  file.path(outdir, "PRAD_final_sample_table.csv"),
  row.names = FALSE
)
