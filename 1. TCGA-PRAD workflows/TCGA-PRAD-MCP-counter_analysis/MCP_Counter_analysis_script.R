# =========================================================
# MCP-counter analysis for TCGA-PRAD primary tumors
# Input required: gene-symbol expression matrix with genes in rows
# Uses the existing logcpm_symbol matrix and expr_group from your pipeline
# =========================================================

# If needed:
# install.packages("devtools")
# devtools::install_github("ebecht/MCPcounter", subdir = "Source")

library(MCPcounter)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# -----------------------------
# 1) Prepare MCP-counter input
# -----------------------------
# logcpm_symbol must already exist from your GSVA workflow
# It should have gene symbols as rownames and samples as columns
mcp_input <- as.matrix(logcpm_symbol)

# Keep only finite rows
mcp_input <- mcp_input[rowSums(is.na(mcp_input)) < ncol(mcp_input), , drop = FALSE]

# -----------------------------
# 2) Run MCP-counter
# -----------------------------
mcp_scores <- MCPcounter.estimate(
  expression = mcp_input,
  featuresType = "HUGO_symbols"
)

# MCPcounter returns cell populations in rows and samples in columns
mcp_scores_df <- as.data.frame(t(mcp_scores))
mcp_scores_df <- mcp_scores_df %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(
    meta_prad %>%
      tibble::rownames_to_column("sample_id") %>%
      dplyr::select(sample_id, expr_group, target_gene, target_expr, target_ensembl),
    by = "sample_id"
  )

write.csv(
  mcp_scores_df,
  file.path(outdir, paste0("MCPcounter_scores_", candidate, ".csv")),
  row.names = FALSE
)

# -----------------------------
# 3) Long format for plots
# -----------------------------
mcp_long <- mcp_scores_df %>%
  pivot_longer(
    cols = -c(sample_id, expr_group, target_gene, target_expr, target_ensembl),
    names_to = "cell_population",
    values_to = "score"
  )

# -----------------------------
# 4) Group comparison plots
# -----------------------------
p_mcp_box <- ggplot(mcp_long, aes(x = expr_group, y = score, fill = expr_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.7) +
  facet_wrap(~ cell_population, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = paste0(candidate, " expression group"),
    y = "MCP-counter score",
    title = paste0("Tumor microenvironment by ", candidate, " High vs Low")
  ) +
  scale_fill_manual(values = c("Low" = "#4575B4", "High" = "#D73027"))

ggsave(
  file.path(outdir, paste0("MCPcounter_boxplots_", candidate, ".png")),
  p_mcp_box, width = 14, height = 10, dpi = 300
)

# -----------------------------
# 5) Statistical testing
# -----------------------------
mcp_stats <- mcp_long %>%
  group_by(cell_population) %>%
  summarise(
    p.value = wilcox.test(score ~ expr_group)$p.value,
    median_low = median(score[expr_group == "Low"], na.rm = TRUE),
    median_high = median(score[expr_group == "High"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(adj.P.Val = p.adjust(p.value, method = "BH"))

write.csv(
  mcp_stats,
  file.path(outdir, paste0("MCPcounter_stats_", candidate, ".csv")),
  row.names = FALSE
)

# -----------------------------
# 6) Heatmap of MCP-counter scores
# -----------------------------
mcp_mat <- as.matrix(mcp_scores)
mcp_mat <- mcp_mat[rownames(mcp_mat) %in% unique(mcp_long$cell_population), , drop = FALSE]

ann_col <- data.frame(expr_group = meta_prad$expr_group)
rownames(ann_col) <- rownames(meta_prad)
ann_col <- ann_col[colnames(mcp_mat), , drop = FALSE]

png(
  file.path(outdir, paste0("MCPcounter_heatmap_", candidate, ".png")),
  width = 10, height = 7, units = "in", res = 300
)

pheatmap::pheatmap(
  mcp_mat,
  scale = "row",
  annotation_col = ann_col,
  main = paste0("MCP-counter cell populations: ", candidate, " High vs Low")
)

dev.off()

# -----------------------------
# 7) Optional correlation with FOLH1 expression
# -----------------------------
mcp_cor <- mcp_long %>%
  group_by(cell_population) %>%
  summarise(
    cor_spearman = cor(score, target_expr, method = "spearman", use = "complete.obs"),
    p.value = cor.test(score, target_expr, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(adj.P.Val = p.adjust(p.value, method = "BH"))

write.csv(
  mcp_cor,
  file.path(outdir, paste0("MCPcounter_correlation_", candidate, ".csv")),
  row.names = FALSE
)