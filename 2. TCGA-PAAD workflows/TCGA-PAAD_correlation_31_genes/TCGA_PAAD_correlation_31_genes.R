library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(writexl)

outdir <- "output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load data
data <- readRDS("TCGA_panCancer_FOLH1.rds")

count_mat <- assay(data)
meta <- as.data.frame(colData(data))

# Extract PAAD samples
meta_paad <- meta[meta$project == "TCGA-PAAD", ]

# Get PAAD count matrix
paad_counts <- count_mat[, rownames(meta_paad)]

# Get gene annotation
gene_annot <- as.data.frame(rowData(data))
gene_annot_paad <- gene_annot[rownames(paad_counts), ]

dds <- DESeqDataSetFromMatrix(
  countData = paad_counts,
  colData = meta_paad,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)

norm_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# Save annotation
write.csv(gene_annot_paad,
          file = paste0(outdir, "/PAAD_gene_annotation.csv"),
          row.names = FALSE)

# Make a lookup table: Ensembl IDs (with version) -> gene symbols
gene_map <- data.frame(
  ensembl_id = rownames(gene_annot_paad),
  gene_name = gene_annot_paad$gene_name,
  stringsAsFactors = FALSE
)

gene_map <- gene_map[!is.na(gene_map$gene_name) & gene_map$gene_name != "", ]
gene_map <- gene_map[!duplicated(gene_map$gene_name), ]

# Gene list
gene_list <- c(
  "RIT1", "KRAS", "NF1", "PIK3CA", "VEGFA", "AKT3", "FLT1", "MAP2K1", "EGFR", "BRAF",
  "TSC1", "PIK3R1", "ROS1", "MET", "NRAS", "CCNA1", "CCND1", "PTEN", "KDR", "MMP9",
  "EGF", "AKT2", "ALK", "RET", "TP53", "ERBB2", "TSC2", "AKT1", "HRAS", "STK11"
)

target <- "PIP5K1A"
all_genes <- unique(c(target, gene_list))

# Match requested symbols
keep_map <- gene_map[gene_map$gene_name %in% all_genes, ]

# Subset VST matrix by Ensembl IDs
vsd_sub <- vsd_mat[keep_map$ensembl_id, , drop = FALSE]

# Replace row names with gene symbols
rownames(vsd_sub) <- keep_map$gene_name

# Check which genes were found/missing
found_genes <- rownames(vsd_sub)
missing_genes <- setdiff(all_genes, found_genes)

write.csv(data.frame(found_genes = found_genes),
          file = paste0(outdir, "/found_genes.csv"),
          row.names = FALSE)

write.csv(data.frame(missing_genes = missing_genes),
          file = paste0(outdir, "/missing_genes.csv"),
          row.names = FALSE)

# Correlation against PIP5K1A using VST values
target_vec <- as.numeric(vsd_sub[target, ])

corr_res <- do.call(rbind, lapply(gene_list, function(g) {
  if (!g %in% rownames(vsd_sub)) {
    return(data.frame(gene = g, rho = NA, p.value = NA))
  }
  x <- as.numeric(vsd_sub[g, ])
  ct <- cor.test(x, target_vec, method = "spearman", use = "complete.obs")
  data.frame(gene = g, rho = unname(ct$estimate), p.value = ct$p.value)
}))

corr_res <- corr_res[order(corr_res$p.value), ]

write.csv(corr_res,
          file = paste0(outdir, "/PIP5K1A_31gene_correlations_VST.csv"),
          row.names = FALSE)
corr_res

plot_df <- corr_res %>%
  dplyr::filter(!is.na(rho)) %>%
  dplyr::arrange(rho) %>%
  dplyr::mutate(gene = factor(gene, levels = gene))

p <- ggplot(plot_df, aes(x = rho, y = gene)) +
  geom_col(fill = "#2C7FB8") +
  geom_vline(xintercept = 0, color = "steelblue", linewidth = 0.6) +
  labs(
    title = "PIP5K1A correlations with available genes",
    x = "Spearman rho with PIP5K1A",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 11)
  )

ggsave(
  filename = file.path(outdir, "PIP5K1A_corr_plot.png"),
  plot = p,
  width = 8,
  height = 7,
  dpi = 300
)

p

# individual plots

plot_dir <- file.path(outdir, "pairwise_corr_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

x_gene <- "PIP5K1A"
expr_mat <- vsd_sub

for (g in gene_list) {
  if (!all(c(x_gene, g) %in% rownames(expr_mat))) next
  
  df_plot <- data.frame(
    x = as.numeric(expr_mat[x_gene, ]),
    y = as.numeric(expr_mat[g, ])
  )
  
  ct <- cor.test(df_plot$x, df_plot$y, method = "spearman", use = "complete.obs")
  rho_val <- unname(ct$estimate)
  p_val <- ct$p.value
  
  p <- ggplot(df_plot, aes(x = x, y = y)) +
    geom_point(color = "#2C7FB8", alpha = 0.85, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.9) +
    labs(
      title = paste0(x_gene, " and ", g, " correlation"),
      x = paste0(x_gene, " (VST)"),
      y = paste0(g, " (VST)")
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("R = ", sprintf("%.4f", rho_val), "\np = ", signif(p_val, 3)),
      hjust = 1.1, vjust = 1.2, size = 4
    ) +
    theme_classic(base_size = 14)
  
  ggsave(
    filename = file.path(plot_dir, paste0(x_gene, "_vs_", g, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

# Combined plot

plot_list <- list()

for (g in gene_list) {
  if (!all(c(x_gene, g) %in% rownames(expr_mat))) next
  
  df_plot <- data.frame(
    x = as.numeric(expr_mat[x_gene, ]),
    y = as.numeric(expr_mat[g, ])
  )
  
  ct <- cor.test(df_plot$x, df_plot$y, method = "spearman", use = "complete.obs")
  rho_val <- unname(ct$estimate)
  p_val <- ct$p.value
  
  p <- ggplot(df_plot, aes(x = x, y = y)) +
    geom_point(color = "#2C7FB8", alpha = 0.75, size = 0.9) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 0.5) +
    labs(
      title = g,
      x = NULL,
      y = NULL
    ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = paste0("R = ", sprintf("%.2f", rho_val), "\np = ", signif(p_val, 2)),
      hjust = 1.05, vjust = 1.15, size = 2.8
    ) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 9),
      axis.text = element_text(size = 7)
    )
  
  plot_list[[g]] <- p
}

n_plots <- length(plot_list)
ncol <- 4
nrow <- ceiling(n_plots / ncol)

combined <- ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow)

ggsave(
  filename = file.path(outdir, "PIP5K1A_all_gene_correlations_grid.png"),
  plot = combined,
  width = 16,
  height = 4 * nrow,
  dpi = 300
)
