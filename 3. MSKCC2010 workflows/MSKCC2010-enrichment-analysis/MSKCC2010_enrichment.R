# =========================================================
# MSKCC2010 Differential Expression + Reactome ORA Pipeline
# Using limma for microarray data (instead of edgeR for RNA-seq)
# =========================================================

# Load packages
library(limma)
library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(biomaRt)
library(ReactomePA)
library(enrichplot)
library(org.Hs.eg.db)
library(tibble)
library(stringr)
library(clusterProfiler)

# =========================================================
# OUTPUT DIRECTORIES
# =========================================================

outdir <- "output_MSKCC2010_e"
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
# METADATA (ONLY use rna-level metadata)
# =========================================================

meta <- as.data.frame(colData(data))   # NOT colData(rna)
meta$sample_id <- rownames(meta)

# Check what columns are available
cat("Metadata columns:\n")
print(colnames(meta))

cat("\nSAMPLE_TYPE values:\n")
print(unique(meta$SAMPLE_TYPE))

# =========================================================
# 1) Extract metadata matching expression matrix columns
# =========================================================

# Keep only samples that have expression data
keep_samples <- rownames(meta) %in% colnames(expr_mat)
meta <- meta[keep_samples, ]

# Verify matching
cat("\nExpression matrix columns:", length(colnames(expr_mat)), "\n")
cat("Metadata rows:", length(rownames(meta)), "\n")

# =========================================================
# 2) Build grouping variable from SAMPLE_TYPE
# =========================================================

# Check the actual values in SAMPLE_TYPE
table(meta$SAMPLE_TYPE)

# Create group factor - adjust levels based on actual values
# Assuming values are "Primary" and "Metastasis"
meta$group <- factor(
  meta$SAMPLE_TYPE,
  levels = c("Primary", "Metastasis"),  # Change if values differ
  labels = c("Primary", "Metastasis")
)

table(meta$group)

# Keep only Primary and Metastasis samples
keep <- meta$group %in% c("Primary", "Metastasis")
meta <- meta[keep, , drop = FALSE]
expr_mat <- expr_mat[, rownames(meta), drop = FALSE]

meta$group <- droplevels(meta$group)

cat("\nFinal sample counts:\n")
table(meta$group)


# =========================================================
# 4) Filter lowly expressed genes
# =========================================================

# For microarray, filter based on mean expression
mean_expr <- rowMeans(expr_mat)
hist(mean_expr, breaks = 50, main = "Mean Expression", xlab = "Mean Expression")

# Filter genes with mean expression > 5 (adjust threshold based on your data)
filtered_expr <- expr_mat[mean_expr > 5, ]
dim(filtered_expr)

cat("\nGenes after filtering:", dim(filtered_expr)[1], "\n")

# =========================================================
# 5) Differential expression with limma
# =========================================================

# Create design matrix (0 + group = no intercept, means model)
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)
design

# Fit linear model
fit <- lmFit(filtered_expr, design)

# Define contrast: Metastasis - Primary
if (all(c("Metastasis", "Primary") %in% colnames(design))) {
  contrast <- makeContrasts(Metastasis - Primary, levels = design)
} else {
  stop("Need Primary and Metastasis groups, or change the contrast to your actual comparison.")
}

contrast

# Apply contrast
fit_contrast <- contrasts.fit(fit, contrast)

# Apply eBayes (moderated t-statistics)
fit_ebayes <- eBayes(fit_contrast)

# Get top tags (de results)
top_res <- topTable(fit_ebayes, number = Inf, adjust.method = "BH")
top_res <- top_res %>%
  rownames_to_column("gene_id")

head(top_res)

# Filter significant DEGs
sig_res <- top_res %>%
  filter(abs(logFC) > 1, adj.P.Val < 0.05)

cat("\nNumber of significant DEGs:", nrow(sig_res), "\n")
head(sig_res)

# Save DEG table
write.csv(top_res, "MSKCC2010_all_DE_results.csv", row.names = FALSE)
write.csv(sig_res, "MSKCC2010_sig_DEGs.csv", row.names = FALSE)

## ---------------------------------------------------------
## 6) Your gene IDs are already Entrez IDs!
##    Convert Entrez → SYMBOL for better readability
## ---------------------------------------------------------

gene_entrez <- sig_res$gene_id  # Already Entrez IDs

cat("First 10 gene IDs (Entrez):\n")
print(head(gene_entrez, 10))

# Convert Entrez IDs to Gene Symbols
anno <- bitr(
  gene_entrez,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

anno$ENTREZID <- as.character(anno$ENTREZID)

sig_anno <- sig_res %>%
  left_join(anno, by = c("gene_id" = "ENTREZID"))

# For ORA, use Entrez IDs directly (they're already Entrez!)
genes_entrez <- unique(gene_entrez[!is.na(gene_entrez)])

cat("\nNumber of genes for ORA:", length(genes_entrez), "\n")

sig_anno <- sig_anno %>%
  mutate(gene_symbol = SYMBOL)

write.csv(sig_anno, "MSKCC2010_sig_DEGs_annotated.csv", row.names = FALSE)

## ---------------------------------------------------------
## 7) Reactome ORA (using Entrez IDs directly)
## ---------------------------------------------------------

ora_res <- enrichPathway(
  gene = genes_entrez,
  organism = "human",
  pvalueCutoff = 0.5,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  readable = TRUE  # This will convert Entrez → Symbols for display
)
# =========================================================
# 8) Visualizations
# =========================================================

# Dotplot
p1 <- dotplot(ora_res, showCategory = 15, x = "GeneRatio") +
  ggtitle("MSKCC2010 Reactome ORA dotplot (Metastasis vs Primary)") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 8)
  )

# Enrichment map
p2 <- {
  sim <- pairwise_termsim(ora_res)
  emapplot(sim, showCategory = 15) +
    ggtitle("MSKCC2010 Reactome enrichment map") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14))
}

print(p1)
print(p2)

# Save plots
ggsave("MSKCC2010_reactome_dotplot.png", p1, width = 10, height = 7, dpi = 300)
ggsave("MSKCC2010_reactome_emap.png", p2, width = 10, height = 7, dpi = 300)

# =========================================================
# 9) Optional: top genes table
# =========================================================

top10 <- sig_anno %>%
  arrange(adj.P.Val) %>%
  select(gene_id, everything()) %>%
  head(10)

write.csv(top10, "MSKCC2010_top10_sig_genes.csv", row.names = FALSE)
cat("\nTop 10 significant genes:\n")
print(top10)

# =========================================================
# 10) Additional visualizations for DE results
# =========================================================

# Volcano plot
volcano_plot <- volcanoplot(
  fit_ebayes, 
  coef = 1,
  holdrange = 5,
  main = "MSKCC2010 Volcano Plot (Metastasis vs Primary)",
  xlab = "Log2 Fold Change",
  ylab = "Moderated t-statistic"
)

print(volcano_plot)
ggsave("MSKCC2010_volcano_plot.png", volcano_plot, width = 10, height = 7, dpi = 300)

# MA plot
ma_plot <- maPlot(fit_ebayes, coef = 1, 
                  main = "MSKCC2010 MA Plot (Metastasis vs Primary)")

print(ma_plot)
ggsave("MSKCC2010_ma_plot.png", ma_plot, width = 10, height = 7, dpi = 300)

cat("\n=== Pipeline complete! ===\n")
cat("Output files saved to:", outdir, "\n")

