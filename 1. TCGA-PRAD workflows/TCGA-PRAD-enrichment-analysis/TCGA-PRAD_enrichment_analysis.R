# =================================================================================================
# Differential Expression Analysis and Enrichment Analysis Between Normal Tissue and Primary Tumor
# =================================================================================================

library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyr)
library(openxlsx)
library(tidyverse)
library(ggpubr)
library(edgeR)
library(DESeq2)
library(biomaRt)
library(ReactomePA)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)

outdir <- "output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

excel_sheets <- list()

tcga_obj <- readRDS("TCGA_panCancer_FOLH1.rds")

count_mat <- assay(tcga_obj)
meta <- as.data.frame(colData(tcga_obj))
meta$barcode <- rownames(meta)

# Extract ONLY TCGA-PRAD samples
meta_prad <- meta[meta$project == "TCGA-PRAD", ]
prad_count_mat <- count_mat[, meta_prad$barcode]
expr_prad <- prad_count_mat

## Build a grouping variable
meta_prad$group <- factor(
  meta_prad$sample_type,
  levels = c("Solid Tissue Normal", "Primary Tumor"),
  labels = c("Normal", "Tumor")
)

table(meta_prad$group)

## Keep only Normal and Tumor samples
keep <- meta_prad$group %in% c("Normal", "Tumor")
meta_prad <- meta_prad[keep, , drop = FALSE]
expr_prad <- expr_prad[, rownames(meta_prad), drop = FALSE]

meta_prad$group <- droplevels(meta_prad$group)

# ----------------------------------
## 1) Filter lowly expressed genes
## ---------------------------------
dge <- DGEList(counts = expr_prad)
dge <- calcNormFactors(dge)

cpm_mat <- cpm(dge)
mean_log2_cpm <- rowMeans(log2(cpm_mat + 1))
hist(mean_log2_cpm, breaks = 50, main = "Mean log2 CPM", xlab = "Mean log2 CPM")

filtered_expr <- expr_prad[mean_log2_cpm > 1, ]
dim(filtered_expr)

# ---------------------------------------
## 2) Differential expression with edgeR
## --------------------------------------
dge <- DGEList(counts = filtered_expr)
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + group, data = meta_prad)
colnames(design) <- levels(meta_prad$group)
design

dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)

if (all(c("Tumor", "Normal") %in% colnames(design))) {
  contrast <- makeContrasts(Tumor - Normal, levels = design)
}

qlf <- glmQLFTest(fit, contrast = contrast)
top_res <- topTags(qlf, n = Inf)$table

sig_res <- top_res %>%
  rownames_to_column("gene_id") %>%
  filter(abs(logFC) > 1, FDR < 0.05)

sig_res$gene_id <- sub("\\..*", "", sig_res$gene_id)

nrow(sig_res)
head(sig_res)

excel_sheets[["TCGA_PRAD_all_DE_results"]] <- top_res
excel_sheets[["TCGA_PRAD_sig_DEGs"]] <- sig_res

# -------------------------------
## 3) Map gene IDs to Entrez IDs
## ------------------------------

gene_ids <- sig_res$gene_id

if (all(grepl("^ENSG", gene_ids))) {
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  anno <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
    filters = "ensembl_gene_id",
    values = gene_ids,
    mart = mart
  )
  
  anno$ensembl_gene_id <- as.character(anno$ensembl_gene_id)
  
  sig_anno <- sig_res %>%
    left_join(anno, by = c("gene_id" = "ensembl_gene_id"))
  
  genes_entrez <- sig_anno$entrezgene_id
  
} else {
  anno <- bitr(
    gene_ids,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  anno$SYMBOL <- as.character(anno$SYMBOL)
  
  sig_anno <- sig_res %>%
    left_join(anno, by = c("gene_id" = "SYMBOL"))
  
  genes_entrez <- sig_anno$ENTREZID
}

genes_entrez <- genes_entrez[!is.na(genes_entrez)]
genes_entrez <- unique(genes_entrez)

length(genes_entrez)

excel_sheets[["TCGA_PRAD_sig_DEGs_annotated"]] <- sig_anno

# ------------------
## 4) Reactome ORA
## -----------------
ora_res <- enrichPathway(
  gene = genes_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

ora_df <- as.data.frame(ora_res)
excel_sheets[["TCGA_PRAD_reactome_ORA"]] <- ora_df

head(ora_df[, c("Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue")])

# --------------------
## 5) Visualizations
## -------------------
p1 <- dotplot(ora_res, showCategory = 15, x = "GeneRatio") +
  ggtitle("TCGA-PRAD Reactome ORA dotplot") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 8)
  )

p2 <- {
  sim <- pairwise_termsim(ora_res)
  emapplot(sim, showCategory = 15) +
    ggtitle("TCGA-PRAD Reactome enrichment map") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14))
}

print(p1)
print(p2)

ggsave(file.path(outdir, "TCGA_PRAD_reactome_dotplot.png"), p1, width = 10, height = 7, dpi = 300)
ggsave(file.path(outdir, "TCGA_PRAD_reactome_emap.png"), p2, width = 10, height = 7, dpi = 300)

# -------------------------------
## 6) Optional: top genes table
## ------------------------------
top10 <- sig_anno %>%
  arrange(FDR) %>%
  dplyr::select(gene_id, dplyr::everything()) %>%
  head(10)

excel_sheets[["TCGA_PRAD_top10_sig_genes"]] <- top10
top10

# ---------------------------
## 7) Write Excel workbook
# ---------------------------

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

wb <- createWorkbook()
existing_names <- character()

for (nm in names(excel_sheets)) {
  sn <- safe_sheet_name(nm, existing_names)
  addWorksheet(wb, sn)
  writeData(wb, sheet = sn, x = excel_sheets[[nm]], rowNames = FALSE)
  existing_names <- c(existing_names, sn)
}

saveWorkbook(
  wb,
  file = file.path(outdir, "TCGA_PRAD_DEG_and_Reactome_results.xlsx"),
  overwrite = TRUE
)