# =========================================================
# MSKCC2010: FOLH1 vs Biomarkers (Metastasis Only)
# =========================================================

library(MultiAssayExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(survival)
library(survminer)
library(pheatmap)

# =========================================================
# OUTPUT
# =========================================================

outdir <- "MSKCC_FOLH1_metastasis"
km_dir <- file.path(outdir, "KM_plots")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(km_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# LOAD DATA
# =========================================================

data <- readRDS("MSKCC2010.rds")

rna <- experiments(data)[["mrna_agilent_microarray"]]

expr_mat <- assay(rna)
mode(expr_mat) <- "numeric"

meta <- as.data.frame(colData(data))
meta$sample_id <- rownames(meta)

# =========================================================
# METASTASIS FILTER
# =========================================================

meta$group <- case_when(
  grepl("met", meta$SAMPLE_TYPE, ignore.case = TRUE) ~ "Metastatic",
  grepl("primary", meta$SAMPLE_TYPE, ignore.case = TRUE) ~ "Primary",
  TRUE ~ NA_character_
)

meta_meta <- meta %>%
  filter(group == "Metastatic")

samples <- intersect(meta_meta$sample_id, colnames(expr_mat))

expr_mat <- expr_mat[, samples, drop = FALSE]
meta_meta <- meta_meta[match(samples, meta_meta$sample_id), ]

# =========================================================
# GENE MAPPING (ENTREZ → SYMBOL)
# =========================================================

gene_map <- c(
  FOLH1="2346", KLK3="354", AR="367", NKX3_1="4824", TMPRSS2="7113",
  AMACR="23600", ERG="2078", MKI67="4288", CDH1="999", VIM="7431",
  AURKA="6790", BRCA1="672", BRCA2="675", MYC="4609", TP53="7157",
  PTEN="5728", RB1="5925", TTF1="7080", INSM1="3642", NKX2_1="7080",
  ACP3="55", CHGA="1113", CHGB="1114", TFRC="7037", NCAM1="4684",
  SCG2="7857", SYP="6855"
)

gene_map <- gene_map[gene_map %in% rownames(expr_mat)]

expr_sub <- expr_mat[gene_map, , drop = FALSE]

rownames(expr_sub) <- names(gene_map)

# sanity check
stopifnot("FOLH1" %in% rownames(expr_sub))

# =========================================================
# STEP 1: CORRELATION
# =========================================================

target <- expr_sub["FOLH1", ]

cor_df <- lapply(rownames(expr_sub), function(g) {
  
  x <- as.numeric(expr_sub[g, ])
  y <- as.numeric(expr_sub["FOLH1", ])
  
  ct <- cor.test(x, y, method = "spearman")
  
  data.frame(
    gene = g,
    rho = unname(ct$estimate),
    p = ct$p.value
  )
})

cor_df <- bind_rows(cor_df)
cor_df$p.adj <- p.adjust(cor_df$p, "BH")

cor_df <- cor_df %>%
  arrange(desc(rho))

write.csv(cor_df,
          file.path(outdir, "FOLH1_correlations.csv"),
          row.names = FALSE)

# =========================================================
# STEP 2: BARPLOT
# =========================================================

plot_df <- cor_df %>%
  filter(!is.na(rho)) %>%
  arrange(rho) %>%
  mutate(
    gene = factor(gene, levels = gene),
    direction = ifelse(rho >= 0, "Positive", "Negative")
  )

p_bar <- ggplot(plot_df, aes(x = gene, y = rho, fill = direction)) +
  geom_col(width = 0.82, color = NA) +
  coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.8, colour = "black") +
  scale_fill_manual(
    values = c(
      "Positive" = "#E64B35FF",
      "Negative" = "#4DBBD5FF"
    )
  ) +
  scale_y_continuous(
    breaks = pretty(plot_df$rho, n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = "Correlation of FOLH1 with known biomarkers",
    x = NULL,
    y = "Spearman correlation with FOLH1"
  ) +
  theme_classic(base_size = 16, base_family = "Arial") +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.line.y = element_line(linewidth = 0.8, colour = "black"),
    axis.line.x = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.7, colour = "black"),
    axis.text = element_text(colour = "black", size = 13),
    axis.title = element_text(face = "bold", colour = "black", size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    legend.position = "none"
  )

ggsave(
  file.path(outdir, "correlation_barplot.png"),
  p_bar,
  width = 9,
  height = 7,
  dpi = 300
)

# =========================================================
# STEP 3: SCATTER PLOTS
# =========================================================

for (g in rownames(expr_sub)) {
  
  df <- data.frame(
    FOLH1 = as.numeric(expr_sub["FOLH1", ]),
    biomarker = as.numeric(expr_sub[g, ])
  )
  
  ct <- suppressWarnings(
    cor.test(df$FOLH1, df$biomarker, method = "spearman")
  )
  
  label_txt <- paste0(
    "rho = ", round(ct$estimate, 3),
    ", p = ", signif(ct$p.value, 3)
  )
  
  p <- ggplot(df, aes(x = FOLH1, y = biomarker)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = paste("FOLH1 vs", g),
      subtitle = label_txt,
      x = "FOLH1 expression",
      y = paste(g, "expression")
    )
  
  ggsave(
    file.path(outdir, paste0("scatter_", g, ".png")),
    p,
    width = 5,
    height = 5
  )
}

# =========================================================
# STEP 4: HEATMAP
# =========================================================

heat_mat <- as.matrix(expr_sub)
mode(heat_mat) <- "numeric"

pheatmap(
  heat_mat,
  scale = "row",
  show_colnames = FALSE,
  main = "FOLH1 and biomarker expression (MSKCC metastasis)",
  filename = file.path(outdir, "heatmap.png"),
  width = 8,
  height = 6
)

# =========================================================
# STEP 5: DFS SURVIVAL
# =========================================================

meta_meta <- meta_meta %>%
  mutate(
    time = suppressWarnings(as.numeric(DFS_MONTHS)),
    event = case_when(
      DFS_STATUS %in% c("1:Recurred", "1:Recurred/Progressed") ~ 1,
      DFS_STATUS %in% c("0:DiseaseFree") ~ 0,
      TRUE ~ NA_real_
    )
  )

cat("DFS missing time:", sum(is.na(meta_meta$time)), "\n")
cat("DFS missing event:", sum(is.na(meta_meta$event)), "\n")

km_results <- list()

for (g in rownames(expr_sub)) {
  
  df <- data.frame(
    expr = as.numeric(expr_sub[g, ]),
    time = meta_meta$time,
    event = meta_meta$event
  ) %>%
    filter(is.finite(expr), is.finite(time), is.finite(event))
  
  if (nrow(df) < 5) next
  
  med <- median(df$expr, na.rm = TRUE)
  df$group <- factor(ifelse(df$expr >= med, "High", "Low"), levels = c("Low", "High"))
  
  if (length(unique(df$group)) < 2) next
  
  fit <- survfit(Surv(time, event) ~ group, data = df)
  sdiff <- survdiff(Surv(time, event) ~ group, data = df)
  pval <- 1 - pchisq(sdiff$chisq, df = 1)
  
  cox_fit <- coxph(Surv(time, event) ~ group, data = df)
  cox_sum <- summary(cox_fit)
  
  hr <- unname(cox_sum$coefficients[1, "exp(coef)"])
  hr_low <- unname(cox_sum$conf.int[1, "lower .95"])
  hr_high <- unname(cox_sum$conf.int[1, "upper .95"])
  
  p <- ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    risk.table = TRUE,
    title = paste0(g, " (DFS)"),
    legend.labs = c("Low", "High")
  )
  
  ggsave(
    file.path(km_dir, paste0(g, "_KM.png")),
    p$plot,
    width = 6,
    height = 6
  )
  
  km_results[[g]] <- data.frame(
    gene = g,
    n = nrow(df),
    cutoff = med,
    n_low = sum(df$group == "Low"),
    n_high = sum(df$group == "High"),
    logrank_p = pval,
    HR_high_vs_low = hr,
    CI95_lower = hr_low,
    CI95_upper = hr_high,
    stringsAsFactors = FALSE
  )
}

km_df <- bind_rows(km_results)

if (nrow(km_df) > 0) {
  km_df <- km_df %>% arrange(logrank_p)
}

write.csv(
  km_df,
  file.path(outdir, "KM_results.csv"),
  row.names = FALSE
)

# =========================================================
# STEP 6: EXCEL EXPORT
# =========================================================

wb <- createWorkbook()

addWorksheet(wb, "Correlation")
writeData(wb, "Correlation", cor_df)

addWorksheet(wb, "KM_results")
writeData(wb, "KM_results", km_df)

saveWorkbook(
  wb,
  file.path(outdir, "MSKCC_FOLH1_results.xlsx"),
  overwrite = TRUE
)

