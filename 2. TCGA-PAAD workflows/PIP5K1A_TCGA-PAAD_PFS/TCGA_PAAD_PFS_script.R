library(SummarizedExperiment)

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

#==============================================
# Normalizing and variance stabalizing PAAD data
# ==============================================

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = paad_counts,
                              colData   = meta_paad,
                              design    = ~ 1)

# Filter low count genes
keep <- rowSums(counts(dds) >= 10) >= 0.2 * ncol(dds)
dds  <- dds[keep, ]

dds  <- estimateSizeFactors(dds)
vsd  <- varianceStabilizingTransformation(dds, blind = TRUE)
expr_vst <- assay(vsd)  # genes x samples

# ==================
# Extracting PIP5K1A
# ==================

# Find PIP5K1A row
pip_row <- which(gene_annot_paad$gene_name == "PIP5K1A")
if (length(pip_row) != 1) {
  stop("PIP5K1A not found or not unique in gene_annot_paad$gene_name")
}

pip5k1a_expr <- expr_vst[pip_row, ]
# Make sure alignment is sample-wise
pip5k1a_expr <- pip5k1a_expr[colnames(expr_vst)]  # just to be explicit

meta_paad$PIP5K1A_expr <- pip5k1a_expr[rownames(meta_paad)]


sum(!is.na(meta_paad$paper_ABSOLUTE.Purity.1))
summary(meta_paad$paper_ABSOLUTE.Purity.1)
range(meta_paad$paper_ABSOLUTE.Purity.1, na.rm = TRUE)


# ===============================================
# Create 4 groups for a single parameter
# ==============================================

library(survival)
library(survminer)

# Choose parameter column
param_name <- "paper_ABSOLUTE.Purity.1"   # <-- edit to match your meta_paad column
param_vec  <- meta_paad[[param_name]]

# Remove samples with missing values in either expression or parameter or PFS
valid_idx <- complete.cases(meta_paad$PIP5K1A_expr,
                            param_vec,
                            meta_paad$PFS_time,
                            meta_paad$PFS_event)
meta_sub  <- meta_paad[valid_idx, ]
pip_vec   <- meta_sub$PIP5K1A_expr
par_vec   <- meta_sub[[param_name]]

# Define cutoffs (median split; you can change to quantile, etc.)
pip_cut <- median(pip_vec, na.rm = TRUE)
par_cut <- median(par_vec, na.rm = TRUE)

pip_group <- ifelse(pip_vec >= pip_cut, "HighPIP5K1A", "LowPIP5K1A")
par_group <- ifelse(par_vec >= par_cut, "HighParam",   "LowParam")

# Combine into 4 groups: 1, 2, 3, 4 as you described
combo_group <- character(length(pip_group))
combo_group[pip_group == "LowPIP5K1A"  & par_group == "LowParam"]  <- "1_LowGene_LowParam"
combo_group[pip_group == "LowPIP5K1A"  & par_group == "HighParam"] <- "2_LowGene_HighParam"
combo_group[pip_group == "HighPIP5K1A" & par_group == "LowParam"]  <- "3_HighGene_LowParam"
combo_group[pip_group == "HighPIP5K1A" & par_group == "HighParam"] <- "4_HighGene_HighParam"

meta_sub$combo_group <- factor(combo_group,
                               levels = c("1_LowGene_LowParam",
                                          "2_LowGene_HighParam",
                                          "3_HighGene_LowParam",
                                          "4_HighGene_HighParam"))

# ==================================
# 3. Build PFS object and KM plot
# =================================

# PFS survival object (make sure column names match your data)
Surv_PFS <- with(meta_sub, Surv(PFS_time, PFS_event))

fit_km <- survfit(Surv_PFS ~ combo_group, data = meta_sub)

p <- ggsurvplot(
  fit_km,
  data        = meta_sub,
  pval        = TRUE,
  risk.table  = TRUE,
  legend.title = paste0("PIP5K1A + ", param_name),
  legend.labs  = c("Low gene + low param",
                   "Low gene + high param",
                   "High gene + low param",
                   "High gene + high param"),
  palette     = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
)

# Save
ggsave(file.path(outdir, paste0("KM_PFS_PIP5K1A_", param_name, "_4groups.png")),
       p$plot, width = 6, height = 5, dpi = 300)





cols <- c(
  "days_to_collection",
  "days_to_diagnosis",
  "days_to_last_follow_up",
  "follow_ups_disease_response",
  "days_to_birth",
  "cause_of_death",
  "days_to_death",
  "days_to_sample_procurement",
  "progression_or_recurrence",
  "year_of_death",
  "cause_of_death_source",
  "alcohol_days_per_week",
  "paper_Survival..months.",
  "paper_days_to_birth",
  "paper_days_to_death",
  "paper_days_to_last_followup",
  "paper_Max.Days.of.Follow.up.from.Initial.Surgery",
  "paper_Event.Free.Days..Aggressive.disease.",
  "paper_Event.Free.Days..Distant.Metastasis.",
  "paper_new_tumor_event_after_initial_treatment",
  "paper_days.to.death.or.last.contact",
  "paper_CLIN.days_to_last_followup",
  "paper_CLIN.days_to_last_known_alive",
  "paper_CLIN.days_to_death",
  "paper_Days.to.last.known.alive",
  "paper_New.tumor.event.after.initial.treatment",
  "paper_Days.from.surgery.to.last.followup",
  "paper_Days.from.surgery.to.death",
  "paper_Max.Follow.up",
  "paper_New.Tumor.After.Initial.Treatment..days.",
  "paper_days_to_last_known_alive",
  "paper_Days.to.last.followup",
  "paper_Days.until.death",
  "paper_Combined.days.to.last.followup.or.death",
  "paper_time_of_follow.up",
  "paper_Height..at.time.of.diagnosis...cm.",
  "paper_Weight..at.time.of.diagnosis...kg.",
  "paper_MET.splicing.event..Brad.Murray..Angela.Brooks.",
  "paper_CURATED_DAYS_TO_DEATH_OR_LAST_FU",
  "paper_CURATED_TCGA_DAYS_TO_DEATH_OR_LAST_FU",
  "paper_OS.days",
  "paper_new.tumor.event.dx.indicator",
  "paper_days.to.new.tumor.event",
  "paper_RFS.days",
  "paper_treatment.outcome.at.tcga.followup",
  "paper_local.recurrence.days",
  "paper_distant.recurrence.days",
  "paper_Survival",
  "paper_Follow.up.cause.of.death",
  "paper_Follow.up.vital.status",
  "paper_Days.to.death",
  "paper_Follow.up.tumor.status",
  "paper_Follow.up.days",
  "paper_Death..Metastasis",
  "paper_Follow_up",
  "paper_Most_Recent_Days_to_Follow_up",
  "paper_Follow_up_New_Tumor_Event",
  "paper_os_days",
  "paper_recurred_progressed",
  "paper_pfs_days",
  "paper_Days.to.Last.Follow.up",
  "paper_Days.to.Death",
  "paper_Days.to.Recurrence"
)

sapply(meta_paad[, cols], function(x) sum(!is.na(x)))
