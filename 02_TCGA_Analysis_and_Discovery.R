# SCRIPT 02: TCGA Selectivity & Case Studies (ACC & TGCT)
# PURPOSE: Run Selectivity Score, analyze mir-202 in ACC, and mir-302b in TGCT.

# --- 1. LOAD LIBRARIES ---
library(dplyr)
library(DESeq2)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- 2. LOAD NORMALIZED DATA ---
# (Assumes you ran Script 01)
if(file.exists("vsd_miRNA_Normalized.rds")){
  vsd <- readRDS("vsd_miRNA_Normalized.rds")
} else {
  stop("Please run Script 01 first to generate 'vsd_miRNA_Normalized.rds'")
}

# --- 3. SELECTIVITY ANALYSIS ---
print("Starting Selectivity Score analysis...")
clean_metadata <- as.data.frame(colData(vsd))
# Remove LAML for solid tumor comparison
laml_project <- clean_metadata$project_id[clean_metadata$sample_type == "Primary Blood Derived Cancer - Peripheral Blood"][1]
vsd_solid <- vsd[, vsd$project_id != laml_project]

# Create Mean Expression Matrix
normalized_counts <- assay(vsd_solid)
metadata_solid <- as.data.frame(colData(vsd_solid))
cancer_types <- unique(metadata_solid$project_id)
mean_matrix <- matrix(NA, nrow = nrow(normalized_counts), ncol = length(cancer_types))
rownames(mean_matrix) <- rownames(normalized_counts)
colnames(mean_matrix) <- cancer_types

for (cancer in cancer_types) {
  samples <- rownames(metadata_solid[metadata_solid$project_id == cancer, ])
  mean_matrix[, cancer] <- rowMeans(normalized_counts[, samples, drop = FALSE], na.rm = TRUE)
}

# Calculate Score
selectivity_scores <- data.frame(miRNA = rownames(mean_matrix), max_expr = NA, max_cancer = NA, second_max_expr = NA, selectivity_score = NA)
for (i in 1:nrow(mean_matrix)) {
  sorted <- sort(mean_matrix[i, ], decreasing = TRUE)
  selectivity_scores[i, 2:5] <- c(sorted[1], names(sorted)[1], sorted[2], sorted[1] - sorted[2])
}
final_list <- selectivity_scores[order(selectivity_scores$selectivity_score, decreasing = TRUE), ]
saveRDS(final_list, "final_selective_list_TCGA.rds")
print("Selectivity analysis complete.")

# --- 4. CASE STUDY 1: ACC (mir-202) ---
print("--- Running ACC Case Study ---")
# 4.1 RNA-Seq Download & Norm
query_acc <- GDCquery(project = "TCGA-ACC", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(query_acc, method = "api", files.per.chunk = 10)
acc_rna_prep <- GDCprepare(query_acc)
vsd_rna_acc <- vst(DESeqDataSet(acc_rna_prep, design = ~1), blind = TRUE)

# 4.2 Correlation
# (Helper function for correlation)
run_cor <- function(mir_vec, rna_vst) {
  rna_dat <- as.data.frame(t(assay(rna_vst)))
  rna_dat$pid <- substr(colData(rna_vst)$barcode, 1, 12)
  mir_dat <- data.frame(pid = names(mir_vec), mir = mir_vec)
  merged <- merge(mir_dat, rna_dat, by = "pid")
  
  res <- apply(merged[, !grepl("pid|mir", colnames(merged))], 2, function(x) {
    test <- try(cor.test(merged$mir, x), silent=TRUE)
    if(inherits(test, "try-error")) c(NA, NA) else c(test$estimate, test$p.value)
  })
  df <- as.data.frame(t(res)); colnames(df) <- c("cor", "p"); df$id <- rownames(df)
  return(df)
}

acc_samples <- rownames(subset(clean_metadata, project_id == "TCGA-ACC"))
mir202 <- assay(vsd)["hsa-mir-202", acc_samples]; names(mir202) <- substr(acc_samples, 1, 12)
acc_res <- run_cor(mir202, vsd_rna_acc)

# Filter & Save
acc_suppressed <- subset(acc_res, cor < 0 & p < 0.01)
acc_activated <- subset(acc_res, cor > 0 & p < 0.01)
# (Mapping IDs code skipped for brevity, assume standard mapIds usage)
saveRDS(list(suppressed=acc_suppressed, activated=acc_activated), "ACC_results.rds")

# 4.3 Survival
clin_acc <- GDCquery_clinic("TCGA-ACC", "clinical")
clin_acc$time <- ifelse(clin_acc$vital_status == "Dead", clin_acc$days_to_death, clin_acc$days_to_last_follow_up)
clin_acc$event <- ifelse(clin_acc$vital_status == "Dead", 1, 0)
surv_dat <- merge(data.frame(pid=names(mir202), expr=mir202), clin_acc, by.x="pid", by.y="submitter_id")
surv_dat$group <- ifelse(surv_dat$expr >= median(surv_dat$expr), "High", "Low")
pdf("ACC_Survival.pdf")
print(ggsurvplot(survfit(Surv(time, event) ~ group, data=surv_dat), pval=TRUE))
dev.off()

# --- 5. CASE STUDY 2: TGCT (mir-302b) ---
print("--- Running TGCT Case Study ---")
# 5.1 RNA-Seq Download & Norm (Manual fix for large download may be needed locally)
# Note: In the paper we used method="api" with small chunks or "client"
query_tgct <- GDCquery(project = "TCGA-TGCT", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
# GDCdownload(query_tgct, method = "api", files.per.chunk = 1) # Uncomment to run
# tgct_rna_prep <- GDCprepare(query_tgct)
# vsd_rna_tgct <- vst(DESeqDataSet(tgct_rna_prep, design = ~1), blind = TRUE)

# 5.2 Correlation
# tgct_samples <- rownames(subset(clean_metadata, project_id == "TCGA-TGCT"))
# mir302 <- assay(vsd)["hsa-mir-302b", tgct_samples]; names(mir302) <- substr(tgct_samples, 1, 12)
# tgct_res <- run_cor(mir302, vsd_rna_tgct)
# Save results...
print("Script 02 Complete")
