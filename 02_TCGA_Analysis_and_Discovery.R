# SCRIPT 02: TCGA Selectivity Analysis & Downstream Case Studies (FINAL VERSION)
#
# PURPOSE: 
#   1. Perform Selectivity Score analysis on TCGA miRNA data (Solid Tumors).
#   2. Case Study I: ACC (mir-202) - Mechanism (Correlation) & Clinical (Survival).
#   3. Case Study II: TGCT (mir-302b) - Mechanism (Correlation).
#
# INPUT: vsd_miRNA_Normalized.rds (Output from Script 01)
# OUTPUTS: final_selective_list_TCGA.rds, ACC_results.rds, TGCT_results.rds, Plots

# --- 1. LOAD LIBRARIES ---
print("Loading libraries...")
library(dplyr)
library(DESeq2)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- 2. LOAD NORMALIZED TCGA DATA ---
print("Loading 'vsd_miRNA_Normalized.rds'...")
if(file.exists("vsd_miRNA_Normalized.rds")){
  vsd <- readRDS("vsd_miRNA_Normalized.rds")
} else {
  stop("Error: 'vsd_miRNA_Normalized.rds' not found. Please run Script 01 first.")
}

# --- 3. SELECTIVITY SCORE ANALYSIS ---
print("Starting Selectivity Analysis...")

# 3.1: Remove LAML (Leukemia) to focus on solid tumors
# (Crucial step to avoid skewing selectivity scores)
print("Filtering out LAML samples...")
vsd_solid <- vsd[, vsd$project_id != "TCGA-LAML"]

# 3.2: Create Mean Expression Matrix
print("Calculating Mean Expression Matrix...")
normalized_counts <- assay(vsd_solid)
metadata_solid <- as.data.frame(colData(vsd_solid))
cancer_types <- unique(metadata_solid$project_id)

mean_matrix <- matrix(NA, nrow = nrow(normalized_counts), ncol = length(cancer_types))
rownames(mean_matrix) <- rownames(normalized_counts)
colnames(mean_matrix) <- cancer_types

for (cancer in cancer_types) {
  samples <- rownames(metadata_solid[metadata_solid$project_id == cancer, ])
  # drop=FALSE ensures it works even if a cancer has only 1 sample
  mean_matrix[, cancer] <- rowMeans(normalized_counts[, samples, drop = FALSE], na.rm = TRUE)
}
# Save mean matrix (Useful for Heatmaps/Figure S1)
saveRDS(mean_matrix, "02_mean_expression_matrix.rds")

# 3.3: Calculate Selectivity Score
print("Calculating Scores (Max - SecondMax)...")
selectivity_scores <- data.frame(
  miRNA = rownames(mean_matrix), max_expr = NA, max_cancer = NA,
  second_max_expr = NA, selectivity_score = NA
)

for (i in 1:nrow(mean_matrix)) {
  expr_row <- mean_matrix[i, ]
  sorted_expr <- sort(expr_row, decreasing = TRUE, na.last = TRUE)
  
  max_val <- sorted_expr[1]
  second_max_val <- sorted_expr[2]
  max_cancer_name <- names(sorted_expr)[1]
  
  selectivity_scores[i, 2:5] <- c(max_val, max_cancer_name, second_max_val, max_val - second_max_val)
}

# 3.4: Rank and Save
final_selective_list_tcga <- selectivity_scores[order(selectivity_scores$selectivity_score, decreasing = TRUE), ]
print("Top 10 Selective miRNAs (TCGA):")
print(head(final_selective_list_tcga, 10))
saveRDS(final_selective_list_tcga, "final_selective_list_TCGA.rds")


# --- [HELPER FUNCTION: Correlation Analysis & Gene Mapping] ---
run_cor_analysis <- function(mir_vec, rna_vst_obj) {
  # Prepare Gene Data
  rna_dat <- as.data.frame(t(assay(rna_vst_obj)))
  rna_dat$pid <- substr(colData(rna_vst_obj)$barcode, 1, 12)
  
  # Prepare miRNA Data
  mir_dat <- data.frame(pid = names(mir_vec), mir = mir_vec)
  
  # Merge by Patient ID
  merged <- merge(mir_dat, rna_dat, by = "pid")
  print(paste("  -> Common samples found:", nrow(merged)))
  
  # Run Correlation
  print("  -> Computing Pearson correlations (this may take a moment)...")
  gene_cols <- merged[, !grepl("pid|mir", colnames(merged))]
  
  res <- apply(gene_cols, 2, function(x) {
    test <- try(cor.test(merged$mir, x, method="pearson"), silent=TRUE)
    if(inherits(test, "try-error")) c(NA, NA) else c(test$estimate, test$p.value)
  })
  
  df <- as.data.frame(t(res))
  colnames(df) <- c("cor", "p")
  df$ensembl_id <- rownames(df)
  
  # Map Ensembl IDs to Gene Symbols
  print("  -> Mapping Gene Symbols...")
  clean_ids <- gsub("\\..*", "", df$ensembl_id) # Remove version numbers
  df$Symbol <- mapIds(org.Hs.eg.db, keys=clean_ids, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  
  return(list(results = df, merged_data = merged))
}


# --- 4. CASE STUDY I: ACC (mir-202) ---
print("--- STARTING CASE STUDY I: ACC (mir-202) ---")

# 4.1: Download & Normalize ACC RNA-Seq
# (Checks if data exists to avoid re-downloading)
query_acc <- GDCquery(project = "TCGA-ACC", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
if(!file.exists("GDCdata/TCGA-ACC")) {
    print("  -> Downloading ACC RNA-Seq...")
    GDCdownload(query_acc, method = "api", files.per.chunk = 10)
}
print("  -> Preparing & Normalizing ACC Data...")
acc_rna_prep <- GDCprepare(query_acc)
vsd_rna_acc <- vst(DESeqDataSet(acc_rna_prep, design = ~1), blind = TRUE)

# 4.2: Correlation Analysis
acc_samples <- rownames(subset(as.data.frame(colData(vsd)), project_id == "TCGA-ACC"))
mir202 <- assay(vsd)["hsa-mir-202", acc_samples]
names(mir202) <- substr(acc_samples, 1, 12)

acc_analysis <- run_cor_analysis(mir202, vsd_rna_acc)
acc_res <- acc_analysis$results

# Save Results
acc_suppressed <- subset(acc_res, cor < 0 & p < 0.01) %>% arrange(cor)
acc_activated <- subset(acc_res, cor > 0 & p < 0.01) %>% arrange(desc(cor))

print("  -> Top Suppressed Targets (ACC):")
print(head(acc_suppressed[, c("Symbol", "cor", "p")], 5))
saveRDS(list(suppressed=acc_suppressed, activated=acc_activated), "ACC_results_with_symbols.rds")
saveRDS(acc_analysis$merged_data, "ACC_final_correlation_data.rds")

# 4.3: Survival Analysis
print("  -> Running Survival Analysis...")
clin_acc <- GDCquery_clinic("TCGA-ACC", "clinical")
clin_acc$time <- ifelse(clin_acc$vital_status == "Dead", clin_acc$days_to_death, clin_acc$days_to_last_follow_up)
clin_acc$event <- ifelse(clin_acc$vital_status == "Dead", 1, 0)

surv_dat <- merge(data.frame(pid=names(mir202), expr=mir202), clin_acc, by.x="pid", by.y="submitter_id")
surv_dat$group <- ifelse(surv_dat$expr >= median(surv_dat$expr), "High", "Low")
surv_dat <- na.omit(surv_dat)

pdf("Figure4_ACC_Survival.pdf")
print(ggsurvplot(survfit(Surv(time, event) ~ group, data=surv_dat), pval=TRUE, title="ACC mir-202 Survival"))
dev.off()
print("  -> Survival Plot Saved.")


# --- 5. CASE STUDY II: TGCT (mir-302b) ---
print("--- STARTING CASE STUDY II: TGCT (mir-302b) ---")

# 5.1: Download & Normalize TGCT RNA-Seq
query_tgct <- GDCquery(project = "TCGA-TGCT", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")

# Note: If downloading for the first time, uncomment GDCdownload
# GDCdownload(query_tgct, method = "api", files.per.chunk = 5) 

# Assuming data is downloaded in specific folder (based on our troubleshooting):
tgct_data_dir <- if(dir.exists("GDCdata_TGCT_RNA")) "GDCdata_TGCT_RNA" else "GDCdata"
print(paste("  -> Loading TGCT data from:", tgct_data_dir))

tgct_rna_prep <- GDCprepare(query_tgct, directory = tgct_data_dir)
vsd_rna_tgct <- vst(DESeqDataSet(tgct_rna_prep, design = ~1), blind = TRUE)

# 5.2: Correlation Analysis
tgct_samples <- rownames(subset(as.data.frame(colData(vsd)), project_id == "TCGA-TGCT"))
mir302 <- assay(vsd)["hsa-mir-302b", tgct_samples]
names(mir302) <- substr(tgct_samples, 1, 12)

tgct_analysis <- run_cor_analysis(mir302, vsd_rna_tgct)
tgct_res <- tgct_analysis$results

# Save Results
tgct_activated <- subset(tgct_res, cor > 0 & p < 0.01) %>% arrange(desc(cor))
print("  -> Top Activated Targets (TGCT):")
print(head(tgct_activated[, c("Symbol", "cor", "p")], 5))

saveRDS(tgct_res, "TGCT_results_with_symbols.rds")
saveRDS(tgct_analysis$merged_data, "TGCT_final_correlation_data.rds")

print("--- SCRIPT 02 COMPLETE ---")
