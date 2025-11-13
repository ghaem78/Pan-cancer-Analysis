# SCRIPT 02: TCGA Selectivity Analysis & Downstream Case Studies
# PURPOSE: 1. Run Selectivity Score on TCGA data
#          2. Run ACC Case Study (mir-202)
#          3. Run TGCT Case Study (mir-302b)
# RUNTIME: ~15-25 minutes
# INPUT: vsd_miRNA_Normalized.rds (from Script 01)

# --- 1. LOAD LIBRARIES ---
print("Loading libraries...")
library(dplyr)
library(DESeq2)
library(TCGAbiolinks)
library.package(survival)
library(survminer)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- 2. LOAD NORMALIZED TCGA DATA (from Script 01) ---
print("Loading 'vsd_miRNA_Normalized.rds'...")
vsd <- readRDS("vsd_miRNA_Normalized.rds")

# --- 3. [SELECTIVITY ANALYSIS] ---
print("Starting Selectivity Score analysis...")

# 3.1: Remove LAML (blood cancer)
print("Removing LAML...")
clean_metadata <- as.data.frame(colData(vsd))
laml_project <- clean_metadata$project_id[clean_metadata$sample_type == "Primary Blood Derived Cancer - Peripheral Blood"][1]
vsd_solid <- vsd[, vsd$project_id != laml_project]

# 3.2: Create Mean Expression Matrix
print("Creating mean expression matrix...")
normalized_counts <- assay(vsd_solid)
metadata_solid <- as.data.frame(colData(vsd_solid))
cancer_types <- unique(metadata_solid$project_id)
mean_matrix <- matrix(NA, nrow = nrow(normalized_counts), ncol = length(cancer_types))
rownames(mean_matrix) <- rownames(normalized_counts)
colnames(mean_matrix) <- cancer_types

for (cancer in cancer_types) {
  samples_for_this_cancer <- rownames(metadata_solid[metadata_solid$project_id == cancer, ])
  mean_values <- rowMeans(normalized_counts[, samples_for_this_cancer, drop = FALSE], na.rm = TRUE)
  mean_matrix[, cancer] <- mean_values
}

# 3.3: Calculate Selectivity Score
print("Calculating Selectivity Scores...")
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
  selectivity_scores[i, 2:6] <- c(max_val, max_cancer_name, second_max_val, max_val - second_max_val)
}

# 3.4: Create and save the final ranked list
final_selective_list_tcga <- selectivity_scores[order(selectivity_scores$selectivity_score, decreasing = TRUE), ]
print("Top 10 selective miRNAs (TCGA):")
print(head(final_selective_list_tcga, 10))
saveRDS(final_selective_list_tcga, "final_selective_list_TCGA.rds")

# --- [Helper Function for Correlation] ---
run_correlation_analysis <- function(mirna_expr_vector, rna_vst_object) {
  # 1. Extract Gene Data
  rna_norm_counts <- assay(rna_vst_object)
  rna_metadata <- as.data.frame(colData(rna_vst_object))
  rna_data <- as.data.frame(t(rna_norm_counts))
  rna_data$patient_id <- substr(rna_metadata$barcode, 1, 12)
  
  # 2. Extract miRNA Data (passed as argument)
  mirna_data <- data.frame(patient_id = names(mirna_expr_vector), mir_expr = mirna_expr_vector)
  
  # 3. Merge
  cor_data <- merge(mirna_data, rna_data, by = "patient_id")
  
  # 4. Run Correlation
  print(paste("Correlating miRNA against", (ncol(cor_data)-3), "genes..."))
  mir_vector <- cor_data$mir_expr
  gene_columns <- !grepl("patient_id|mir_expr", colnames(cor_data))
  gene_matrix <- cor_data[, gene_columns]
  
  cor_results <- apply(gene_matrix, 2, function(gene_col) {
    test <- try(cor.test(mir_vector, gene_col, method = "pearson"), silent = TRUE)
    if (inherits(test, "try-error")) return(c(estimate = NA, p.value = NA))
    else return(c(estimate = test$estimate, p.value = test$p.value))
  })
  
  # 5. Format Results
  cor_results_df <- as.data.frame(t(cor_results))
  cor_results_df$gene_id <- rownames(cor_results_df)
  colnames(cor_results_df)[1:2] <- c("estimate.cor", "p.value")
  
  # 6. Suppressed List
  suppressed_list <- subset(cor_results_df, estimate.cor < 0 & p.value < 0.01)
  suppressed_list <- suppressed_list[order(suppressed_list$estimate.cor), ]
  
  # 7. Activated List
  activated_list <- subset(cor_results_df, estimate.cor > 0 & p.value < 0.01)
  activated_list <- activated_list[order(activated_list$estimate.cor, decreasing = TRUE), ]
  
  # 8. Translate Names
  suppressed_list$GeneSymbol <- mapIds(org.Hs.eg.db, keys = gsub("\\..*", "", suppressed_list$gene_id), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  activated_list$GeneSymbol <- mapIds(org.Hs.eg.db, keys = gsub("\\..*", "", activated_list$gene_id), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  
  return(list(Suppressed = suppressed_list, Activated = activated_list))
}

# --- 4. [CASE STUDY 1: hsa-mir-202 in ACC] ---
print("--- STARTING CASE STUDY 1: ACC (hsa-mir-202) ---")

# 4.1: Download & Prepare ACC RNA-Seq
print("Querying ACC RNA-Seq data...")
query_rna_acc <- GDCquery(project = "TCGA-ACC", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(query = query_rna_acc, directory = "GDCdata_ACC_RNA", method = "api", files.per.chunk = 10)
rna_acc_prepared <- GDCprepare(query = query_rna_acc, directory = "GDCdata_ACC_RNA")

# 4.2: Normalize ACC RNA-Seq
print("Normalizing ACC RNA-Seq data...")
dds_rna_acc <- DESeqDataSet(rna_acc_prepared, design = ~ 1)
vsd_rna_acc <- vst(dds_rna_acc, blind = TRUE)

# 4.3: Prepare ACC mir-202 vector
acc_samples <- rownames(subset(clean_metadata, project_id == "TCGA-ACC"))
mir_202_vector <- assay(vsd)["hsa-mir-202", acc_samples]
names(mir_202_vector) <- substr(acc_samples, 1, 12) # Name vector with Patient IDs

# 4.4: Run Correlation
acc_correlation_results <- run_correlation_analysis(mir_202_vector, vsd_rna_acc)
print("Top 10 Suppressed Targets in ACC:")
print(head(acc_correlation_results$Suppressed, 10))
print("Top 10 Activated Targets in ACC:")
print(head(acc_correlation_results$Activated, 10))
saveRDS(acc_correlation_results, "ACC_target_lists.rds")

# 4.5: Run Survival Analysis
print("Running ACC Survival Analysis...")
clinical_data_acc <- GDCquery_clinic(project = "TCGA-ACC", type = "clinical")
survival_info <- clinical_data_acc[, c("submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up")]
survival_info$Time <- ifelse(survival_info$vital_status == "Dead", as.numeric(survival_info$days_to_death), as.numeric(survival_info$days_to_last_follow_up))
survival_info$Event <- ifelse(survival_info$vital_status == "Dead", 1, 0)

mirna_data_acc_survival <- data.frame(patient_id = names(mir_202_vector), mir_202_expr = mir_202_vector)
survival_analysis_data <- merge(mirna_data_acc_survival, survival_info, by.x = "patient_id", by.y = "submitter_id")
survival_analysis_data <- na.omit(survival_analysis_data)
median_expr <- median(survival_analysis_data$mir_202_expr)
survival_analysis_data$mir_group <- ifelse(survival_analysis_data$mir_202_expr >= median_expr, "High_mir-202", "Low_mir-202")

fit <- survfit(Surv(Time, Event) ~ mir_group, data = survival_analysis_data)
km_plot <- ggsurvplot(fit, data = survival_analysis_data, pval = TRUE, risk.table = TRUE, title = "Kaplan-Meier Plot for mir-202 in TCGA-ACC")
pdf("ACC_survival_plot.pdf")
print(km_plot)
dev.off()
print(paste("ACC Survival Plot (p=", surv_pvalue(fit, data = survival_analysis_data)$pval, ") saved."))


# --- 5. [CASE STUDY 2: hsa-mir-302b in TGCT] ---
print("--- STARTING CASE STUDY 2: TGCT (hsa-mir-302b) ---")

# 5.1: Download & Prepare TGCT RNA-Seq
# (This is the ~661MB download)
print("Querying TGCT RNA-Seq data...")
query_rna_tgct <- GDCquery(project = "TCGA-TGCT", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(query = query_rna_tgct, directory = "GDCdata_TGCT_RNA", method = "api", files.per.chunk = 1) # Use chunk=1 for stability
rna_tgct_prepared <- GDCprepare(query = query_rna_tgct, directory = "GDCdata_TGCT_RNA")

# 5.2: Normalize TGCT RNA-Seq
print("Normalizing TGCT RNA-Seq data...")
dds_rna_tgct <- DESeqDataSet(rna_tgct_prepared, design = ~ 1)
vsd_rna_tgct <- vst(dds_rna_tgct, blind = TRUE)

# 5.3: Prepare TGCT mir-302b vector
tgct_samples <- rownames(subset(clean_metadata, project_id == "TCGA-TGCT"))
mir_302b_vector <- assay(vsd)["hsa-mir-302b", tgct_samples]
names(mir_302b_vector) <- substr(tgct_samples, 1, 12) # Name vector with Patient IDs

# 5.4: Run Correlation
tgct_correlation_results <- run_correlation_analysis(mir_302b_vector, vsd_rna_tgct)
print("Top 10 Suppressed Targets in TGCT:")
print(head(tgct_correlation_results$Suppressed, 10))
print("Top 10 Activated Targets in TGCT:")
print(head(tgct_correlation_results$Activated, 10))
saveRDS(tgct_correlation_results, "TGCT_target_lists.rds")

print("--- SCRIPT 02 COMPLETE ---")
