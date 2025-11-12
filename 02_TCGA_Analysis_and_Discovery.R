# SCRIPT 02: TCGA Selectivity & ACC Case Study
# PURPOSE: Loads normalized data from Script 01, runs Selectivity Score,
#          and performs the downstream analysis for the top hit (ACC).
# RUNTIME: ~5-10 minutes
# INPUT: vsd_miRNA_Normalized.rds

# --- 1. LOAD LIBRARIES ---
# (Install if missing)
# install.packages("dplyr")
# BiocManager::install("survival")
# BiocManager::install("survminer")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
library(dplyr)
library(DESeq2) # (Needed for colData/assay)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- 2. LOAD NORMALIZED TCGA DATA (from Script 01) ---
print("Loading 'vsd_miRNA_Normalized.rds'...")
vsd <- readRDS("vsd_miRNA_Normalized.rds")

# --- 3. [SELECTIVITY ANALYSIS] ---
print("Starting Selectivity Score analysis...")

# 3.1: Remove LAML (blood cancer) for solid tumor comparison
print("Removing LAML...")
clean_metadata <- as.data.frame(colData(vsd))
laml_project <- clean_metadata$project_id[clean_metadata$sample_type == "Primary Blood Derived Cancer - Peripheral Blood"][1]
vsd_solid <- vsd[, vsd$project_id != laml_project]

# 3.2: Create Mean Expression Matrix (like step 19)
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

# 3.3: Calculate Selectivity Score (like step 23)
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
  selectivity_scores$max_expr[i] <- max_val
  selectivity_scores$max_cancer[i] <- max_cancer_name
  selectivity_scores$second_max_expr[i] <- second_max_val
  selectivity_scores$selectivity_score[i] <- max_val - second_max_val
}

# 3.4: Create and save the final ranked list
final_selective_list_tcga <- selectivity_scores[order(selectivity_scores$selectivity_score, decreasing = TRUE), ]
print("Top 10 selective miRNAs (TCGA):")
print(head(final_selective_list_tcga, 10))
saveRDS(final_selective_list_tcga, "final_selective_list_TCGA.rds")

# --- 4. [ACC CASE STUDY: RNA-SEQ CORRELATION] ---
print("Starting ACC Case Study: Downloading RNA-Seq data...")

# 4.1: Download RNA-Seq (Gene) data for TCGA-ACC
query_rna_seq <- GDCquery(
    project = "TCGA-ACC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)
GDCdownload(query = query_rna_seq)
rna_data_prepared <- GDCprepare(query = query_rna_seq)

# 4.2: Normalize RNA-Seq data
print("Normalizing ACC RNA-Seq data...")
dds_rna <- DESeqDataSet(rna_data_prepared, design = ~ 1)
vsd_rna <- vst(dds_rna, blind = TRUE)
rna_norm_counts <- assay(vsd_rna)
rna_metadata <- as.data.frame(colData(vsd_rna))

# 4.3: Extract mir-202 data (from vsd)
acc_mirna_metadata <- subset(clean_metadata, project_id == "TCGA-ACC")
mirna_data <- data.frame(
  sample_barcode = rownames(acc_mirna_metadata),
  mir_202_expr = assay(vsd)["hsa-mir-202", rownames(acc_mirna_metadata)]
)
mirna_data$patient_id <- substr(mirna_data$sample_barcode, 1, 12)

# 4.4: Extract RNA-Seq data (genes)
rna_data <- as.data.frame(t(rna_norm_counts))
rna_data$patient_id <- substr(rna_metadata$barcode, 1, 12)

# 4.5: Merge miRNA and RNA data
final_correlation_data <- merge(mirna_data, rna_data, by = "patient_id")

# 4.6: Run Correlation Analysis (like step 35)
print("Running correlation analysis for mir-202 vs 60,000 genes...")
mir_vector <- final_correlation_data$mir_202_expr
gene_columns <- !grepl("patient_id|sample_barcode|mir_202_expr", colnames(final_correlation_data))
gene_matrix <- final_correlation_data[, gene_columns]

cor_results <- apply(gene_matrix, 2, function(gene_col) {
  test <- try(cor.test(mir_vector, gene_col, method = "pearson"), silent = TRUE)
  if (inherits(test, "try-error")) return(c(estimate = NA, p.value = NA))
  else return(c(estimate = test$estimate, p.value = test$p.value))
})

cor_results_df <- as.data.frame(t(cor_results))
cor_results_df$gene_id <- rownames(cor_results_df)
colnames(cor_results_df)[1:2] <- c("estimate.cor", "p.value") # Fix column names

# 4.7: Get Suppressed Targets (Negative Correlation)
final_target_list <- subset(cor_results_df, estimate.cor < 0 & p.value < 0.01)
final_target_list <- final_target_list[order(final_target_list$estimate.cor), ]

# 4.8: Get Activated Targets (Positive Correlation)
activated_gene_list <- subset(cor_results_df, estimate.cor > 0 & p.value < 0.01)
activated_gene_list <- activated_gene_list[order(activated_gene_list$estimate.cor, decreasing = TRUE), ]

# 4.9: Translate Gene IDs and Save
print("Translating gene IDs...")
gene_ids_cleaned_suppressed <- gsub("\\..*", "", final_target_list$gene_id)
final_target_list$GeneSymbol <- mapIds(org.Hs.eg.db, keys = gene_ids_cleaned_suppressed, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

gene_ids_cleaned_activated <- gsub("\\..*", "", activated_gene_list$gene_id)
activated_gene_list$GeneSymbol <- mapIds(org.Hs.eg.db, keys = gene_ids_cleaned_activated, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

print("Top 10 Suppressed Targets:")
print(head(final_target_list, 10))
print("Top 10 Activated Targets:")
print(head(activated_gene_list, 10))
saveRDS(final_target_list, "ACC_target_list_suppressed.rds")
saveRDS(activated_gene_list, "ACC_target_list_activated.rds")


# --- 5. [ACC CASE STUDY: SURVIVAL ANALYSIS] ---
print("Running ACC Survival Analysis...")

# 5.1: Download Clinical Data
clinical_data_acc <- GDCquery_clinic(project = "TCGA-ACC", type = "clinical")

# 5.2: Prepare Survival Info
survival_info <- clinical_data_acc[, c("submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up")]
survival_info$Time <- ifelse(survival_info$vital_status == "Dead", as.numeric(survival_info$days_to_death), as.numeric(survival_info$days_to_last_follow_up))
survival_info$Event <- ifelse(survival_info$vital_status == "Dead", 1, 0)

# 5.3: Merge with mir-202 data
survival_analysis_data <- merge(
  mirna_data, # (from step 4.3)
  survival_info,
  by.x = "patient_id",
  by.y = "submitter_id"
)
survival_analysis_data <- na.omit(survival_analysis_data)

# 5.4: Create Kaplan-Meier Plot
median_expr <- median(survival_analysis_data$mir_202_expr)
survival_analysis_data$mir_group <- ifelse(survival_analysis_data$mir_202_expr >= median_expr, "High_mir-202", "Low_mir-202")

fit <- survfit(Surv(Time, Event) ~ mir_group, data = survival_analysis_data)
km_plot <- ggsurvplot(
  fit, data = survival_analysis_data, pval = TRUE, risk.table = TRUE,
  legend.title = "Groups", legend.labs = c("High mir-202", "Low mir-202"),
  palette = c("red", "blue"), title = "Kaplan-Meier Plot for mir-202 in TCGA-ACC"
)

# 5.5: Save plot
pdf("ACC_survival_plot.pdf")
print(km_plot)
dev.off()

print(paste("Survival Plot p-value:", surv_pvalue(fit, data = survival_analysis_data)$pval))
print("--- SCRIPT 02 COMPLETE ---")