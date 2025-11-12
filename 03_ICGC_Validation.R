# SCRIPT 03: ICGC (PCAWG) Validation
# PURPOSE: Repeats the Selectivity Score analysis on an independent
#          pan-cancer cohort (ICGC/PCAWG) to validate the methodology.
# RUNTIME: ~5 minutes
# NOTE: This script requires MANUAL download of data files first.

# --- 1. LOAD LIBRARIES ---
library(dplyr)

# --- 2. MANUAL DOWNLOAD (User Action Required) ---
# Due to server-side restrictions (403 Forbidden errors),
# these files MUST be downloaded manually from a browser.
#
# INSTRUCTIONS:
# 1. Open this URL: https://pcawg.xenahubs.net/datapages/?dataset=x3t2m1.mature.TMM.mirna.matrix.donor.log&host=https%3A%2F%2Fpcawg.xenahubs.net
# 2. Click "download" and save "x3t2m1.mature.TMM.mirna.matrix.donor.log" to this project folder.
#
# 3. Open this URL: https://pcawg.xenahubs.net/datapages/?dataset=project_code_donor&host=https%3A%2F%2Fpcawg.xenahubs.net
# 4. Click "download" and save "project_code_donor" to this project folder.
#
# (We are not using the survival data for ICGC in this example)

# --- 3. LOAD ICGC MIRNA MATRIX ---
print("Loading ICGC miRNA matrix (manual download)...")
file_name_mirna <- "x3t2m1.mature.TMM.mirna.matrix.donor.log"
icgc_mirna_matrix_raw <- read.table(
  file_name_mirna, header = TRUE, sep = "\t", fill = TRUE
)
# Clean and set rownames
rownames(icgc_mirna_matrix_raw) <- icgc_mirna_matrix_raw[, 1]
icgc_mirna_matrix <- icgc_mirna_matrix_raw[, -1]
colnames(icgc_mirna_matrix) <- gsub("^X", "", colnames(icgc_mirna_matrix))

# --- 4. LOAD ICGC METADATA (CANCER CODES) ---
print("Loading ICGC metadata (manual download)...")
file_name_meta <- "project_code_donor"
icgc_project_data <- read.table(
  file_name_meta, header = TRUE, sep = "\t", fill = TRUE, quote = ""
)

# --- 5. ALIGN DATA ---
print("Aligning ICGC matrix and metadata...")
common_samples <- intersect(
  colnames(icgc_mirna_matrix),
  icgc_project_data$icgc_donor_id
)
icgc_matrix_final <- icgc_mirna_matrix[, common_samples]
icgc_meta_final <- icgc_project_data[icgc_project_data$icgc_donor_id %in% common_samples, ]
icgc_meta_final <- icgc_meta_final[match(colnames(icgc_matrix_final), icgc_meta_final$icgc_donor_id), ]

# --- 6. RUN SELECTIVITY ANALYSIS (on ICGC) ---
print("Calculating Selectivity Scores for ICGC data...")
cancer_types_icgc <- unique(icgc_meta_final$dcc_project_code)
mean_matrix_icgc <- matrix(NA, nrow = nrow(icgc_matrix_final), ncol = length(cancer_types_icgc))
rownames(mean_matrix_icgc) <- rownames(icgc_matrix_final)
colnames(mean_matrix_icgc) <- cancer_types_icgc

for (cancer in cancer_types_icgc) {
  samples_for_this_cancer <- icgc_meta_final$icgc_donor_id[icgc_meta_final$dcc_project_code == cancer]
  mean_values <- rowMeans(icgc_matrix_final[, samples_for_this_cancer, drop = FALSE], na.rm = TRUE)
  mean_matrix_icgc[, cancer] <- mean_values
}

# 6.2: Calculate Score
selectivity_scores_icgc <- data.frame(
  miRNA = rownames(mean_matrix_icgc), max_expr = NA, max_cancer = NA,
  second_max_expr = NA, selectivity_score = NA
)
for (i in 1:nrow(mean_matrix_icgc)) {
  expr_row <- mean_matrix_icgc[i, ]
  sorted_expr <- sort(expr_row, decreasing = TRUE, na.last = TRUE)
  max_val <- sorted_expr[1]
  second_max_val <- sorted_expr[2]
  max_cancer_name <- names(sorted_expr)[1]
  selectivity_scores_icgc$max_expr[i] <- max_val
  selectivity_scores_icgc$max_cancer[i] <- max_cancer_name
  selectivity_scores_icgc$second_max_expr[i] <- second_max_val
  selectivity_scores_icgc$selectivity_score[i] <- max_val - second_max_val
}

# 6.3: Create and save the final ranked list
final_selective_list_icgc <- selectivity_scores_icgc[order(selectivity_scores_icgc$selectivity_score, decreasing = TRUE), ]
print("Top 10 selective miRNAs (ICGC):")
print(head(final_selective_list_icgc, 10))
saveRDS(final_selective_list_icgc, "final_selective_list_ICGC.rds")

print("--- SCRIPT 03 COMPLETE ---")