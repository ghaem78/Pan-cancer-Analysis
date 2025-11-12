# SCRIPT 04: Final TCGA vs ICGC Comparison
# PURPOSE: Loads the final results from both projects
#          and creates the cross-cohort validation table.
# RUNTIME: < 1 minute
# INPUTS: final_selective_list_TCGA.rds, final_selective_list_ICGC.rds

# --- 1. LOAD LIBRARIES ---
library(dplyr)
library(purrr)

# --- 2. LOAD FINAL RESULTS FROM BOTH PROJECTS ---
print("Loading final results from TCGA and ICGC...")
final_selective_list_tcga <- readRDS("final_selective_list_TCGA.rds")
final_selective_list_icgc <- readRDS("final_selective_list_ICGC.rds")

# --- 3. DEFINE CLEANING FUNCTION ---
clean_mir_name <- function(names) {
  names_clean <- tolower(names)               # 1. lowercase
  names_clean <- gsub("-\\d*p$", "", names_clean)  # 2. remove -3p/-5p
  names_clean <- gsub("-\\d$", "", names_clean)    # 3. remove -1/-2
  return(names_clean)
}

# --- 4. FIND COMMON CANCERS ---
tcga_codes <- unique(final_selective_list_tcga$max_cancer)
tcga_codes_short <- gsub("TCGA-", "", tcga_codes)
icgc_codes <- unique(final_selective_list_icgc$max_cancer)
icgc_codes_short <- gsub("-US$|-AU$|-UK$|-EU$|-JP$|-CA$", "", icgc_codes)
common_cancers <- intersect(tcga_codes_short, icgc_codes_short)

print(paste("Found", length(common_cancers), "common cancer types to compare."))

# --- 5. BUILD VALIDATION TABLE (like step 48) ---
print("Building final validation table...")
validation_results <- list()

for (cancer_code in common_cancers) {
  tcga_full_name <- tcga_codes[tcga_codes_short == cancer_code]
  icgc_full_name <- icgc_codes[icgc_codes_short == cancer_code]
  
  if (length(tcga_full_name) == 0 || length(icgc_full_name) == 0) next
  
  top_tcga <- final_selective_list_tcga[final_selective_list_tcga$max_cancer == tcga_full_name, ]
  top_5_tcga <- head(top_tcga$miRNA, 5)
  
  top_icgc <- final_selective_list_icgc[final_selective_list_icgc$max_cancer == icgc_full_name, ]
  top_5_icgc <- head(top_icgc$miRNA, 5)
  
  clean_tcga <- clean_mir_name(top_5_tcga)
  clean_icgc <- clean_mir_name(top_5_icgc)
  validated_mirnas <- intersect(clean_tcga, clean_icgc)
  
  validation_results[[cancer_code]] <- list(
    TCGA_Top_5 = top_5_tcga,
    ICGC_Top_5 = top_5_icgc,
    Validated_Families = validated_mirnas
  )
}

# --- 6. CONVERT LIST TO FINAL TABLE & SAVE ---
final_summary_table <- map_dfr(validation_results, ~data.frame(
  TCGA_Top_Hits = paste(.x$TCGA_Top_5, collapse = ", "),
  ICGC_Top_Hits = paste(.x$ICGC_Top_5, collapse = ", "),
  Validated_Families = paste(.x$Validated_Families, collapse = ", ")
), .id = "Cancer_Type")
final_summary_table$Validated_Families[final_summary_table$Validated_Families == ""] <- "None"

print("--- FINAL VALIDATION TABLE ---")
print(final_summary_table)
write.csv(final_summary_table, "Final_Validation_Table.csv", row.names = FALSE)

print("--- SCRIPT 04 COMPLETE ---")