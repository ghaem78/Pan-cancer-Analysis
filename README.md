# A Validated Pan-Cancer Atlas of Selective microRNAs
# Code for our article...

This repository contains the complete R code for the bioinformatics analysis described in our paper. The analysis is divided into two main projects (TCGA Discovery and ICGC Validation) and split into four sequential scripts.

## Analysis Workflow

The scripts are numbered and designed to be run in order:

### 1. TCGA Discovery Project

* **`01_TCGA_miRNA_Processing.R`**:
    * **Purpose:** To download, clean, and normalize the miRNA-Seq data for all 33 TCGA Pan-Cancer cohorts.
    * **Runtime:** Very long (approx. 30-45 minutes).
    * **Key Output:** `vsd_miRNA_Normalized.rds` (The main normalized TCGA data object).

* **`02_TCGA_Analysis_and_Discovery.R`**:
    * **Purpose:** Loads the normalized data from Script 01. Performs the primary "Selectivity Score" analysis to identify selective miRNAs (like `hsa-mir-202`). Then, performs the downstream case-study on `TCGA-ACC` (correlation with RNA-Seq and Survival Analysis).
    * **Runtime:** Long (approx. 5-10 minutes, includes RNA-Seq download).
    * **Key Outputs:** `final_selective_list_TCGA.rds`, `ACC_target_list_suppressed.rds`, `ACC_target_list_activated.rds`, `ACC_survival_plot.pdf`.

### 2. ICGC Validation Project

* **`03_ICGC_Validation.R`**:
    * **Purpose:** To perform an independent validation of our "Selectivity Score" methodology on the ICGC/PCAWG pan-cancer dataset.
    * **Runtime:** Fast (approx. 2-5 minutes, includes data download).
    * **Note:** Requires manual download of `pcawg.miRNA.rpm.gz` and metadata files as described in the script comments, due to server-side download restrictions.
    * **Key Output:** `final_selective_list_ICGC.rds`.

### 3. Final Comparison

* **`04_Final_Comparison.R`**:
    * **Purpose:** Loads the final result lists from both Script 02 and Script 03. Compares them to generate the final validation table (Table 1 in the manuscript).
    * **Runtime:** Very fast (< 1 minute).
    * **Key Output:** `Final_Validation_Table.csv`.

## Requirements

All required R packages (e.g., `TCGAbiolinks`, `DESeq2`, `survival`, `AnnotationDbi`) are listed at the top of each script.
