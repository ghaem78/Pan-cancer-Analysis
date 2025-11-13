# A Validated Pan-Cancer Atlas of Selective microRNAs
# Code for: [The Selectivity Map: A Cross-Cohort Atlas Unmasks...]

This repository contains the complete R code for the bioinformatics analysis. The analysis is divided into two main projects (TCGA Discovery and ICGC Validation) and split into four sequential scripts.

## Analysis Workflow

The scripts are numbered and designed to be run in order:

### 1. TCGA Discovery Project

* **`01_TCGA_miRNA_Processing.R`**:
    * **Purpose:** Downloads, cleans, and normalizes the miRNA-Seq data for all 33 TCGA Pan-Cancer cohorts.
    * **Runtime:** Very long (approx. 30-45 minutes).
    * **Key Output:** `vsd_miRNA_Normalized.rds` (The main normalized TCGA data object).

* **`02_TCGA_Analysis_and_CaseStudies.R`**:
    * **Purpose:** Loads the normalized data from Script 01.
    * 1. Runs the primary "Selectivity Score" analysis (saves `final_selective_list_TCGA.rds`).
    * 2. Runs the **ACC Case Study** (downloads `TCGA-ACC` RNA-Seq, correlates with `hsa-mir-202`, runs survival analysis).
    * 3. Runs the **TGCT Case Study** (downloads `TCGA-TGCT` RNA-Seq, correlates with `hsa-mir-302b`).
    * **Runtime:** Long (approx. 15-25 minutes, includes two RNA-Seq downloads).
    * **Key Outputs:** `final_selective_list_TCGA.rds`, `ACC_target_lists.rds`, `TGCT_target_lists.rds`, `ACC_survival_plot.pdf`.

### 2. ICGC Validation Project

* **`03_ICGC_Validation.R`**:
    * **Purpose:** Repeats the "Selectivity Score" analysis on the independent ICGC/PCAWG cohort to validate the methodology.
    * **Runtime:** Fast (approx. 2-5 minutes).
    * **Note:** Requires *manual* download of data files (see script comments).
    * **Key Output:** `final_selective_list_ICGC.rds`.

### 3. Final Comparison

* **`04_Final_Comparison.R`**:
    * **Purpose:** Loads the final lists from Script 02 and 03. Compares them to generate the final validation table (Table 1 in the manuscript).
    * **Runtime:** Very fast (< 1 minute).
    * **Key Output:** `Final_Validation_Table.csv`.
