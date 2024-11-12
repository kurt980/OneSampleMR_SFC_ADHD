#

# Sensitivity Analysis

# 1. Association Analysis (Pearson's, Kendall's Tau, Spearman, Distance Correlation, Linear Regression)

## Read Data
library(snpStats)
library(tibble)
library(dplyr)
library(AER)
library(readxl)
library(ggplot2)
library(energy)     # For distance correlation
library(infotheo)   # For mutual information
library(Hmisc)      # For Hoeffding's D
library(minerva)    # For MIC

### Load Utilities, data, tool and define directory
utilities_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_SCFC_utilities/"
bfile_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_bfile/"
cov_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Covariates/"
harm_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Harmonized_scfc/"
pheno_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Phenotype/"
tool_path <- "/Users/kurt/Documents/Data_Science_Project/Tools/"
directory <- "/Users/kurt/Documents/Data_Science_Project/1003_KendallTau/"
load(paste0(utilities_path, "coupling_cbcl.RData"))
load(paste0(utilities_path, "scfc_coupling_matched.RData"))

### Compute association
#### Load dsm5 (ADHD and others) scores
dsm5 <- demo_cbcl_base_scfc %>%
  select(subjectkey, cbcl_scr_dsm5_depress_t, cbcl_scr_dsm5_anxdisord_t, cbcl_scr_dsm5_somaticpr_t, cbcl_scr_dsm5_adhd_t, cbcl_scr_dsm5_opposit_t, cbcl_scr_dsm5_conduct_t) %>%
  rename(iid = subjectkey) %>%
  mutate(iid = sub("^(.{4})", "\\1_", iid))
#### Load SCFC ROIs
csv_files <- list.files(path = paste0(harm_path, ""), pattern = "^scfc.*\\.csv$", full.names = FALSE)
combined_results <- list()
#### Compute correlations iteratively
for (scfc_csv in csv_files) {
  scfc <- read.csv(paste0(harm_path, "", scfc_csv))
  df <- dsm5 %>%
    select(iid, cbcl_scr_dsm5_adhd_t) %>%
    rename(adhd = cbcl_scr_dsm5_adhd_t) %>%
    inner_join(scfc, by = "iid")
  
  adhd <- as.vector(df$adhd)
  ROIs <- as.matrix(df[, 4:382])  # starting from column 4
  # create ROI number corresponding to name
  roi_numbers <- 1:379
  roi_df <- data.frame(
    roi = roi_numbers,      # The numeric sequence
    ROI = colnames(ROIs)        # The list of ROI names
  )
  #=============================================================================
  # Pearson's Correlation
  pearson_results <- apply(ROIs, 2, function(x) {
    test <- cor.test(x, adhd, use = "complete.obs")  # Perform correlation test
    c(correlation = test$estimate, p_value = test$p.value)  # Extract correlation and p-value
  })
  pearson_df <- data.frame(
    ROI = colnames(ROIs),
    coeff_pearson = pearson_results["correlation.cor", ],
    p_value_pearson = pearson_results["p_value", ]
  ) %>% inner_join(roi_df, by = "ROI") %>% select(roi, ROI, coeff_pearson, p_value_pearson)
  if (!dir.exists(paste0(directory, "pearson/"))) {
    dir.create(paste0(directory, "pearson/"), recursive = F) }
  write.csv(pearson_df, paste0(directory, "pearson/adhd_pearson_results_", sub("^scfc_(.*)\\.csv$", "\\1", basename(scfc_csv)), ".csv"), row.names = F)
  #=============================================================================
  # Linear Regression
  model <- lm(adhd ~ ROIs)  # 没用X因为有bug
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients  # This returns a matrix
  reg_df <- data.frame(
    ROI = sub("ROIs", "", rownames(coefficients)),  # Variable names
    coeff_reg = coefficients[, "Estimate"],  # Coefficients
    p_value_reg = coefficients[, "Pr(>|t|)"]
  ) %>% inner_join(roi_df, by = "ROI") %>% select(roi, ROI, coeff_reg, p_value_reg)
  if (!dir.exists(paste0(directory, "linear_reg/"))) {
    dir.create(paste0(directory, "linear_reg/"), recursive = F) }
  write.csv(reg_df, paste0(directory, "linear_reg/adhd_linear_reg_results_", sub("^scfc_(.*)\\.csv$", "\\1", basename(scfc_csv)), ".csv"), row.names = F)
  #=============================================================================
  # Kendall's Tau
  kendall_results <- apply(ROIs, 2, function(x) {
    test <- cor.test(x, adhd, method = "kendall", use = "complete.obs")  # Perform Kendall's tau test
    c(tau = test$estimate, p_value = test$p.value)  # Extract Kendall's tau and p-value
  })
  kendall_df <- data.frame(
    ROI = colnames(ROIs),
    coeff_kendall = kendall_results["tau.tau", ],  # Kendall's tau coefficients
    p_value_kendall = kendall_results["p_value", ]
  ) %>% inner_join(roi_df, by = "ROI") %>% select(roi, ROI, coeff_kendall, p_value_kendall)
  if (!dir.exists(paste0(directory, "kendall/"))) {
    dir.create(paste0(directory, "kendall/"), recursive = F) }
  write.csv(kendall_df, paste0(directory, "kendall/adhd_kendall_results_", sub("^scfc_(.*)\\.csv$", "\\1", basename(scfc_csv)), ".csv"), row.names = F)
  #=============================================================================
  # Spearman's Rank Correlation
  spearman_results <- apply(ROIs, 2, function(x) {
    test <- cor.test(x, adhd, method = "spearman", use = "complete.obs", exact = FALSE)  # Perform Spearman's rank correlation
    c(spearman = test$estimate, p_value = test$p.value)  # Extract Spearman's rho and p-value
  })
  spearman_df <- data.frame(
    ROI = colnames(ROIs),
    coeff_spearman = spearman_results["spearman.rho", ],  # Spearman's rho coefficients
    p_value_spearman = spearman_results["p_value", ]
  ) %>% inner_join(roi_df, by = "ROI") %>% select(roi, ROI, coeff_spearman, p_value_spearman)
  if (!dir.exists(paste0(directory, "spearman/"))) {
    dir.create(paste0(directory, "spearman/"), recursive = F) }
  write.csv(spearman_df, paste0(directory, "spearman/adhd_spearman_results_", sub("^scfc_(.*)\\.csv$", "\\1", basename(scfc_csv)), ".csv"), row.names = F)
  #=============================================================================
  # Merge all results into one data frame by ROI
  combined_df <- Reduce(function(x, y) merge(x, y, by = c("ROI", "roi"), all = TRUE),
                        list(pearson_df, reg_df, kendall_df, spearman_df)) %>% arrange(roi)
  
  combined_results[[scfc_csv]] <- combined_df %>% arrange(roi) # Add the combined results for this file to the list
  write.csv(combined_df, paste0(directory, "adhd_combined_association_results_", sub("^scfc_(.*)\\.csv$", "\\1", basename(scfc_csv)), ".csv"), row.names = F)
}

# final_combined_df <- do.call(rbind, combined_results) # Combine all results from all files into a single data frame
# final_combined_df <- final_combined_df[ , !duplicated(names(final_combined_df))] %>% arrange(roi)
# write.csv(final_combined_df, paste0(directory, "all_adhd_combined_association_results.csv"), row.names = F) # save the final combined results to a single CSV file
