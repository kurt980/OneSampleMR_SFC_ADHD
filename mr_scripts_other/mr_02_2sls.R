#

# Mendelian Randomization Analysis (Sensitivity Analysis Forward MR)

# 2. IV Regression (2-Stage Least Square)

################################################################################
## Readme: this code is for running IV regression using IVs selected by C+T
################################################################################


library(snpStats)
library(tibble)
library(dplyr)
library(AER)
library(readxl)

trait_name <- 'adhd'                    # trait name
pvalthrs <- c(5e-08)                    # thresholds 5e-06, 5e-07,                   
qc <- "_qc"                             # whether bfile used is QCed
unrel <- "_unrelated" # _unrelated"
method <- "gwasplink20_covsex_qc_unrelated_sd_harm2"                    # corresponds to 1009 GWAS results (gwasplink2.0 and fastGWA using sex as covariates)
r2 <- 0.001

### Load Utilities, data, tool and define directory
utilities_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_SCFC_utilities/"
bfile_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_bfile/bfile_for_gwas_comparison/"
cov_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Covariates/"
harm_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Harmonized_scfc/"
pheno_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Phenotype/"
tool_path <- "/Users/kurt/Documents/Data_Science_Project/Tools/"

# if (qc == "_qc") { bfile_path <- paste0(bfile_path,"scfc_all_qc", unrel) # read bfile
# } else { bfile_path <- paste0(bfile_path,"scfc_all_noqc", unrel) }
bfile_path <- paste0(bfile_path, "ABCD_all")
directory <- paste0("/Users/kurt/Documents/Data_Science_Project/1009_MR/", trait_name, "/", method, "_MR/")                     # directory for task
ct_path <- paste0(directory, "C+T/C+Ted_", method, "/")                                                                         # directory for C+T results
#mr_dir <- paste0(directory, "MR_results/")                                                                                      # directory for storing mr results
work_dir <- paste0(directory, "MR_results/")                                                                                # directory for storing current mr task results
#if (!dir.exists(mr_dir)) { dir.create(mr_dir, recursive = F) }
if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }

load(paste0(utilities_path, "coupling_cbcl.RData"))
load(paste0(utilities_path, "scfc_coupling_matched.RData"))

## This step is not entirely necessary, only aiming to proof read the id correspondence
#-------------------------------------------------------------------------------
new_fam <- read.table(paste0(utilities_path, "acspsw03.txt"), header=T) %>%       # get new fam id; this step is not necessary as the .fam files already uses new fam ids, just verifying
  filter(eventname=="baseline_year_1_arm_1") %>%
  select(rel_family_id, subjectkey) %>%
  rename(fid = rel_family_id, iid = subjectkey) %>%
  mutate(fid = as.integer(fid))
cbcl_fam <- read.csv(paste0(utilities_path, "demo_cbcl_baseline_n11876.csv")) %>% # get cbcl iid;
  select(subjectkey) %>%
  rename(iid = subjectkey)
old_bfile_fam <- read.table(paste0(paste0("/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_bfile/bfile_for_gwas_comparison/","ABCD_all_old.fam"))) %>% #substring(bfile_path, 1, (nchar(bfile_path) - nchar(paste0("scfc_all", qc, unrel)))), "ABCD_all_old.fam")) %>%           # get old bfile fam id
  select(V1, V2) %>%
  rename(fid_old = V1, iid = V2)
scfc_fam <- demo_cbcl_base_scfc %>%                                               # get scfc id
  select(subjectkey) %>%  
  mutate(subjectkey = sub("^(.{4})", "\\1_", subjectkey)) %>%
  rename(iid = subjectkey)
bfile_fam <- read.table(paste0(substring(bfile_path, 1, (nchar(bfile_path) - nchar("ABCD_all"))), "scfc_all", qc, unrel, ".fam")) %>%                           # get current bfile fam id
  select(V1, V2) %>%
  rename(fid = V1, iid = V2) %>%
  mutate(fid = as.integer(fid))

# get id correspondence
id_list <- new_fam %>%
  inner_join(cbcl_fam, by = "iid") %>%
  inner_join(old_bfile_fam, by = "iid") %>%
  inner_join(scfc_fam, by = "iid") %>%
  inner_join(bfile_fam, by = c("iid", "fid"))
#-------------------------------------------------------------------------------

## ------------------------read outcome and covariates--------------------------
df_outcome <- demo_cbcl_base_scfc %>% # Read outcome
  select(subjectkey, cbcl_scr_dsm5_depress_t, cbcl_scr_dsm5_anxdisord_t, cbcl_scr_dsm5_somaticpr_t, cbcl_scr_dsm5_adhd_t, cbcl_scr_dsm5_opposit_t, cbcl_scr_dsm5_conduct_t) %>%
  rename(iid = subjectkey, outcome = paste0("cbcl_scr_dsm5_",trait_name,"_t")) %>%
  mutate(iid = sub("^(.{4})", "\\1_", iid)) %>%
  inner_join(id_list, by = c("iid")) %>%
  select(fid, iid, fid_old, outcome)

CFDRs <- demo_cbcl_base_scfc %>% # get covariates
  select(subjectkey, site_id_l, interview_age, sex, race, ethn) %>%
  mutate(
    iid = sub("^(.{4})", "\\1_", subjectkey)  # Add '_' after the 4th character
    #ethn = as.factor(ethn),  # Convert ethnicity factor if needed
    #race = as.factor(race)   # Convert race factor if needed
  ) %>%
  inner_join(id_list, by = c("iid")) %>%
  select(fid, iid, site_id_l, interview_age, sex, race, ethn) %>%
  mutate(  # deal with spaces in variables
    race = gsub(" ", "_", race),                      # Replace spaces in 'race'
    ethn = gsub(" ", "_", ethn)                       # Replace spaces in 'ethnicity'
  )
# Fill NAs in race and ethnicity with mode
calculate_mode <- function(x) {
  ux <- unique(na.omit(x))  # Remove NA values before finding the mode
  ux[which.max(tabulate(match(x, ux)))]
}
CFDRs$race <- ifelse(is.na(CFDRs$race), calculate_mode(CFDRs$race), CFDRs$race)
CFDRs$ethn <- ifelse(is.na(CFDRs$ethn), calculate_mode(CFDRs$ethn), CFDRs$ethn)

df_out_cov <- left_join(df_outcome, CFDRs, by = c("fid", "iid")) # Merge with covariates
#-------------------------------------------------------------------------------



# Load necessary libraries for parallel computing
library(doParallel)
library(foreach)

# Set up parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

############################## Compute IVReg in Parallel ##############################
for (pvalthr in pvalthrs) {
  param <- paste0(method, "_p=", pvalthr) # Set parameters
  dir.create(paste0(work_dir, param, "/"), recursive = TRUE)
  
  # Create directory to save individual IV summaries
  iv_summary_dir <- paste0(work_dir, param, "/iv_summaries/")
  dir.create(iv_summary_dir, recursive = TRUE)
  
  # Parallel computing loop over ROI
  results <- foreach(roi = 1:379, .packages = c("dplyr", "tibble", "data.table", "snpStats", "AER")) %dopar% {
    tryCatch({
      #=========================== Read Exposure data ==============================
      exposure <- read.table(paste0(pheno_path, "scfc_sd_harm2_fam_roi", roi, ".phen"), header = FALSE)
      names(exposure)[1:3] <- c("fid", "iid", "exposure")
      
      df_out_cov_exp <- left_join(df_out_cov, exposure, by = c("fid", "iid"))
      
      #=========================== Read SNP data ===================================
      ct_pvalthr_path <- paste0(ct_path, "p=", pvalthr, "/")
      if (!file.exists(paste0(ct_pvalthr_path, 'roi', roi, '.', "_p=", pvalthr, "_r2=", r2, ".txt"))) {
        cat(paste0("\n\nNo SNP to use for roi ", roi, "\n"))
        return(NULL)
      } else {
        clumped_roi <- read.table(paste0(ct_pvalthr_path, 'roi', roi, '.', "_p=", pvalthr, "_r2=", r2, ".txt"), header = TRUE) %>% distinct()
        num_snps <- nrow(clumped_roi)
        cat(paste0("\n\n################### Doing ivreg for roi ", roi, "; number of SNPs to be used as IVs: ", num_snps, "; for pvalthr:", pvalthr, " ########################\n\n"))
      }
      
      #======================= Filter SNPs Using PLINK and Prepare Data =======================
      dir.create(paste0(work_dir, param, "/iv_bfile_c+t/"), recursive = FALSE)
      write.table(data.frame(clumped_roi$SNP), paste0(work_dir, param, "/iv_bfile_c+t/snp_id_roi", roi, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
      
      system(paste0(tool_path, "plink --bfile ", bfile_path,
                    " --extract ", paste0(work_dir, param, "/iv_bfile_c+t/snp_id_roi", roi, ".txt"),
                    " --make-bed --out ", work_dir, param, "/iv_bfile_c+t/iv_bfile_roi", roi),
             intern = FALSE)
      
      bfile <- read.plink(paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi", roi, ".bed"), 
                          paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi", roi, ".bim"), 
                          paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi", roi, ".fam"))
      sub_matrix <- bfile$genotypes
      SNPs <- as.data.frame(as.matrix(sub_matrix, type = "numeric"))
      SNPs <- tibble::rownames_to_column(SNPs, "iid")
      
      #======================== Merge all data =====================================
      df_all <- left_join(df_out_cov_exp, SNPs, by = "iid")
      snp_names <- names(df_all)[11:ncol(df_all)]
      df_all[snp_names] <- lapply(df_all[snp_names], function(x) as.numeric(as.character(x)) - 1)
      df_all <- df_all %>%
        filter(outcome != -9) %>%
        filter(!if_any(everything(), is.na))
      
      #=========================== Perform IV Regression ===========================
      # First, do OLS
      lm_fit <- lm(outcome ~ exposure + sex + interview_age + ethn + race, data = df_all)
      lm_fit_summary <- summary(lm_fit)
      ols_coef <- coef(lm_fit_summary)["exposure", "Estimate"]
      ols_pvalue <- coef(lm_fit_summary)["exposure", "Pr(>|t|)"]
      ols_n <- nobs(lm_fit)
      
      # Then, do IV regression
      Y <- df_all$outcome
      X <- model.matrix(~ exposure + sex + interview_age + ethn + race, data = df_all)
      IV_names <- names(df_all)[11:ncol(df_all)]
      Z <- cbind(model.matrix(~ sex + interview_age + ethn + race, data = df_all), 
                 as.matrix(df_all[, IV_names, drop = FALSE]))
      
      iv_fit <- ivreg(Y ~ X | Z)  # Using intercept by default
      iv_fit_summary <- summary(iv_fit)
      iv_coef <- coef(iv_fit_summary)["Xexposure", "Estimate"]
      iv_pvalue <- coef(iv_fit_summary)["Xexposure", "Pr(>|t|)"]
      iv_stderror <- coef(iv_fit_summary)["Xexposure", "Std. Error"]
      iv_n <- nobs(iv_fit)
      
      # Save full summary for inspection
      summary_file <- paste0(iv_summary_dir, "iv_summary_roi", roi, ".txt")
      capture.output(iv_fit_summary, file = summary_file)
      
      # Return results for this ROI
      data.frame(
        roi = roi,
        ROI = mmp_subcor[roi],
        OLS_Coefficient = ols_coef,
        OLS_P_Value = ols_pvalue,
        OLS_N = ols_n,
        IV_Coefficient = iv_coef,
        IV_StdError = iv_stderror,
        IV_P_Value = iv_pvalue,
        Num_snps = length(IV_names),
        IV_N = iv_n,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      return(NULL)  # Return NULL if any error occurs for this ROI
    })
  }
  
  results_df <- do.call(rbind, results)
  
  # Save results for this pvalthr value
  write.csv(results_df, paste0(work_dir, "iv_results_", param, "_intercept.csv"))
  cat(paste0("\n\n################### Successfully saved ivreg with intercept for pvalthr:", pvalthr, " ########################\n\n"))
}

# Stop the parallel backend
stopCluster(cl)








