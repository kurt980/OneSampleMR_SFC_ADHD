#

# Mendelian Randomization Analysis (Case/Control Analysis Forward MR)

# 2. IV Regression (2-Stage Residual Inclusion Estimation)

################################################################################
## Readme: this code is for running IV regression using IVs selected by C+T
################################################################################

library(snpStats)
library(tibble)
library(dplyr)
library(AER)
library(readxl)

trait_name <- 'ADHD'                    # trait name
pvalthrs <- c(5e-08)                    # thresholds                 
qc <- "_qc"                             # whether bfile used is QCed
gwas <- "gwasplink19"                   # method

### Load Utilities, data, tool and define directory
utilities_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_SCFC_utilities/"
bfile_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_bfile/"
cov_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Covariates/"
harm_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Harmonized_scfc/"
pheno_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Phenotype/"
tool_path <- "/Users/grc8mu/Desktop/DS/Tools/"

if (qc == "_qc") { bfile_path <- paste0(bfile_path,"scfc_all_qc") # read bfile
} else { bfile_path <- paste0(bfile_path,"scfc_all_noqc") }

directory <- "/Users/grc8mu/Desktop/DS/1102_MR/gwasplink19_qc_sd_harm2_MR/"                                        # directory for task
ct_path <- paste0(directory, "C+T/C+Ted/")                                                                         # directory for C+T results
mr_dir <- paste0("/Users/grc8mu/Desktop/DS/1102_MR/mr_one_sample_", trait_name, "/") # directory for storing mr results
work_dir <- paste0(mr_dir, "ivreg_gwasplink19/")                                                                   # directory for storing current mr task results
if (!dir.exists(mr_dir)) { dir.create(mr_dir, recursive = F) }
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
old_bfile_fam <- read.table(paste0(substring(bfile_path, 1, (nchar(bfile_path) - 11)), "ABCD_all_old.fam")) %>%           # get old bfile fam id
  select(V1, V2) %>%
  rename(fid_old = V1, iid = V2)
scfc_fam <- demo_cbcl_base_scfc %>%                                               # get scfc id
  select(subjectkey) %>%  
  mutate(subjectkey = sub("^(.{4})", "\\1_", subjectkey)) %>%
  rename(iid = subjectkey)
bfile_fam <- read.table(paste0(bfile_path, ".fam")) %>%                           # get current bfile fam id
  select(V1, V2) %>%
  rename(fid = V1, iid = V2)

# get id correspondence
id_list <- new_fam %>%
  inner_join(cbcl_fam, by = "iid") %>%
  inner_join(old_bfile_fam, by = "iid") %>%
  inner_join(scfc_fam, by = "iid") %>%
  inner_join(bfile_fam, by = c("iid", "fid"))
#-------------------------------------------------------------------------------

## ------------------------read outcome and covariates--------------------------
# Load the diagnosis dataset
dsm5_diag <- read_excel(paste0(utilities_path, "abcd_diagnosis_sum_baseline_n11976.xlsx"))
dsm5_diag <- dsm5_diag %>%
  select(subjectkey, ADHD, Dep, Bip, DMDD, Schi, Anx, OCD, Eat, ODD, Cond, PTSD) %>%
  rename(iid = subjectkey) %>%
  inner_join((id_list %>% select(fid, iid)), by = "iid") %>%
  select(fid, iid, ADHD, Dep, Bip, DMDD, Schi, Anx, OCD, Eat, ODD, Cond, PTSD)

df_outcome <- dsm5_diag %>% # Read outcome
  #select(subjectkey, cbcl_scr_dsm5_depress_t, cbcl_scr_dsm5_anxdisord_t, cbcl_scr_dsm5_somaticpr_t, cbcl_scr_dsm5_adhd_t, cbcl_scr_dsm5_opposit_t, cbcl_scr_dsm5_conduct_t) %>%
  rename(outcome = paste0("",trait_name,"")) %>%
  #mutate(iid = sub("^(.{4})", "\\1_", iid)) %>%
  #inner_join(id_list, by = c("iid")) %>%
  select(fid, iid, outcome)

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


############################## Compute 2SRI for Binary Outcome #####################################
for (pvalthr in pvalthrs) {
  param <- paste0(gwas, qc, "_p=", pvalthr) # Set parameters
  dir.create(paste0(work_dir, param, "/"), recursive = TRUE)
  
  # Initialize DataFrame to store results
  results_df <- data.frame(
    roi = integer(),
    ROI = character(),
    IV_Coefficient = numeric(),
    IV_StdError = numeric(),
    IV_P_Value = numeric(),
    Num_snps = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (roi in 1:1) {
    #===========================Read Exposure data==============================
    exposure <- read.table(paste0(pheno_path, "scfc_sd_harm2_fam_roi", roi, ".phen"), header = F)
    names(exposure)[1:3] <- c("fid", "iid", "exposure")
    
    df_out_cov_exp <- left_join(df_out_cov, exposure, by = c("fid", "iid"))
    
    #===========================Read SNP data===================================
    ct_pvalthr_path <- paste0(ct_path, "p=", pvalthr, "/") 
    if (!file.exists(paste0(ct_pvalthr_path, 'roi', roi, '.', "_p=", pvalthr, "_r2=", r2, ".txt"))) {
      print(paste0("\n\nNo SNP to use for roi ", roi, "\n"))
      next
    } else {
      clumped_roi <- read.table(paste0(ct_pvalthr_path, 'roi', roi, '.', "_p=", pvalthr, "_r2=", r2, ".txt"), header = T) %>% distinct()
      num_snps <- nrow(clumped_roi)
      cat(paste0("\n\n###################Performing 2SRI for roi ", roi, "; number of SNPs: ", num_snps, "; for pvalthr:", pvalthr, "########################\n\n"))
    }
    
    # Filter SNPs in PLINK
    dir.create(paste0(work_dir, param, "/iv_bfile_c+t/"), recursive = FALSE)
    write.table(data.frame(clumped_roi$SNP), paste0(work_dir, param, "/iv_bfile_c+t/snp_id_roi", roi, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    
    system(paste0(tool_path, "plink --bfile ", bfile_path,
                  " --extract ", paste0(work_dir, param, "/iv_bfile_c+t/snp_id_roi", roi, ".txt"),
                  " --make-bed --out ", work_dir, param, "/iv_bfile_c+t/iv_bfile_roi", roi),
           intern = FALSE)
    
    # Read bfiles and get them ready
    bfile <- read.plink(paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi", roi, ".bed"), 
                        paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi", roi, ".bim"), 
                        paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi", roi, ".fam"))
    sub_matrix <- bfile$genotypes
    SNPs <- as.data.frame(as.matrix(sub_matrix, type = "numeric"))
    SNPs <- tibble::rownames_to_column(SNPs, "iid")
    
    #========================Merge all data=====================================
    df_all <- left_join(df_out_cov_exp, SNPs, by = "iid")
    snp_names <- names(df_all)[10:ncol(df_all)]
    df_all[snp_names] <- lapply(df_all[snp_names], function(x) {
      as.numeric(as.character(x)) - 1
    })
    df_all <- df_all %>%
      filter(outcome %in% c(0, 1)) %>% # Only keep binary outcome values (0 or 1)
      filter(!if_any(everything(), is.na))
    
    #===========================2SRI Method===========================
    # Stage 1: Linear regression of exposure on SNPs and covariates to obtain residuals
    stage1 <- lm(exposure ~ . - outcome, data = df_all)  # Include SNPs and covariates
    df_all$residuals_stage1 <- residuals(stage1)
    
    # Stage 2: Logistic regression on outcome with exposure, residuals, and covariates
    stage2 <- glm(outcome ~ exposure + residuals_stage1 + sex + interview_age + ethn + race, 
                  family = binomial(link = "logit"), data = df_all)
    
    # Extract coefficient, standard error, and p-value for the exposure variable
    stage2_summary <- summary(stage2)
    iv_coef <- stage2_summary$coefficients["exposure", "Estimate"]
    iv_stderror <- stage2_summary$coefficients["exposure", "Std. Error"]
    iv_pvalue <- stage2_summary$coefficients["exposure", "Pr(>|z|)"]
    
    # Append results to the DataFrame
    results_df <- rbind(results_df, data.frame(
      roi = roi,
      ROI = mmp_subcor[roi],
      IV_Coefficient = iv_coef,
      IV_StdError = iv_stderror,
      IV_P_Value = iv_pvalue,
      Num_snps = num_snps,
      stringsAsFactors = FALSE
    ))
  }
  
  # Save results
  write.csv(results_df, paste0(work_dir, "2SRI_results_", param, ".csv"))
  cat(paste0("\n\n###################Successfully saved 2SRI results for roi ", roi, "; number of SNPs: ", num_snps, "; for pvalthr:", pvalthr, "########################\n\n"))
}
