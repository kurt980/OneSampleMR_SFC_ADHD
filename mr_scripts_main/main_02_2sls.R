#

# Mendelian Randomization Analysis (Main Analysis Forward MR)

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
pvalthrs <- c(5e-06, 5e-07, 5e-08)      # thresholds                 
qc <- "_qc"                             # whether bfile used is QCed
gwas <- "gwasplink19"                   # corresponds to 0928 GWAS using PLINK1.9 with fisher transformed SCFC and without covariates (redo 0600 GWAS)
r2 <-0.001

### Load Utilities, data, tool and define directory
utilities_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_SCFC_utilities/"
bfile_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_bfile/"
cov_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Covariates/"
harm_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Harmonized_scfc/"
pheno_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Phenotype/"
tool_path <- "/Users/kurt/Documents/Data_Science_Project/Tools/"

if (qc == "_qc") { bfile_path <- paste0(bfile_path,"scfc_all_qc") # read bfile
} else { bfile_path <- paste0(bfile_path,"scfc_all_noqc") }

directory <- "/Users/kurt/Documents/Data_Science_Project/1004_MR_GWASplink19/"                                     # directory for task
ct_path <- paste0(directory, "C+T/C+Ted/")                                                                         # directory for C+T results
mr_dir <- paste0("/Users/kurt/Documents/Data_Science_Project/1004_MR_GWASplink19/mr_one_sample_", trait_name, "/") # directory for storing mr results
work_dir <- paste0(mr_dir, "ivreg_gwasplink19/")                                                                   # directory for storing current mr task results
if (!dir.exists(mr_dir)) { dir.create(mr_dir, recursive = F) }
if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }

load(paste0(utilities_path, "coupling_cbcl.RData"))
load(paste0(utilities_path, "scfc_coupling_matched.RData"))

### This step is not entirely necessary, only aiming to proof read the id correspondence
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

### ------------------------read outcome and covariates-------------------------
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

##############################Compute IVReg#####################################
for (pvalthr in pvalthrs) {
  param <- paste0(gwas, qc, "_p=", pvalthr) #=============set parameters===========
  dir.create(paste0(work_dir, param, "/"), recursive = T)
  # Initialize DataFrame to store results
  results_df <- data.frame(
    roi = integer(),
    ROI = character(),
    OLS_Coefficient = numeric(),
    OLS_P_Value = numeric(),
    OLS_N = numeric(),
    IV_Coefficient = numeric(),
    IV_StdError = numeric(),
    IV_P_Value = numeric(),
    IV_N = numeric(),
    Num_snps = numeric(),
    stringsAsFactors = FALSE
  )
  for (roi in 1:379) {
    #===========================Read Exposure data==============================
    exposure <- read.table(paste0(pheno_path, "scfc_sd_harm2_fam_roi",roi,".phen"), header = F) # read scfc ROI
    names(exposure)[1] <- "fid"
    names(exposure)[2] <- "iid"
    names(exposure)[3] <- "exposure"
    
    df_out_cov_exp <- left_join(df_out_cov, exposure, by = c("fid","iid")) # join with outcome and covariates
    
    #===========================Read SNP data===================================
    ## Need to read bfile for SNP information after clumping+thresholding
    #### Since we already done the clumping, we just fetch our ids and filter the binary files:
    ct_pvalthr_path <- paste0(ct_path, "p=", pvalthr, "/") 
    if (!file.exists(paste0(ct_pvalthr_path, 'roi', roi, '.', "_p=", pvalthr, "_r2=", r2, ".txt"))) {
      print(paste0("\n\nNo SNP to use for roi ", roi, "\n"))
      next
    } else {
      clumped_roi <- read.table(paste0(ct_pvalthr_path, 'roi', roi, '.', "_p=", pvalthr, "_r2=", r2, ".txt"), header = T) %>% distinct() # in case there is duplicated values
      num_snps <- nrow(clumped_roi)
      cat(paste0("\n\n###################Doing ivreg for roi ", roi, "; number of SNPs to be used as IVs: ", num_snps, "; for pvalthr:", pvalthr, "########################\n\n"))
    }
    #### The lists of IDs are in the C+T folder
    dir.create(paste0(work_dir, param, "/iv_bfile_c+t/"), recursive = F)
    write.table(data.frame(clumped_roi$SNP), paste0(work_dir, param, "/iv_bfile_c+t/snp_id_roi",roi,".txt"), row.names = F, quote = FALSE, sep = "\t")
    
    #### Filter the bfile SNPs using PLINK
    snp_id_list <- paste0(work_dir, param, "/iv_bfile_c+t/snp_id_roi",roi,".txt")
    system(paste0(tool_path, "plink --bfile ", bfile_path,
                  " --extract ", snp_id_list,
                  " --make-bed --out ", work_dir, param, "/iv_bfile_c+t/iv_bfile_roi",roi),
           intern = F)
    
    #### Read bfiles and get them ready
    bfile <- read.plink(paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi",roi,".bed"), 
                        paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi",roi,".bim"), 
                        paste0(work_dir, param, "/iv_bfile_c+t/iv_bfile_roi",roi,".fam"))
    sub_matrix <- bfile$genotypes
    SNPs <- as.data.frame(as.matrix(sub_matrix, type = "numeric"))
    SNPs <- tibble::rownames_to_column(SNPs, "ID")
    names(SNPs)[1] <- "iid"
    
    #========================Merge all data=====================================
    df_all <- left_join(df_out_cov_exp, SNPs, by = "iid")
    #### encode alleles as 1,2,3
    snp_names <- names(df_all)[11:ncol(df_all)] # REMEMBER TO CHANGE THE VOLUMN INDEX IF DATA STRUCTURE CHANGES!!! (NOT NECESSARILY 11)
    df_all[snp_names] <- lapply(df_all[snp_names], function(x) {
      as.numeric(as.character(x)) - 1
    })
    #### remove NAs in outcome and any observations **把ADHD=-9的去掉**
    df_all <- df_all %>%
      filter(outcome != -9) %>%
      filter(!if_any(everything(), is.na))
    
    ############################################################################
    
    #===========================Perform IV Regression===========================
    #############
    #----------------------------First do OLS-----------------------------------
    lm_fit <- lm(outcome ~ exposure + sex + interview_age + ethn + race, data = df_all)
    lm_fit_summary <- summary(lm_fit)
    ols_coef <- coef(lm_fit_summary)["exposure", "Estimate"]
    ols_pvalue <- coef(lm_fit_summary)["exposure", "Pr(>|t|)"]
    ols_n <- nobs(lm_fit)
    
    ##---------------------------Do 2SLS using ivreg----------------------------
    Y <- df_all$outcome
    X <- model.matrix(~ exposure + sex + interview_age + ethn + race, data = df_all)
    #### get SNP names
    IV_names <- names(df_all)[11:ncol(df_all)]
    #### formulate Z matrix
    Z <- cbind(model.matrix(~ sex + interview_age + ethn + race - 1, data = df_all),
               as.matrix(df_all[, IV_names, drop = FALSE]))
    
    ## Fit the IV model using ivreg with matrices
    iv_fit <- ivreg(Y ~ X-1 | Z-1)  # '-1' to exclude the intercept, as model.matrix includes it by default
    iv_fit_summary <- summary(iv_fit)
    iv_coef <- coef(iv_fit_summary)["Xexposure", "Estimate"]
    iv_pvalue <- coef(iv_fit_summary)["Xexposure", "Pr(>|t|)"]
    iv_stderror <- coef(iv_fit_summary)["Xexposure", "Std. Error"]
    iv_n <- nobs(iv_fit)
    
    #### Append results to the DataFrame
    results_df <- rbind(results_df, data.frame(
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
    ))
    
  }
  ## save
  results_df <- results_df %>%
    mutate(roi = paste0("roi", roi))
  write.csv(results_df, paste0(work_dir, "iv_results_",param,".csv"))
  cat(paste0("\n\n###################Successfully saved ivreg for roi ", roi, "; number of SNPs to be used as IVs: ", num_snps, "; for pvalthr:", pvalthr, "########################\n\n"))
}