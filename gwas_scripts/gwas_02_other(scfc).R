#

# Running GWAS 

# 2. Script for running GWAS for all other methods

###############################################################################################
## Readme: this code runs and saves GWAS in different settings for SCFC
## The GWAS used for sensitivity analysis is included (sex as covariates, on unrelated samples)
###############################################################################################


## I. create covariates, filter out unrelated observations and create bfiles 

library(snpStats)
library(tibble)
library(dplyr)

# Sample some ROIs or just select all 379
set.seed(2024)  # Optional: to make the result reproducible
rois <- sample(1:379, 10)
rois <- 1:379
rois <- c(1,10,31,100,301,305,308,309)
print(rois)

### Load Utilities, data, tool and define directory
utilities_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_SCFC_utilities/"
bfile_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/ABCD_bfile/bfile_for_gwas_comparison/"
cov_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Covariates/cov_for_gwas_comparison/"
harm_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Harmonized_scfc/"
pheno_path <- "/Users/kurt/Documents/Data_Science_Project/ABCD_data/Phenotype/"
tool_path <- "/Users/kurt/Documents/Data_Science_Project/Tools/"

directory <- "/Users/kurt/Documents/Data_Science_Project/1006_GWAS_comparison/"

load(paste0(utilities_path, "coupling_cbcl.RData"))
load(paste0(utilities_path, "scfc_coupling_matched.RData"))


# ########################Create bfile, covariates################################
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
old_bfile_fam <- read.table(paste0(bfile_path, "ABCD_all_old.fam")) %>%           # get old bfile fam id
  select(V1, V2) %>%
  rename(fid_old = V1, iid = V2)
scfc_fam <- demo_cbcl_base_scfc %>%                                               # get scfc id
  select(subjectkey) %>%
  mutate(subjectkey = sub("^(.{4})", "\\1_", subjectkey)) %>%
  rename(iid = subjectkey)

# get id correspondence
id_list <- new_fam %>%
  inner_join(cbcl_fam, by = "iid") %>%
  inner_join(old_bfile_fam, by = "iid") %>%
  inner_join(scfc_fam, by = "iid")
#-------------------------------------------------------------------------------

##-------------------------read outcome and covariates--------------------------
# df_outcome <- demo_cbcl_base_scfc %>% # Read outcome
#   select(subjectkey, cbcl_scr_dsm5_depress_t, cbcl_scr_dsm5_anxdisord_t, cbcl_scr_dsm5_somaticpr_t, cbcl_scr_dsm5_adhd_t, cbcl_scr_dsm5_opposit_t, cbcl_scr_dsm5_conduct_t) %>%
#   rename(iid = subjectkey, outcome = paste0("cbcl_scr_dsm5_",trait_name,"_t")) %>%
#   mutate(iid = sub("^(.{4})", "\\1_", iid)) %>%
#   inner_join(id_list, by = c("iid")) %>%
#   select(fid, iid, fid_old, outcome)

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
#-------------------------------------------------------------------------------
#-----------------------Create ID files-----------------------------------------
# Create id file for all observations
fid_iid_scfc <- CFDRs %>% # write new id list to local
  select(fid, iid) %>%
  arrange(iid)
write.table(fid_iid_scfc, paste0(bfile_path, "fid_iid_scfc.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
# Create id file for unrelated observations
set.seed(2024)
fid_iid_scfc_unrelated <- fid_iid_scfc %>%
  group_by(fid) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  arrange(iid)

write.table(fid_iid_scfc_unrelated, paste0(bfile_path, "fid_iid_scfc_unrelated.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
# inspect, should be 5628(5487 after removing covar, but not here)
cat(length(unique(fid_iid_scfc$fid)), length(unique(read.table(paste0(bfile_path, "fid_iid_scfc_unrelated.txt"))[,1])))


##======================Create covariate files================================== **Technically covariates do not need to be subset as the software automatically picks the ids in the .fam file**
##---------------------------All obs (6381)-------------------------------------
qcovar_age <- CFDRs %>%
  select(fid, iid, interview_age) %>%
  arrange(iid)
covar_sex <- CFDRs %>%
  select(fid, iid, sex) %>%
  arrange(iid)
covar_sex_race_ethn <- CFDRs %>%
  select(fid, iid, sex, race, ethn) %>%
  arrange(iid)
covar_all <- qcovar_age %>%
  inner_join(covar_sex_race_ethn, by = c("fid", "iid")) %>%
  arrange(iid)

write.table(qcovar_age, paste0(cov_path, "qcovar_age.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(covar_sex, paste0(cov_path, "covar_sex.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(covar_sex_race_ethn, paste0(cov_path, "covar_sex_race_ethn.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
write.table(covar_all, paste0(cov_path, "covar_all.txt"), col.names = F, row.names = F, quote = F, sep = "\t")


##------------------------Unrelated obs (5628)---------------------------------- **Technically covariates do not need to be subset for unrelated obs as the software automatically picks the ids in the .fam file**
# qcovar_age_unrelated <- qcovar_age %>%
#   inner_join(fid_iid_scfc_unrelated, by = c("fid", "iid"))
# covar_sex_race_ethn_unrelated <- covar_sex_race_ethn %>%
#   inner_join(fid_iid_scfc_unrelated, by = c("fid", "iid"))
# covar_all_unrelated <- covar_all %>%
#   inner_join(fid_iid_scfc_unrelated, by = c("fid", "iid"))
#
# write.table(qcovar_age_unrelated, paste0(cov_path,"qcovar_age_unrelated.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(covar_sex_race_ethn_unrelated, paste0(cov_path,"covar_sex_race_ethn_unrelated.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
# write.table(covar_all_unrelated, paste0(cov_path,"covar_all_unrelated.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
#

##===========================Create bfiles======================================
#----------------------------Create new fam-------------------------------------
old_fam <- read.table(paste0(bfile_path, "ABCD_all_old.fam"))
colnames(old_fam) <- c("fid_old", "iid", "pid", "mid", "sex", "phenotype")

new_bfile_fam <- old_fam %>%
  left_join(new_fam, by = "iid") %>%
  mutate(fid = ifelse(is.na(fid), fid_old, fid)) %>%
  select(fid, iid, pid, mid, sex, phenotype)

write.table(new_bfile_fam, paste0(bfile_path, "ABCD_all.fam"), row.names = F, col.names = F, quote = F, sep = "\t")

#---------------------Create QC & No QC bfiles (All obs)------------------------
# Subset by id
system(paste0(tool_path, "plink --bfile ", bfile_path, "ABCD_all --keep ", bfile_path, "fid_iid_scfc.txt --make-bed --out ", bfile_path, "scfc_all_noqc"), intern = FALSE)
# QC
system(paste0(tool_path, "plink --bfile ", bfile_path, "scfc_all_noqc --geno 0.1 --maf 0.05 --mind 0.1 --make-bed --out ", bfile_path, "scfc_all_qc"), intern = FALSE)
print("############################Made noQC and QC bfile for all obs#######################################################################")

#---------------------Create QC & No QC bfiles (Unrelated)----------------------
# Subset by id
system(paste0(tool_path, "plink --bfile ", bfile_path, "ABCD_all --keep ", bfile_path, "fid_iid_scfc_unrelated.txt --make-bed --out ", bfile_path, "scfc_all_noqc_unrelated"), intern = FALSE)
# QC
system(paste0(tool_path, "plink --bfile ", bfile_path, "scfc_all_noqc_unrelated --geno 0.1 --maf 0.05 --mind 0.1 --make-bed --out ", bfile_path, "scfc_all_qc_unrelated"), intern = FALSE)
print("############################Made noQC and QC bfile for unrelated obs#######################################################################")


##===============Create GRM matrices and Sparse GRM for GCTA====================
#---------------------------Create GRM matrix for All obs-----------------------
# No QC
lib_path <- 'export DYLD_LIBRARY_PATH="/Users/kurt/Documents/Data_Science_Project/gcta-1.94.2-MacOS-ARM-x86_64/gcclib:$DYLD_LIBRARY_PATH"'
system(lib_path, intern = TRUE)  # Set library path
system(
  paste0(tool_path, "gcta64 ", "--bfile ", bfile_path, "scfc_all_noqc ",
      "--autosome  --make-grm  ", "--thread-num 16 ",
      "--out ", bfile_path, "scfc_all_noqc_out"), intern = F)
system(
  paste0(tool_path, "gcta64 ", "--grm ", bfile_path, "scfc_all_noqc_out ",
      "--make-bK-sparse 0.05 ", "--thread-num 16 ",
      "--out ", bfile_path, "scfc_all_noqc_sp_out"), intern = F)
# QC
lib_path <- 'export DYLD_LIBRARY_PATH="/Users/kurt/Documents/Data_Science_Project/gcta-1.94.2-MacOS-ARM-x86_64/gcclib:$DYLD_LIBRARY_PATH"'
system(lib_path, intern = TRUE)  # Set library path
system(
  paste0(tool_path, "gcta64 ", "--bfile ", bfile_path, "scfc_all_qc ",
      "--autosome  --make-grm  ", "--thread-num 16 ",
      "--out ", bfile_path, "scfc_all_qc_out"), intern = F)
system(
  paste0(tool_path, "gcta64 ", "--grm ", bfile_path, "scfc_all_qc_out ",
      "--make-bK-sparse 0.05 ", "--thread-num 16 ",
      "--out ", bfile_path, "scfc_all_qc_sp_out"), intern = F)

#---------------------------Create GRM matrix for unrelated obs-----------------
# No QC
lib_path <- 'export DYLD_LIBRARY_PATH="/Users/kurt/Documents/Data_Science_Project/gcta-1.94.2-MacOS-ARM-x86_64/gcclib:$DYLD_LIBRARY_PATH"'
system(lib_path, intern = TRUE)  # Set library path
system(
  paste0(tool_path, "gcta64 ", "--bfile ", bfile_path, "scfc_all_noqc_unrelated ",
      "--autosome  --make-grm  ", "--thread-num 16 ",
      "--out ", bfile_path, "scfc_all_noqc_unrelated_out"), intern = F)
system(
  paste0(tool_path, "gcta64 ", "--grm ", bfile_path, "scfc_all_noqc_unrelated_out ",
      "--make-bK-sparse 0.05 ", "--thread-num 16 ",
      "--out ", bfile_path, "scfc_all_noqc_unrelated_sp_out"), intern = F)
# QC
lib_path <- 'export DYLD_LIBRARY_PATH="/Users/kurt/Documents/Data_Science_Project/gcta-1.94.2-MacOS-ARM-x86_64/gcclib:$DYLD_LIBRARY_PATH"'
system(lib_path, intern = TRUE)  # Set library path
system(
  paste0(tool_path, "gcta64 ", "--bfile ", bfile_path, "scfc_all_qc_unrelated ",
      "--autosome  --make-grm  ", "--thread-num 16 ",
      "--out ", bfile_path, "scfc_all_qc_unrelated_out"), intern = F)
system(
  paste0(tool_path, "gcta64 ", "--grm ", bfile_path, "scfc_all_qc_out ",
      "--make-bK-sparse 0.05 ", "--thread-num 16 ",
      "--out ", bfile_path, "scfc_all_qc_unrelated_sp_out"), intern = F)


## II. Run GWAS commands
# Modify these parameters below to run the GWAS commands you want (name your methods first)
# for reference: gwas_methods <- c("gwasplink19", "gwasplink20_nocov", "gwasplink20_cov", "fastGWA_nocov", "fastGWA_cov")
QCs <- c("_noqc", "_qc") # qc and no qc
unrels <- c("", "_unrelated") # all obs and unrelated obs
harms <- c("_harm2", "_sd_harm2") # harm2 and fisher+harm2

for (qc in QCs) {
  for (unrel in unrels) {
    #-----------get bfile---------------------------------------------------------
    bfile <- paste0(bfile_path, "scfc_all", qc, unrel)
    #-------------get grm matrix(for fastGWA)-----------------------------------
    bgrm <- paste0(bfile_path, "scfc_all", qc, unrel, "_out")
    bspr <- paste0(bfile_path, "scfc_all", qc, unrel, "_sp_out")
    
    for (harm in harms) {
      ##########For both all and unrelated obs, 4 methods for each combination##
      
      
      #=========1. GWAS with PLINK1.9, no covariates (--assoc)==================
      method <- "gwasplink19"
      #+++++++++Set work directory++++++++++++++++++++++++++++++++++++++++++++++
      work_dir <- paste0(directory, method, qc, unrel, harm, "/")
      if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(work_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("Commands run for GWAS analysis:", method, qc, unrel, harm), readme_conn)
      start_time <- Sys.time()
      #---------Run GWAS commands for ROIs--------------------------------------
      cat(paste0("\n***************Running GWAS: ", method, qc, unrel, harm, "***************\n\n"))
      for (roi in rois) {
        #+++++++++get SCFC ROI phenotypes+++++++++++++++++++++++++++++++++++++++
        scfc_roi <- paste0(pheno_path, "scfc", harm, "_fam_roi", roi, ".phen")
        #+++++++++Get covariates++++++++++++++++++++++++++++++++++++++++++++++++
        # No covar is needed for current method
        #+++++++++Name output file name+++++++++++++++++++++++++++++++++++++++++
        outfile <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi)
        #+++++++++Run GWAS commands+++++++++++++++++++++++++++++++++++++++++++++
        command <- paste0(tool_path, "plink --bfile ", bfile, " ",
                          "--pheno ", scfc_roi, " ", 
                          "--assoc --allow-no-sex ",  # Added hide-covar to avoid covariate effect lines in the output
                          "--out ", outfile)
        system(command)
        cat(method, qc, unrel, harm, " executed successfully for:", "scfc", harm, "_fam_roi", roi, "\n")
        writeLines(command, readme_conn, useBytes = TRUE)
      }
      end_time <- Sys.time()
      writeLines(sprintf("Total time to run all commands: %s minutes", difftime(end_time, start_time, units = "mins")), readme_conn)
      close(readme_conn)
      
      #---------Save GWAS results-----------------------------------------------
      #+++++++++Set save directory++++++++++++++++++++++++++++++++++++++++++++++
      save_dir <- paste0(work_dir, "gwas_results_csv/")
      if (!dir.exists(save_dir)) { dir.create(save_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(save_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("GWAS using ", method, qc, harm, unrel, "\n\n"), readme_conn)
      writeLines(paste0("Commands run: ", command, "\n\n"), readme_conn)
      start_time <- Sys.time()
      missing <- 0
      #+++++++++Create SNP counts+++++++++++++++++++++++++++++++++++++++++++++++
      count_path <- paste0(save_dir, "SNP_counts.txt")
      count_conn <- file(count_path, "w")
      writeLines(paste0("SNPs P<5e-08 for GWAS using ", method, qc, harm, unrel, "\n\n"), count_conn)
      #---------Read results for ROIs-------------------------------------------
      cat(paste0("\n^^^^^^^^^^^^Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n"))
      for (roi in rois) {
        gwas_result_path <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".qassoc")
        if (file.exists(gwas_result_path)) {
          gwas_result <- read.table(gwas_result_path, header = TRUE, sep = "")
          write.csv(gwas_result, paste0(save_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".csv"), row.names = FALSE, quote = FALSE)
          writeLines(paste0("\n\nSignificant SNPs count for ROI ", roi, ": ", nrow(gwas_result %>% filter(P < 5e-08)), " for GWAS using ", method, qc, harm, unrel), count_conn)
        } else {
          missing = missing + 1
          writeLines(paste0("Missing file for ROI ", roi, ": ", gwas_result_path), readme_conn)
          cat("No GWAS results: ", gwas_result_path, "\n")  # Print message to the console as well
        }
      }
      cat(paste0("^^^^^^^^^^^^Finished Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n\n"))
      
      end_time <- Sys.time()
      total_duration <- difftime(end_time, start_time, units = "mins")
      writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
      writeLines(paste("Total missing files:", missing), readme_conn)
      close(readme_conn)
      close(count_conn)
      cat(paste0("***************Finished running commands:", method, qc, unrel, harm, "***************\n\n\n"))
      #=========================================================================
      
      
      #=========2. GWAS with PLINK2.0, no covariates (--glm)====================
      method <- "gwasplink20_nocov"
      #+++++++++Set work directory++++++++++++++++++++++++++++++++++++++++++++++
      work_dir <- paste0(directory, method, qc, unrel, harm, "/")
      if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(work_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("Commands run for GWAS analysis:", method, qc, unrel, harm), readme_conn)
      start_time <- Sys.time()
      #---------Run GWAS commands for ROIs--------------------------------------
      cat(paste0("\n***************Running GWAS: ", method, qc, unrel, harm, "***************\n\n"))
      for (roi in rois) {
        #+++++++++get SCFC ROI phenotypes+++++++++++++++++++++++++++++++++++++++
        scfc_roi <- paste0(pheno_path, "scfc", harm, "_fam_roi", roi, ".phen")
        #+++++++++Get covariates++++++++++++++++++++++++++++++++++++++++++++++++
        # No covar is needed for current method
        #+++++++++Name output file name+++++++++++++++++++++++++++++++++++++++++
        outfile <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi)
        #+++++++++Run GWAS commands+++++++++++++++++++++++++++++++++++++++++++++
        command <- paste0(tool_path, "plink2 --bfile ", bfile, " ",
                          "--pheno ", scfc_roi, " ",
                          "--glm allow-no-covars no-x-sex hide-covar ",  # Added hide-covar to avoid covariate effect lines in the output
                          "--covar-variance-standardize ",
                          "--out ", outfile)
        system(command)
        cat(method, qc, unrel, harm, " executed successfully for:", "scfc", harm, "_fam_roi", roi, "\n")
        writeLines(command, readme_conn, useBytes = TRUE)
      }
      end_time <- Sys.time()
      writeLines(sprintf("Total time to run all commands: %s minutes", difftime(end_time, start_time, units = "mins")), readme_conn)
      close(readme_conn)
      
      #---------Save GWAS results-----------------------------------------------
      #+++++++++Set save directory++++++++++++++++++++++++++++++++++++++++++++++
      save_dir <- paste0(work_dir, "gwas_results_csv/")
      if (!dir.exists(save_dir)) { dir.create(save_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(save_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("GWAS using ", method, qc, harm, unrel, "\n\n"), readme_conn)
      writeLines(paste0("Commands run: ", command, "\n\n"), readme_conn)
      start_time <- Sys.time()
      missing <- 0
      #+++++++++Create SNP counts+++++++++++++++++++++++++++++++++++++++++++++++
      count_path <- paste0(save_dir, "SNP_counts.txt")
      count_conn <- file(count_path, "w")
      writeLines(paste0("SNPs P<5e-08 for GWAS using ", method, qc, harm, unrel, "\n\n"), count_conn)
      #---------Read results for ROIs-------------------------------------------
      cat(paste0("\n^^^^^^^^^^^^Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n"))
      for (roi in rois) {
        gwas_result_path <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".PHENO1.glm.linear")
        if (file.exists(gwas_result_path)) {
          gwas_result <- read.table(gwas_result_path, header = F, sep = "\t")
          colnames(gwas_result) <- c("#CHROM", "POS", "SNP", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
          write.csv(gwas_result, paste0(save_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".csv"), row.names = FALSE, quote = FALSE)
          writeLines(paste0("\n\nSignificant SNPs count for ROI ", roi, ": ", nrow(gwas_result %>% filter(P < 5e-08)), " for GWAS using ", method, qc, harm, unrel), count_conn)
        } else {
          missing = missing + 1
          writeLines(paste0("Missing file for ROI ", roi, ": ", gwas_result_path), readme_conn)
          cat("No GWAS results: ", gwas_result_path, "\n")  # Print message to the console as well
        }
      }
      cat(paste0("^^^^^^^^^^^^Finished Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n\n"))
      
      end_time <- Sys.time()
      total_duration <- difftime(end_time, start_time, units = "mins")
      writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
      writeLines(paste("Total missing files:", missing), readme_conn)
      close(readme_conn)
      close(count_conn)
      cat(paste0("***************Finished running commands:", method, qc, unrel, harm, "***************\n\n\n"))
      #=========================================================================  
      
      
      #=========3. GWAS with PLINK2.0, with covariates (--glm)==================
      method <- "gwasplink20_cov"
      #+++++++++Set work directory++++++++++++++++++++++++++++++++++++++++++++++
      work_dir <- paste0(directory, method, qc, unrel, harm, "/")
      if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(work_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("Commands run for GWAS analysis:", method, qc, unrel, harm), readme_conn)
      start_time <- Sys.time()
      #---------Run GWAS commands for ROIs--------------------------------------
      cat(paste0("\n***************Running GWAS: ", method, qc, unrel, harm, "***************\n\n"))
      for (roi in rois) {
        #+++++++++get SCFC ROI phenotypes+++++++++++++++++++++++++++++++++++++++
        scfc_roi <- paste0(pheno_path, "scfc", harm, "_fam_roi", roi, ".phen")
        #+++++++++Get covariates++++++++++++++++++++++++++++++++++++++++++++++++
        covar_all <- paste0(cov_path, "covar_all.txt")
        #+++++++++Name output file name+++++++++++++++++++++++++++++++++++++++++
        outfile <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi)
        #+++++++++Run GWAS commands+++++++++++++++++++++++++++++++++++++++++++++
        command <- paste0(tool_path, "plink2 --bfile ", bfile, " ",
                          "--pheno ", scfc_roi, " ",
                          "--covar ", covar_all, " ",
                          "--glm no-x-sex hide-covar ",  # Added hide-covar to avoid covariate effect lines in the output
                          "--covar-variance-standardize ",
                          "--out ", outfile)
        system(command)
        cat(method, qc, unrel, harm, " executed successfully for:", "scfc", harm, "_fam_roi", roi, "\n")
        writeLines(command, readme_conn, useBytes = TRUE)
      }
      end_time <- Sys.time()
      writeLines(sprintf("Total time to run all commands: %s minutes", difftime(end_time, start_time, units = "mins")), readme_conn)
      close(readme_conn)
      
      #---------Save GWAS results-----------------------------------------------
      #+++++++++Set save directory++++++++++++++++++++++++++++++++++++++++++++++
      save_dir <- paste0(work_dir, "gwas_results_csv/")
      if (!dir.exists(save_dir)) { dir.create(save_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(save_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("GWAS using ", method, qc, harm, unrel, "\n\n"), readme_conn)
      writeLines(paste0("Commands run: ", command, "\n\n"), readme_conn)
      start_time <- Sys.time()
      missing <- 0
      #+++++++++Create SNP counts+++++++++++++++++++++++++++++++++++++++++++++++
      count_path <- paste0(save_dir, "SNP_counts.txt")
      count_conn <- file(count_path, "w")
      writeLines(paste0("SNPs P<5e-08 for GWAS using ", method, qc, harm, unrel, "\n\n"), count_conn)
      #---------Read results for ROIs-------------------------------------------
      cat(paste0("\n^^^^^^^^^^^^Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n"))
      for (roi in rois) {
        gwas_result_path <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".PHENO1.glm.linear")
        if (file.exists(gwas_result_path)) {
          gwas_result <- read.table(gwas_result_path, header = F, sep = "\t")
          colnames(gwas_result) <- c("#CHROM", "POS", "SNP", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
          write.csv(gwas_result, paste0(save_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".csv"), row.names = FALSE, quote = FALSE)
          writeLines(paste0("\n\nSignificant SNPs count for ROI ", roi, ": ", nrow(gwas_result %>% filter(P < 5e-08)), " for GWAS using ", method, qc, harm, unrel), count_conn)
        } else {
          missing = missing + 1
          writeLines(paste0("Missing file for ROI ", roi, ": ", gwas_result_path), readme_conn)
          cat("No GWAS results: ", gwas_result_path, "\n")  # Print message to the console as well
        }
      }
      cat(paste0("^^^^^^^^^^^^Finished Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n\n"))
      
      end_time <- Sys.time()
      total_duration <- difftime(end_time, start_time, units = "mins")
      writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
      writeLines(paste("Total missing files:", missing), readme_conn)
      close(readme_conn)
      close(count_conn)
      cat(paste0("***************Finished running commands:", method, qc, unrel, harm, "***************\n\n\n"))
      #=========================================================================
      
      
      #=========4. GWAS with fastGWA, no covariates=============================
      method <- "fastGWA_nocov"
      #+++++++++Set work directory++++++++++++++++++++++++++++++++++++++++++++++
      work_dir <- paste0(directory, method, qc, unrel, harm, "/")
      if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(work_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("Commands run for GWAS analysis:", method, qc, unrel, harm), readme_conn)
      start_time <- Sys.time()
      #---------Run GWAS commands for ROIs--------------------------------------
      cat(paste0("\n***************Running GWAS: ", method, qc, unrel, harm, "***************\n\n"))
      for (roi in rois) {
        #+++++++++get SCFC ROI phenotypes+++++++++++++++++++++++++++++++++++++++
        scfc_roi <- paste0(pheno_path, "scfc", harm, "_fam_roi", roi, ".phen")
        #+++++++++Get covariates++++++++++++++++++++++++++++++++++++++++++++++++
        # No covar is needed for current method
        #+++++++++Name output file name+++++++++++++++++++++++++++++++++++++++++
        outfile <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi)
        #+++++++++Run GWAS commands+++++++++++++++++++++++++++++++++++++++++++++
        command <- paste0(tool_path, "gcta64 --bfile ", bfile, " ",
                          "--pheno ", scfc_roi, " ",
                          "--grm-sparse ", bspr, " ",
                          "--fastGWA-mlm ",
                          "--thread-num 16 ", "--out ", outfile)
        system(command)
        cat(method, qc, unrel, harm, " executed successfully for:", "scfc", harm, "_fam_roi", roi, "\n")
        writeLines(command, readme_conn, useBytes = TRUE)
      }
      end_time <- Sys.time()
      writeLines(sprintf("Total time to run all commands: %s minutes", difftime(end_time, start_time, units = "mins")), readme_conn)
      close(readme_conn)
      
      #---------Save GWAS results-----------------------------------------------
      #+++++++++Set save directory++++++++++++++++++++++++++++++++++++++++++++++
      save_dir <- paste0(work_dir, "gwas_results_csv/")
      if (!dir.exists(save_dir)) { dir.create(save_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(save_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("GWAS using ", method, qc, harm, unrel, "\n\n"), readme_conn)
      writeLines(paste0("Commands run: ", command, "\n\n"), readme_conn)
      start_time <- Sys.time()
      missing <- 0
      #+++++++++Create SNP counts+++++++++++++++++++++++++++++++++++++++++++++++
      count_path <- paste0(save_dir, "SNP_counts.txt")
      count_conn <- file(count_path, "w")
      writeLines(paste0("SNPs P<5e-08 for GWAS using ", method, qc, harm, unrel, "\n\n"), count_conn)
      #---------Read results for ROIs-------------------------------------------
      cat(paste0("\n^^^^^^^^^^^^Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n"))
      for (roi in rois) {
        gwas_result_path <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".fastGWA")
        if (file.exists(gwas_result_path)) {
          gwas_result <- read.table(gwas_result_path, header = TRUE, sep = "")
          write.csv(gwas_result, paste0(save_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".csv"), row.names = FALSE, quote = FALSE)
          writeLines(paste0("\n\nSignificant SNPs count for ROI ", roi, ": ", nrow(gwas_result %>% filter(P < 5e-08)), " for GWAS using ", method, qc, harm, unrel), count_conn)
        } else {
          missing = missing + 1
          writeLines(paste0("Missing file for ROI ", roi, ": ", gwas_result_path), readme_conn)
          cat("No GWAS results: ", gwas_result_path, "\n")  # Print message to the console as well
        }
      }
      cat(paste0("^^^^^^^^^^^^Finished Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n\n"))
      
      end_time <- Sys.time()
      total_duration <- difftime(end_time, start_time, units = "mins")
      writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
      writeLines(paste("Total missing files:", missing), readme_conn)
      close(readme_conn)
      close(count_conn)
      cat(paste0("***************Finished running commands:", method, qc, unrel, harm, "***************\n\n\n"))
      #=========================================================================
      
      
      #=========5. GWAS with fastGWA, with covariates===========================
      method <- "fastGWA_cov"
      #+++++++++Set work directory++++++++++++++++++++++++++++++++++++++++++++++
      work_dir <- paste0(directory, method, qc, unrel, harm, "/")
      if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(work_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("Commands run for GWAS analysis:", method, qc, unrel, harm), readme_conn)
      start_time <- Sys.time()
      #---------Run GWAS commands for ROIs--------------------------------------
      cat(paste0("\n***************Running GWAS: ", method, qc, unrel, harm, "***************\n\n"))
      for (roi in rois) {
        #+++++++++get SCFC ROI phenotypes+++++++++++++++++++++++++++++++++++++++
        scfc_roi <- paste0(pheno_path, "scfc", harm, "_fam_roi", roi, ".phen")
        #+++++++++Get covariates++++++++++++++++++++++++++++++++++++++++++++++++
        qcovar_age <- paste0(cov_path, "qcovar_age.txt")
        covar_sex_race_ethn <- paste0(cov_path, "covar_sex_race_ethn.txt")
        #+++++++++Name output file name+++++++++++++++++++++++++++++++++++++++++
        outfile <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi)
        #+++++++++Run GWAS commands+++++++++++++++++++++++++++++++++++++++++++++
        command <- paste0(tool_path, "gcta64 --bfile ", bfile, " ",
                          "--pheno ", scfc_roi, " ",
                          "--grm-sparse ", bspr, " ",
                          "--fastGWA-mlm ",
                          "--qcovar ", qcovar_age, " ",
                          "--covar ", covar_sex_race_ethn, " ",
                          "--thread-num 16 ", "--out ", outfile)
        system(command)
        cat(method, qc, unrel, harm, " executed successfully for:", "scfc", harm, "_fam_roi", roi, "\n")
        writeLines(command, readme_conn, useBytes = TRUE)
      }
      end_time <- Sys.time()
      writeLines(sprintf("Total time to run all commands: %s minutes", difftime(end_time, start_time, units = "mins")), readme_conn)
      close(readme_conn)
      
      #---------Save GWAS results-----------------------------------------------
      #+++++++++Set save directory++++++++++++++++++++++++++++++++++++++++++++++
      save_dir <- paste0(work_dir, "gwas_results_csv/")
      if (!dir.exists(save_dir)) { dir.create(save_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(save_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("GWAS using ", method, qc, harm, unrel, "\n\n"), readme_conn)
      writeLines(paste0("Commands run: ", command, "\n\n"), readme_conn)
      start_time <- Sys.time()
      missing <- 0
      #+++++++++Create SNP counts+++++++++++++++++++++++++++++++++++++++++++++++
      count_path <- paste0(save_dir, "SNP_counts.txt")
      count_conn <- file(count_path, "w")
      writeLines(paste0("SNPs P<5e-08 for GWAS using ", method, qc, harm, unrel, "\n\n"), count_conn)
      #---------Read results for ROIs-------------------------------------------
      cat(paste0("\n^^^^^^^^^^^^Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n"))
      for (roi in rois) {
        gwas_result_path <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".fastGWA")
        if (file.exists(gwas_result_path)) {
          gwas_result <- read.table(gwas_result_path, header = TRUE, sep = "")
          write.csv(gwas_result, paste0(save_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".csv"), row.names = FALSE, quote = FALSE)
          writeLines(paste0("\n\nSignificant SNPs count for ROI ", roi, ": ", nrow(gwas_result %>% filter(P < 5e-08)), " for GWAS using ", method, qc, harm, unrel), count_conn)
        } else {
          missing = missing + 1
          writeLines(paste0("Missing file for ROI ", roi, ": ", gwas_result_path), readme_conn)
          cat("No GWAS results: ", gwas_result_path, "\n")  # Print message to the console as well
        }
      }
      cat(paste0("^^^^^^^^^^^^Finished Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n\n"))
      
      end_time <- Sys.time()
      total_duration <- difftime(end_time, start_time, units = "mins")
      writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
      writeLines(paste("Total missing files:", missing), readme_conn)
      close(readme_conn)
      close(count_conn)
      cat(paste0("***************Finished running commands:", method, qc, unrel, harm, "***************\n\n\n"))
      #=========================================================================
    }
  }
}


# for reference: gwas_methods <- c("gwasplink19", "gwasplink20_nocov", "gwasplink20_cov", "fastGWA_nocov", "fastGWA_cov")



## 1008: Test with covariates

QCs <- c("_noqc", "_qc") # qc and no qc
unrels <- c("", "_unrelated") # all obs and unrelated obs
harms <- c("_sd_harm2") # harm2 and fisher+harm2
# for reference: gwas_methods <- c("gwasplink19", "gwasplink20_nocov", "gwasplink20_cov", "fastGWA_nocov", "fastGWA_cov")

for (qc in QCs) {
  for (unrel in unrels) {
    #-----------get bfile---------------------------------------------------------
    bfile <- paste0(bfile_path, "scfc_all", qc, unrel)
    #-------------get grm matrix(for fastGWA)-----------------------------------
    bgrm <- paste0(bfile_path, "scfc_all", qc, unrel, "_out")
    bspr <- paste0(bfile_path, "scfc_all", qc, unrel, "_sp_out")
    
    for (harm in harms) {
      ##########For both all and unrelated obs, 4 methods for each combination##
      
      
      #=========1. GWAS with PLINK2.0, with sex covariate (--glm)===============
      method <- "gwasplink20_covsex"
      #+++++++++Set work directory++++++++++++++++++++++++++++++++++++++++++++++
      work_dir <- paste0(directory, method, qc, unrel, harm, "/")
      if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(work_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("Commands run for GWAS analysis:", method, qc, unrel, harm), readme_conn)
      start_time <- Sys.time()
      #---------Run GWAS commands for ROIs--------------------------------------
      cat(paste0("\n***************Running GWAS: ", method, qc, unrel, harm, "***************\n\n"))
      for (roi in rois) {
        #+++++++++get SCFC ROI phenotypes+++++++++++++++++++++++++++++++++++++++
        scfc_roi <- paste0(pheno_path, "scfc", harm, "_fam_roi", roi, ".phen")
        #+++++++++Get covariates++++++++++++++++++++++++++++++++++++++++++++++++
        covar_sex <- paste0(cov_path, "covar_sex.txt")
        #+++++++++Name output file name+++++++++++++++++++++++++++++++++++++++++
        outfile <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi)
        #+++++++++Run GWAS commands+++++++++++++++++++++++++++++++++++++++++++++
        command <- paste0(tool_path, "plink2 --bfile ", bfile, " ",
                          "--pheno ", scfc_roi, " ",
                          "--covar ", covar_sex, " ",
                          "--glm no-x-sex hide-covar ",  # Added hide-covar to avoid covariate effect lines in the output
                          "--covar-variance-standardize ",
                          "--out ", outfile)
        system(command)
        cat(method, qc, unrel, harm, " executed successfully for:", "scfc", harm, "_fam_roi", roi, "\n")
        writeLines(command, readme_conn, useBytes = TRUE)
      }
      end_time <- Sys.time()
      writeLines(sprintf("Total time to run all commands: %s minutes", difftime(end_time, start_time, units = "mins")), readme_conn)
      close(readme_conn)
      
      #---------Save GWAS results-----------------------------------------------
      #+++++++++Set save directory++++++++++++++++++++++++++++++++++++++++++++++
      save_dir <- paste0(work_dir, "gwas_results_csv/")
      if (!dir.exists(save_dir)) { dir.create(save_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(save_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("GWAS using ", method, qc, harm, unrel, "\n\n"), readme_conn)
      writeLines(paste0("Commands run: ", command, "\n\n"), readme_conn)
      start_time <- Sys.time()
      missing <- 0
      #+++++++++Create SNP counts+++++++++++++++++++++++++++++++++++++++++++++++
      count_path <- paste0(save_dir, "SNP_counts.txt")
      count_conn <- file(count_path, "w")
      writeLines(paste0("SNPs P<5e-08 for GWAS using ", method, qc, harm, unrel, "\n\n"), count_conn)
      #---------Read results for ROIs-------------------------------------------
      cat(paste0("\n^^^^^^^^^^^^Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n"))
      for (roi in rois) {
        gwas_result_path <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".PHENO1.glm.linear")
        if (file.exists(gwas_result_path)) {
          gwas_result <- read.table(gwas_result_path, header = F, sep = "\t")
          colnames(gwas_result) <- c("#CHROM", "POS", "SNP", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
          write.csv(gwas_result, paste0(save_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".csv"), row.names = FALSE, quote = FALSE)
          writeLines(paste0("\n\nSignificant SNPs count for ROI ", roi, ": ", nrow(gwas_result %>% filter(P < 5e-08)), " for GWAS using ", method, qc, harm, unrel), count_conn)
        } else {
          missing = missing + 1
          writeLines(paste0("Missing file for ROI ", roi, ": ", gwas_result_path), readme_conn)
          cat("No GWAS results: ", gwas_result_path, "\n")  # Print message to the console as well
        }
      }
      cat(paste0("^^^^^^^^^^^^Finished Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n\n"))
      
      end_time <- Sys.time()
      total_duration <- difftime(end_time, start_time, units = "mins")
      writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
      writeLines(paste("Total missing files:", missing), readme_conn)
      close(readme_conn)
      close(count_conn)
      cat(paste0("***************Finished running commands:", method, qc, unrel, harm, "***************\n\n\n"))
      #=========================================================================
      
      
      #=========2. GWAS with fastGWA, with sex covariate========================
      method <- "fastGWA_covsex"
      #+++++++++Set work directory++++++++++++++++++++++++++++++++++++++++++++++
      work_dir <- paste0(directory, method, qc, unrel, harm, "/")
      if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(work_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("Commands run for GWAS analysis:", method, qc, unrel, harm), readme_conn)
      start_time <- Sys.time()
      #---------Run GWAS commands for ROIs--------------------------------------
      cat(paste0("\n***************Running GWAS: ", method, qc, unrel, harm, "***************\n\n"))
      for (roi in rois) {
        #+++++++++get SCFC ROI phenotypes+++++++++++++++++++++++++++++++++++++++
        scfc_roi <- paste0(pheno_path, "scfc", harm, "_fam_roi", roi, ".phen")
        #+++++++++Get covariates++++++++++++++++++++++++++++++++++++++++++++++++
        covar_sex <- paste0(cov_path, "covar_sex.txt")
        #+++++++++Name output file name+++++++++++++++++++++++++++++++++++++++++
        outfile <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi)
        #+++++++++Run GWAS commands+++++++++++++++++++++++++++++++++++++++++++++
        command <- paste0(tool_path, "gcta64 --bfile ", bfile, " ",
                          "--pheno ", scfc_roi, " ",
                          "--grm-sparse ", bspr, " ",
                          "--fastGWA-mlm ",
                          "--covar ", covar_sex, " ",
                          "--thread-num 16 ", "--out ", outfile)
        system(command)
        cat(method, qc, unrel, harm, " executed successfully for:", "scfc", harm, "_fam_roi", roi, "\n")
        writeLines(command, readme_conn, useBytes = TRUE)
      }
      end_time <- Sys.time()
      writeLines(sprintf("Total time to run all commands: %s minutes", difftime(end_time, start_time, units = "mins")), readme_conn)
      close(readme_conn)
      
      #---------Save GWAS results-----------------------------------------------
      #+++++++++Set save directory++++++++++++++++++++++++++++++++++++++++++++++
      save_dir <- paste0(work_dir, "gwas_results_csv/")
      if (!dir.exists(save_dir)) { dir.create(save_dir, recursive = F) }
      #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
      readme_path <- paste0(save_dir, "readme.txt")
      readme_conn <- file(readme_path, "w")
      writeLines(paste0("GWAS using ", method, qc, harm, unrel, "\n\n"), readme_conn)
      writeLines(paste0("Commands run: ", command, "\n\n"), readme_conn)
      start_time <- Sys.time()
      missing <- 0
      #+++++++++Create SNP counts+++++++++++++++++++++++++++++++++++++++++++++++
      count_path <- paste0(save_dir, "SNP_counts.txt")
      count_conn <- file(count_path, "w")
      writeLines(paste0("SNPs P<5e-08 for GWAS using ", method, qc, harm, unrel, "\n\n"), count_conn)
      #---------Read results for ROIs-------------------------------------------
      cat(paste0("\n^^^^^^^^^^^^Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n"))
      for (roi in rois) {
        gwas_result_path <- paste0(work_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".fastGWA")
        if (file.exists(gwas_result_path)) {
          gwas_result <- read.table(gwas_result_path, header = TRUE, sep = "")
          write.csv(gwas_result, paste0(save_dir, method, qc, unrel, "_scfc", harm, "_fam_roi", roi, ".csv"), row.names = FALSE, quote = FALSE)
          writeLines(paste0("\n\nSignificant SNPs count for ROI ", roi, ": ", nrow(gwas_result %>% filter(P < 5e-08)), " for GWAS using ", method, qc, harm, unrel), count_conn)
        } else {
          missing = missing + 1
          writeLines(paste0("Missing file for ROI ", roi, ": ", gwas_result_path), readme_conn)
          cat("No GWAS results: ", gwas_result_path, "\n")  # Print message to the console as well
        }
      }
      cat(paste0("^^^^^^^^^^^^Finished Reading GWAS results: ", method, qc, unrel, harm, "^^^^^^^^^^^^\n\n"))
      
      end_time <- Sys.time()
      total_duration <- difftime(end_time, start_time, units = "mins")
      writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
      writeLines(paste("Total missing files:", missing), readme_conn)
      close(readme_conn)
      close(count_conn)
      cat(paste0("***************Finished running commands:", method, qc, unrel, harm, "***************\n\n\n"))
      #=========================================================================
    }
  }
}
