#

# Running GWAS 

# 4. Script for running GWAS for ADHD scores (for case/control Reverse MR)

###############################################################################################
## Readme: this code runs and saves GWAS for case/control reverse MR analysis, using PLINK
###############################################################################################


library(snpStats)
library(tibble)
library(dplyr)
library(foreach)
library(doParallel)

dsm5diags <- c("depress", "anxdisord", "somaticpr", "adhd", "opposit", "conduct")
dsm5diags <- c("ADHD")

### Load Utilities, data, tool and define directory
utilities_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_SCFC_utilities/"
bfile_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_bfile/bfile_for_gwas_comparison/"
cov_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Covariates/cov_for_gwas_comparison/"
harm_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Harmonized_scfc/"
pheno_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Phenotype/"
tool_path <- "/Users/grc8mu/Desktop/DS/Tools/"

directory <- "/Users/grc8mu/Desktop/DS/1102_GWAS_Diagnosis/"

load(paste0(utilities_path, "coupling_cbcl.RData"))
load(paste0(utilities_path, "scfc_coupling_matched.RData"))
new_fam <- read.table(paste0(utilities_path, "acspsw03.txt"), header = T)
new_fam <- new_fam %>%
  select(rel_family_id, subjectkey, eventname, race_ethnicity)

# read covariates dataset
CFDRs <- demo_cbcl_base_scfc %>%
  select(subjectkey, site_id_l, eventname, interview_age, sex, race, ethn) %>%
  mutate(
    subjectkey = sub("^(.{4})", "\\1_", subjectkey)  # Add '_' after the 4th character
    #ethn = as.factor(ethn),  # Convert ethnicity to a factor
    #race = as.factor(race)   # Convert race to a factor if needed
  ) %>%
  inner_join(new_fam, by = c("subjectkey", "eventname")) %>%
  select(rel_family_id, subjectkey, site_id_l, interview_age, sex, race, ethn, race_ethnicity) %>%
  rename(fid = rel_family_id, iid = subjectkey) %>%
  mutate(  # Add underscore to subjectkey
    race = gsub(" ", "_", race),                      # Replace spaces in 'race'
    ethn = gsub(" ", "_", ethn)                       # Replace spaces in 'ethnicity'
  )


## Currently, only need to create new phenotypes for ADHD (Created in the case/control folder)
qc <- "_qc" # 
QCs <- c("_qc") # qc and no qc
unrel <- "" # 
unrels <- c("", "_unrelated") # all obs and unrelated obs
# for reference: gwas_methods <- c("gwasplink19", "gwasplink20_nocov", "gwasplink20_cov", "fastGWA_nocov", "fastGWA_cov")

for (qc in QCs) {
  for (unrel in unrels) {
    #-----------get bfile---------------------------------------------------------
    bfile <- paste0(bfile_path, "scfc_all", qc, unrel)
    #-------------get grm matrix(for fastGWA)-----------------------------------
    bgrm <- paste0(bfile_path, "scfc_all", qc, unrel, "_out")
    bspr <- paste0(bfile_path, "scfc_all", qc, unrel, "_sp_out")
    
    ##########For both all and unrelated obs, 4 methods for each combination##
    
    ########1. GWAS with PLINK1.9, with no covariate (--assoc)#########
    method <- "dsm5diag_gwasplink19"
    #+++++++++Get covariates++++++++++++++++++++++++++++++++++++++++++++++++
    #covar_sex <- paste0(cov_path, "covar_sex.txt")
    #+++++++++Set work directory++++++++++++++++++++++++++++++++++++++++++++++
    work_dir <- paste0(directory, method, qc, unrel, "/")
    if (!dir.exists(work_dir)) { dir.create(work_dir, recursive = F) }
    #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
    readme_path <- paste0(work_dir, "readme.txt")
    readme_conn <- file(readme_path, "w")
    writeLines(paste0("Commands run for GWAS analysis:", method, qc, unrel), readme_conn)
    start_time <- Sys.time()
    #---------Run GWAS commands for dsm5diags--------------------------------------
    cat(paste0("\n***************Running GWAS: ", method, qc, unrel, "***************\n\n"))
    # Set up parallel backend
    num_cores <- detectCores() - 1
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    foreach(dsm5diag = dsm5diags, .combine = list) %dopar% {
      tryCatch({
        #+++++++++get SCFC dsm5diag phenotypes+++++++++++++++++++++++++++++++++++++++
        dsm5diag_pheno <- paste0(pheno_path, "dsm5diag_", dsm5diag, ".phen")
        #+++++++++Name output file name+++++++++++++++++++++++++++++++++++++++++
        outfile <- paste0(work_dir, method, qc, unrel, "_dsm5diag_", dsm5diag)
        #+++++++++Run GWAS commands+++++++++++++++++++++++++++++++++++++++++++++
        command <- paste0(tool_path, "plink --bfile ", bfile, " ",
                          "--allow-no-sex ",
                          "--pheno ", dsm5diag_pheno, " ",
                          "--logistic ",  # Added hide-covar to avoid covariate effect lines in the output
                          "--out ", outfile)
        system(command)
        cat(method, qc, unrel, " executed successfully for:", "dsm5diag_", dsm5diag, "\n")
        writeLines(command, readme_conn, useBytes = TRUE)}, error = function(e) {
          return(paste("Error:", e$message))
        })
    }
    end_time <- Sys.time()
    writeLines(sprintf("Total time to run all commands: %s minutes", difftime(end_time, start_time, units = "mins")), readme_conn)
    close(readme_conn)
    stopCluster(cl)
    
    #---------Save GWAS results-----------------------------------------------
    #+++++++++Set save directory++++++++++++++++++++++++++++++++++++++++++++++
    save_dir <- paste0(work_dir, "gwas_results_csv/")
    if (!dir.exists(save_dir)) { dir.create(save_dir, recursive = F) }
    #+++++++++Create readme+++++++++++++++++++++++++++++++++++++++++++++++++++
    readme_path <- paste0(save_dir, "readme.txt")
    readme_conn <- file(readme_path, "w")
    writeLines(paste0("GWAS using ", method, qc, unrel, "\n\n"), readme_conn)
    command <- paste0(tool_path, "plink --bfile ", bfile, " ",
                      "--allow-no-sex ", 
                      "--pheno ", paste0(pheno_path, "dsm5diag_", "adhd", ".phen"), " ",
                      "--logistic ",  # Added hide-covar to avoid covariate effect lines in the output
                      "--out ", paste0(work_dir, method, qc, unrel, "_dsm5diag_", "adhd"))
    writeLines(paste0("Commands run: ", command, "\n\n"), readme_conn)
    start_time <- Sys.time()
    missing <- 0
    #+++++++++Create SNP counts+++++++++++++++++++++++++++++++++++++++++++++++
    count_path <- paste0(save_dir, "SNP_counts.txt")
    count_conn <- file(count_path, "w")
    writeLines(paste0("SNPs P<5e-08 for GWAS using ", method, qc, unrel, "\n\n"), count_conn)
    #---------Read results for dsm5diags-------------------------------------------
    cat(paste0("\n^^^^^^^^^^^^Reading GWAS results: ", method, qc, unrel, "^^^^^^^^^^^^\n"))
    num_cores <- detectCores() - 1
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    # parallel read gwas results
    parallel_results <- foreach(dsm5diag = dsm5diags, .packages = c("dplyr")) %dopar% {
      library(dplyr)
      tryCatch({
        gwas_result_path <- paste0(work_dir, method, qc, unrel, "_dsm5diag_", dsm5diag,  ".assoc.logistic")
        if (file.exists(gwas_result_path)) {
          gwas_result <- read.table(gwas_result_path, header = T)
          write.csv(gwas_result, paste0(save_dir, method, qc, unrel, "_dsm5diag_", dsm5diag,  ".csv"), row.names = FALSE, quote = FALSE)
          sig_snp_count <- nrow(gwas_result %>% filter(P < 5e-08))
          return(list(dsm5diag = dsm5diag, missing = FALSE, snp_count = sig_snp_count))
        } else {
          return(list(dsm5diag = dsm5diag, missing = TRUE, snp_count = NA))
        }}, error = function(e) {
          return(list(dsm5diag = dsm5diag, error = paste("Error:", e$message), missing = TRUE, snp_count = NA))
        })
    }
    stopCluster(cl)
    cat(paste0("^^^^^^^^^^^^Finished Reading GWAS results: ", method, qc, unrel, "^^^^^^^^^^^^\n\n"))
    
    for (res in parallel_results) {
      if (res$missing) {
        writeLines(paste0("Missing file for dsm5diag ", res$dsm5diag, ": ", paste0(work_dir, method, qc, unrel, "_dsm5diag_", res$dsm5diag, ".assoc.logistic")), readme_conn)
        cat("No GWAS results: ", paste0(work_dir, method, qc, unrel, "_dsm5diag_", dsm5diag, ".PHENO1.glm.linear"), "\n")
      } else {
        writeLines(paste0("\n\nSignificant SNPs count for dsm5diag ", res$dsm5diag, ": ", res$snp_count, " for GWAS using ", method, qc, unrel), count_conn)
      }
    }
    end_time <- Sys.time()
    total_duration <- difftime(end_time, start_time, units = "mins")
    writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
    writeLines(paste("Total missing files:", missing), readme_conn)
    close(readme_conn)
    close(count_conn)
    cat(paste0("***************Finished running commands:", method, qc, unrel, "***************\n\n\n"))
    #=========================================================================
  }
}