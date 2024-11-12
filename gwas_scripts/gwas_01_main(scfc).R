#

# Running GWAS 

# 1. Script for running GWAS for the main analysis

############################################################
## Readme: this code runs and saves GWAS results for SCFC
############################################################


## Utilities
library(neuroCombat)
library(tidyr)
library(tibble)
library(dplyr)
library(readxl)
library(ggplot2)
library(readr)
path <- "utilities/"
load(paste0(path, "coupling_cbcl.RData"))
load(paste0(path, "scfc_coupling_matched.RData"))
new_fam <- read.table(paste0(path, "acspsw03.txt"), header = T)
new_fam <- new_fam %>%
  select(rel_family_id, subjectkey, eventname, race_ethnicity)
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
### Fill NAs in race and ethnicity with mode
calculate_mode <- function(x) {
  ux <- unique(na.omit(x))  # Remove NA values before finding the mode
  ux[which.max(tabulate(match(x, ux)))]
}
CFDRs$race <- ifelse(is.na(CFDRs$race), calculate_mode(CFDRs$race), CFDRs$race)
CFDRs$ethn <- ifelse(is.na(CFDRs$ethn), calculate_mode(CFDRs$ethn), CFDRs$ethn)
### write new id list to local
fid_iid_scfc <- CFDRs %>%
  select(fid, iid)
write.table(fid_iid_scfc, paste0("bfile/fid_iid_scfc.txt"), col.names = F, row.names = F, quote = F, sep = "\t")


## Read Phenotype
SCFC_coup_matrix <- as.data.frame(SCFC_coup %>% unlist() %>% matrix(ncol = 379, byrow = T))
SCFC_coup_matrix[is.na(SCFC_coup_matrix)] = 0 # 1. Fill NAs with 0
colnames(SCFC_coup_matrix) <- mmp_subcor
iid <- c() # Get family id and id; format by adding "-" to scfc ids
### Loop through each ID in the list
for (id in sub_SC_FC) {
    # print(id)
    id_format <- paste0(substring(id, 1, 4), "_", substring(id, 5))
    iid <- append(iid, id_format) # Append modified ID to the list
}
SCFC_DF <- data.frame(iid)
SCFC_DF <- cbind(SCFC_DF, SCFC_coup_matrix)
ALL_DF <- CFDRs %>%
  inner_join(SCFC_DF, by = "iid") # This is the big big DF that has everything

### Do Fisher Z Transformation
# use fisher transformation
fisher_transform <- function(r, n) {
  z <- 0.5 * log((1 + r) / (1 - r))
  z_scaled <- z * sqrt(n - 3)
  return(z_scaled)
}
SCFC_coup_matrix_transformed <- apply(SCFC_coup_matrix, 2, function(x) fisher_transform(x, n = 379))

### HARMONIZATION 
scfc_matrix <- t(as.matrix(SCFC_coup_matrix_transformed))
#### get covariates matrix
interview_age = ALL_DF$interview_age
sex = ALL_DF$sex
race = ALL_DF$race
ethn = ALL_DF$ethn
#### get batch
batch <- ALL_DF$site_id_l
mod <- model.matrix(~interview_age+sex+race+ethn)

#_______________________________________________________________________
#### harm2 is harmonized for age, sex race and ethnicity information
data_trans_harm2 <- neuroCombat(dat=scfc_matrix, batch=batch, mod=mod)
scfc_trans_harm2 <- cbind(fid_iid_scfc, as.data.frame(t(data_trans_harm2$dat.combat)))
write.csv(scfc_trans_harm2, "utilities/scfc_harm2_sd.csv", row.names = FALSE)
for (i in 1:379) {
  phenotype_data <- scfc_trans_harm2[, c(1, 2, i + 2)]
  file_name <- paste0("phenotype_roi/scfc_sd_harm2_family_roi", i, ".phen")
  write.table(phenotype_data, file = file_name, row.names = F, col.names = F, quote = F)
  # cat("File", file_name, "has been created.\n")
}

### Read SNP bfile
old_fam <- read.table("bfile/ABCD_all_old.fam")
colnames(old_fam) <- c("fid", "iid", "pid", "mid", "sex", "phenotype")
new_fam <- old_fam %>%
  left_join(fid_iid_scfc, by = "iid") %>%
  mutate(fid = ifelse(is.na(fid.y), fid.x, fid.y)) %>%
  select(fid, iid, pid, mid, sex, phenotype)

write.table(new_fam, paste0("bfile/ABCD_all.fam"), row.names = F, col.names = F, quote = F, sep = "\t")
### Subset SNP bfile and qc
dirpath <- '/Users/kurt/Documents/Data_Science_Project/0928_GWAS_standardized/'
plink_path <- "/Users/kurt/Documents/Data_Science_Project/0928_GWAS_standardized/"
#### Subset by id
system(paste0(plink_path, "plink --bfile ", dirpath, "bfile/ABCD_all --keep ", dirpath, "bfile/fid_iid_scfc.txt --make-bed --out ", dirpath, "bfile/scfc_all"), intern = FALSE)
#### QC
system(paste0(plink_path, "plink --bfile ", dirpath, "bfile/scfc_all --geno 0.1 --maf 0.05 --mind 0.1 --make-bed --out ", dirpath, "bfile/scfc_all_qc"), intern = FALSE)


## PLINK GWAS using --assoc ()
sample <- 379
library(stringr)
files <- list.files(path = paste0(dirpath, 'phenotype_roi'))
phens <- files[stringr::str_detect(files, "\\.phen$")]
file_pairs <- setNames(paste0('gwasplink19_', sub("\\.phen$", "", phens)), phens)

lib_path <- 'export DYLD_LIBRARY_PATH="/Users/kurt/Documents/Data_Science_Project/gcta-1.94.2-MacOS-ARM-x86_64/gcclib:$DYLD_LIBRARY_PATH"'
system(lib_path, intern = TRUE)  # Set library path

readme_path <- paste0(dirpath, "gwas_results/gwasplink19/readme.txt")
readme_conn <- file(readme_path, "w")
writeLines("Commands run for GWAS --assoc analysis:", readme_conn)
start_time <- Sys.time()

for (idx in seq_along(file_pairs)) {
  if (idx <= sample) {
    phenotype_file <- names(file_pairs)[idx]
    out_file <- file_pairs[[idx]]
    command <- paste0(dirpath,
                      "plink ",
                      "--allow-no-sex ",
                      "--assoc ",
                      "--bfile ", dirpath, "bfile/scfc_all_qc ",
                      "--pheno ", dirpath, "phenotype_roi/", phenotype_file, " ",
                      "--out ", dirpath, "gwas_results/gwasplink19/", out_file
               )
    system(command, intern = FALSE)
    writeLines(command, readme_conn, useBytes = TRUE)
    
    cat("\n\nGWAS(PLINK1.9) command executed successfully for the phenotype file:", phenotype_file, "\n\n\n")
  }
}
end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")

writeLines(sprintf("Total time to run all commands: %s minutes", total_duration), readme_conn)
writeLines("End of command list.", readme_conn)
close(readme_conn)


## Read GWAS data
readme_path <- paste0(dirpath, "gwas_results/gwasplink19/gwasplink19_csv/readme.txt") # create readme.txt
readme_conn <- file(readme_path, "w")
writeLines("GWAS using PLNK1.9\n no covariates\n\n\n", readme_conn)
writeLines(paste0("Commands run for GWAS analysis: ", command, "\n\n"), readme_conn)

start_time <- Sys.time()
missing <- 0
for (i in 1:sample) {
  filename <- paste0(dirpath, "gwas_results/gwasplink19/gwasplink19_scfc_sd_harm2_family_roi", i, ".qassoc")
  
  if (file.exists(filename)) {
    qassoc_data <- read.table(filename, header = TRUE, sep = "")
    output_filename <- paste0(dirpath, "gwas_results/gwasplink19/gwasplink19_csv/gwasplink19_roi", i, ".csv")
    write.csv(qassoc_data, output_filename, row.names = FALSE, quote = FALSE)
  } else {
    missing = missing + 1
    writeLines(paste0("Missing file for ROI ", i, ": ", filename), readme_conn)
    cat("File does not exist: ", filename, "\n")  # Print message to the console as well
  }
}

end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")

writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
writeLines(paste("Total missing files:", missing), readme_conn)
close(readme_conn)


## Additional: check results
results_path <- "/Users/kurt/Documents/Data_Science_Project/0928_GWAS_standardized/gwas_results/gwasplink19/"
count_path <- paste0(results_path, "SNP_counts.txt")
count_conn <- file(count_path, "w")
writeLines("PLINK1.9 --assoc SNP counts under 0.05/379\n\n", readme_conn)

for (i in 1:sample) {
  gwas_assoc <- read.table(paste0(results_path, "gwasplink19_scfc_sd_harm2_family_roi", i, ".qassoc"), header = T) %>%
    filter(P < 0.00013192612)
  writeLines(paste0("\n\nSignificant SNPs count for ROI ", i, ": ", nrow(gwas_assoc)), count_conn)
}
close(count_conn)


## Additional: PLINK GWAS using --linear (same as --assoc for continuous phenotypes)
sample <- 379

library(stringr)
files <- list.files(path = paste0(dirpath, 'phenotype_roi'))
phens <- files[stringr::str_detect(files, "\\.phen$")]
file_pairs <- setNames(paste0('gwasplink19_linear_', sub("\\.phen$", "", phens)), phens)

lib_path <- 'export DYLD_LIBRARY_PATH="/Users/kurt/Documents/Data_Science_Project/gcta-1.94.2-MacOS-ARM-x86_64/gcclib:$DYLD_LIBRARY_PATH"'
system(lib_path, intern = TRUE)  # Set library path

readme_path <- paste0(dirpath, "gwas_results/gwasplink19_linear/readme.txt")
readme_conn <- file(readme_path, "w")
writeLines("Commands run for GWAS --linear analysis:", readme_conn)
start_time <- Sys.time()

for (idx in seq_along(file_pairs)) {
  if (idx <= sample) {
    phenotype_file <- names(file_pairs)[idx]
    out_file <- file_pairs[[idx]]
    
    command <- paste0(dirpath,
                      "plink ",
                      "--allow-no-sex ",
                      "--linear ",
                      "--bfile ", dirpath, "bfile/scfc_all_qc ",
                      "--pheno ", dirpath, "phenotype_roi/", phenotype_file, " ",
                      "--out ", dirpath, "gwas_results/gwasplink19_linear/", out_file
               )
    
    system(command, intern = FALSE)
    writeLines(command, readme_conn, useBytes = TRUE)
    
    cat("\n\nGWAS --linear (PLINK1.9) command executed successfully for the phenotype file:", phenotype_file, "\n\n\n")
  }
}

end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")

writeLines(sprintf("Total time to run all commands: %s minutes", total_duration), readme_conn)
writeLines("End of command list.", readme_conn)
close(readme_conn)


## Read GWAS data
readme_path <- paste0(dirpath, "gwas_results/gwasplink19_linear/gwasplink19_linear_csv/readme.txt")
readme_conn <- file(readme_path, "w")
writeLines("GWAS using PLNK1.9 --linear\n no covariates\n\n\n", readme_conn)
writeLines(paste0("Commands run for GWAS analysis: ", command, "\n\n"), readme_conn)

start_time <- Sys.time()

missing <- 0
for (i in 1:sample) {
  filename <- paste0(dirpath, "gwas_results/gwasplink19_linear/gwasplink19_linear_scfc_sd_harm2_family_roi", i, ".assoc.linear")
  
  if (file.exists(filename)) {
    qassoc_data <- read.table(filename, header = TRUE, sep = "")
    output_filename <- paste0(dirpath, "gwas_results/gwasplink19_linear/gwasplink19_linear_csv/gwasplink19_linear_roi", i, ".csv")
    write.csv(qassoc_data, output_filename, row.names = FALSE, quote = FALSE)
  } else {
    missing = missing + 1
    writeLines(paste0("Missing file for ROI ", i, ": ", filename), readme_conn)
    cat("File does not exist: ", filename, "\n")  # Print message to the console as well
  }
}

end_time <- Sys.time()
total_duration <- difftime(end_time, start_time, units = "mins")

writeLines(sprintf("Total time to process all files: %s minutes", total_duration), readme_conn)
writeLines(paste("Total missing files:", missing), readme_conn)
close(readme_conn)


## Additional: check results
results_path <- "/Users/kurt/Documents/Data_Science_Project/0928_GWAS_standardized/gwas_results/gwasplink19_linear/"
count_path <- paste0(results_path, "SNP_counts.txt")
count_conn <- file(count_path, "w")
writeLines("PLINK1.9 --linear SNP counts under 0.05/379\n\n", readme_conn)
for (i in 1:sample) {
  gwas_assoc <- read.table(paste0(results_path, "gwasplink19_linear_scfc_sd_harm2_family_roi", i, ".assoc.linear"), header = T) %>%
    filter(P < 0.00013192612)
  writeLines(paste0("\n\nSignificant SNPs count for ROI ", i, ": ", nrow(gwas_assoc)), count_conn)
}
close(count_conn)

