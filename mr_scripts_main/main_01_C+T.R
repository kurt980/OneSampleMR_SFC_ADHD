#

# Mendelian Randomization Analysis (Main Analysis Forward MR)

# 1. Instrumental Variable Selection: Clumping and Thresholding (C+T)

############################################################
## Readme: this code is for running C+T for IV selection.
############################################################


library(optparse)
library(parallel)
library(readr)
library(bigreadr)
library(bigsnpr)
library(data.table)
library(dplyr)
library(metafor)

### Define directories
in_dir <- "/Users/kurt/Documents/Data_Science_Project/0928_GWASplink19_standardized/"
directory <- "/Users/kurt/Documents/Data_Science_Project/1004_MR_GWASplink19/"
plink_path <- paste0("/Users/kurt/Documents/Data_Science_Project/Tools/", "plink")
bfile_path <- paste0("/Users/kurt/Documents/Data_Science_Project/ABCD_data/", "ABCD_bfile/scfc_all_qc")

#-------------1st level directories---------------------------------------------
# 1.1: for storing bfile by chr, same for all pvalthr
if (!dir.exists(paste0(directory, "bfile_by_chr/"))) { dir.create(paste0(directory, "bfile_by_chr/"), recursive = F) }
# 1.2: for storing C+T results (contain parallel Thresholding, Clumping and C+T directories)
ct_results_path <- paste0(directory, "C+T/") 
if (!dir.exists(ct_results_path)) { dir.create(ct_results_path, recursive = F) }
# 1.3: beneath C+T/: Thresholding, Clumping and C+T
thresh_path <- paste0(directory, "C+T/Thresholding/") # 1.3.1: Thresholding results
if (!dir.exists(thresh_path)) { dir.create(thresh_path, recursive = F) }
clumping_path <- paste0(directory, "C+T/Clumping/")   # 1.3.2: Clumping results
if (!dir.exists(clumping_path)) { dir.create(clumping_path, recursive = F) }
ct_path <- paste0(directory, "C+T/C+Ted/")            # 1.3.3: Processed C+T results containing SNP and test results
if (!dir.exists(ct_path)) { dir.create(ct_path, recursive = F) }

## I. Break geno files by CHRs (ABCD actually comes separated but QC combined them) and store to "bfile_by_chr/"
for (chr in 1:22) {
  output_path <- paste0('"', directory, "bfile_by_chr/scfc_all_qc_chr", chr, '"')
  command <- paste(plink_path, "--bfile", bfile_path, "--chr", chr, "--make-bed --out", output_path)
  system(command, intern=TRUE)
}


#========================================================================================================================  
#========================================================================================================================  

################################################### C + T ###################################################
for (roi in 1:379) {
  #==============================read GWAS results==============================
  # Either a fastGWA file, a PLINK2.0, or a PLINK1.9 file
  fastGWA_path <- paste0(in_dir, "gwas_results/fastGWA_default_allcovar/fastGWA_default_allcovar_csv/fastGWA_roi", roi, ".csv")
  gwas20_path <- paste0(in_dir, "gwas_unrelated/gwas_results/gwas_unrelated/gwas_unrelated_csv/gwas_scfc_unrelated_sd_harm2_family_roi", roi, ".csv")
  gwas19_path <- paste0(in_dir, "gwas_results/gwasplink19/gwasplink19_csv/gwasplink19_roi", roi, ".csv")
  #-----------------------------------------------------------------------------
  if (file.exists(fastGWA_path)) {
    print(paste0("roi",roi,":fastGWA"))
    gwas_roi <- read.csv(fastGWA_path, header = TRUE) %>%
      rename(A1_fastgwa = A1, A2_fastgwa = A2, AF1_fastgwa = AF1) %>%
      mutate(REF_plink = NA, ALT_plink = NA, PROVISIONAL_REF_plink = NA, A1_plink = NA, OMITTED_plink = NA, AF1_plink = NA) %>% # Add empty columns for PLINK-specific fields
      select(CHR, SNP, POS, REF_plink, ALT_plink, PROVISIONAL_REF_plink, A1_plink, OMITTED_plink, AF1_plink, A1_fastgwa, A2_fastgwa, AF1_fastgwa, N, BETA, SE, P) # only keep necessary columns (both plink and fastgwa columns are kept)
    
  } else if (file.exists(gwas20_path)) {
    print(paste0("roi",roi,":GWAS plink2.0"))
    gwas_roi <- read.csv(gwas20_path, header = TRUE) %>%
      rename(CHR = `#CHROM`, REF_plink = REF, ALT_plink = ALT, PROVISIONAL_REF_plink = `PROVISIONAL_REF?`, A1_plink = A1, OMITTED_plink = OMITTED, AF1_plink = A1_FREQ, N = OBS_CT) %>%
      mutate(A1_fastgwa = NA, A2_fastgwa = NA, AF1_fastgwa = NA) %>% # Add empty columns for FastGWA-specific fields
      select(CHR, SNP, POS, REF_plink, ALT_plink, PROVISIONAL_REF_plink, A1_plink, OMITTED_plink, AF1_plink, A1_fastgwa, A2_fastgwa, AF1_fastgwa, N, BETA, SE, P) # only keep necessary columns (both plink and fastgwa columns are kept)
    
  } else if (file.exists(gwas19_path)) {
    print(paste0("roi",roi,":GWAS plink1.9"))
    gwas_roi <- read.csv(gwas19_path, header = TRUE) %>%
      rename(POS = BP, N = NMISS) %>%
      select(CHR, SNP, POS, N, BETA, SE, P)
    
  } else {
    cat("Neither FastGWA nor GWAS result files are available for ROI", roi, "\n")
    next  # Skip to the next iteration if no files are found
  }
  #-----------------------------------------------------------------------------
  #### IV Selection (using exposure: scfc)
  ### Set Parameters
  # Input:
  NCORES = 1 # Number of cores used for training lassosum2 models under different tuning parameter settings.
  # Parameters for the Clumping step
  kb = 500 # SNPs within 500kb of the index SNP are considered for clumping
  p.ldclump = 0.5 # LD p-value threshold for a SNP to be included as an index SNP
  pvalthrs = c(5e-06, 5e-07, 5e-08) # p-value threshold for thresholding
  r2 = 0.001 # SNPs having r2 higher than 0.1 with the index SNPs will be removed # optional: c(0.1,0.2,0.4,0.6,0.8)
  
  trait_name <- "adhd" ## 
  method <- 'method' ## 
  ite <- 1 ## 
  
  ### Run T+C (Thresholding first and Clumping next)============================
  for (pvalthr in pvalthrs) {
    # -----for inspection purpose, output thresholding results to local---------
    gwas_roi.thresh <- gwas_roi[gwas_roi$P <= pvalthr,]
    thresh_count_path <- paste0(thresh_path, "threshold_count/gwasplink19_qc_", pvalthr, "/")
    if (!dir.exists(thresh_count_path)) { dir.create(thresh_count_path, recursive = TRUE) }
    write.table(gwas_roi.thresh, paste0(thresh_count_path, "gwas_thresh_roi",roi,".txt"), row.names = F, col.names = T, quote = F)
    readme_path <- paste0(thresh_count_path, "readme.txt")
    readme_conn <- file(readme_path, "a")
    writeLines(paste0("roi", roi, " after thresholding: ", nrow(gwas_roi.thresh)), readme_conn)
    close(readme_conn)
    
    ## ---------------------2nd level Directories-------------------------------
    ## 2.1 Thresholding
    thresh_idlist_path <- paste0(thresh_path, "thr_idlists/")                   # 2.1.1 for storing SNP id lists for plink to subset bfiles
    if (!dir.exists(thresh_idlist_path)) { dir.create(thresh_idlist_path, recursive = F) }
    thresh_bfile_by_chr_path <- paste0(thresh_path, "thr_bfile_by_chr/")        # 2.1.2 for storing bfile by chr AFTER thresholding (for clumping)
    if (!dir.exists(thresh_bfile_by_chr_path)) { dir.create(thresh_bfile_by_chr_path, recursive = F) }
    thresh_gwas_path <- paste0(thresh_path, "thr_gwas/")                           # 2.1.3 for storing GWAS results filtered AFTER thresholding (for clumping) 
    if (!dir.exists(thresh_gwas_path)) { dir.create(thresh_gwas_path, recursive = F) }
    
    ## 2.2 Clumping
    clump_clumped_path <- paste0(clumping_path, "clumped/")                     # 2.2.1 for storing raw clumped results, by chr
    if (!dir.exists(clump_clumped_path)) { dir.create(clump_clumped_path, recursive = F) }
    
    ## 2.3 C+T does not need additional sub-directories
    
    ## --------------------3rd level Directories (pvalthr)----------------------
    # create directories for thresholding p-values (pvalthr)
    # 2.1.1 thresholding directories
    thresh_idlist_pvalthr_path <- paste0(thresh_idlist_path, "idlist_p=", pvalthr, "/")                    # 2.1.1.1 for storing SNP id lists for current pvalthr
    if (!dir.exists(thresh_idlist_pvalthr_path)) { dir.create(thresh_idlist_pvalthr_path, recursive = F) }
    thresh_bfile_by_chr_pvalthr_path <- paste0(thresh_bfile_by_chr_path, "bfile_by_chr_p=", pvalthr, "/")  # 2.1.2.1 for storing bfiles for current pvalthr
    if (!dir.exists(thresh_bfile_by_chr_pvalthr_path)) { dir.create(thresh_bfile_by_chr_pvalthr_path, recursive = F) }
    thresh_gwas_pvalthr_path <- paste0(thresh_gwas_path, "gwas_p=", pvalthr, "/")                          # 2.1.3.1 for storing GWAS results for current pvalthr
    if (!dir.exists(thresh_gwas_pvalthr_path)) { dir.create(thresh_gwas_pvalthr_path, recursive = F) }
    
    # 2.2.1 clumping directories
    clump_clumped_pvalthr_path <- paste0(clump_clumped_path, "clumped_p=", pvalthr, "/")                 # 2.2.1.1 for storing raw clumped results for current pvalthr
    if (!dir.exists(clump_clumped_pvalthr_path)) { dir.create(clump_clumped_pvalthr_path, recursive = F) }
    
    # 2.3.1 C+T directories
    ct_pvalthr_path <- paste0(ct_path, "p=", pvalthr, "/")                                                 # 2.3.1.1 for storing C+T results for current pvalthr 
    if (!dir.exists(ct_pvalthr_path)) { dir.create(ct_pvalthr_path, recursive = F) }
    ct_final_path <- paste0(ct_pvalthr_path, 'roi', roi, '.', "_p=", pvalthr, "_r2=", r2, ".txt")          # 2.3.1.2 txt file for storing final C+T SNPs, used if clumping is successful
    #---------------------------------------------------------------------------
    for (chr in 1:22){
      
      sum.chrspecific = gwas_roi[gwas_roi$CHR == chr,]
      # --------------------- Step 2.2: Run C+T --------------------------------
      set.seed(2023)
      
      ## --------------------- Thresholding (T) --------------------------------
      # filter SNPs by p-value
      sum.data <- sum.chrspecific
      sum.data <- sum.data[sum.data$P <= pvalthr,]
      # writes list of filtered SNP id
      snp_id_list <- paste0(thresh_idlist_pvalthr_path, "idlist_roi",roi,'_chr',chr,'_pvalthr',pvalthr,".txt") #这个filtered指的是thresholding
      write.table(sum.data[,'SNP'], snp_id_list, quote = F, col.names = F, row.names = F)
      
      # write filtered SNPs bfile back to SNP list using PLINK1.9
      system(paste0(plink_path, " --bfile ",directory, "bfile_by_chr/scfc_all_qc_chr",chr,
                    " --extract ",snp_id_list,
                    " --make-bed --out ",thresh_bfile_by_chr_pvalthr_path,"scfc_all_qc_thr","_roi",roi,'_chr',chr,'_pvalthr',pvalthr), # filtered snps by chr
             intern = T)
      
      # write the threshold GWAS summary with SNP and P-value to local (to be used by clumping)
      sumstats_input <- sum.data[,c('SNP','P')]
      fwrite(sumstats_input, file=paste0(thresh_gwas_pvalthr_path,'gwas_roi',roi,'_chr',chr,'_pvalthr',pvalthr,'.txt'),row.names = F, quote = F, sep=' ')
      
      cat(paste0('\nTHRESHOLDING FINISHED roi: ',roi,' chr: ',chr," pvalthr: ",pvalthr))
      
      ## --------------------- LD Clumping (C) ---------------------
      # outputfilename
      clumpedfilename <- paste0("roi",roi,'_chr',chr,"_clumped") # C+T files naming: removed "adhd" in the name, removed "ite", put "roi" before "chr"
      # Change the directory to the one where PLINK1.9 was stored: (cannot do this in Windows)
      ldclumpcode <- paste0(plink_path, 
                            ' --bfile ',thresh_bfile_by_chr_pvalthr_path,"scfc_all_qc_thr","_roi",roi,'_chr',chr,'_pvalthr',pvalthr,
                            ' --clump ',thresh_gwas_pvalthr_path,'gwas_roi',roi,'_chr',chr,'_pvalthr',pvalthr,'.txt',
                            ' --clump-p1 ',p.ldclump,
                            ' --clump-r2 ',r2,
                            ' --clump-kb ',kb,
                            ' --threads 1',
                            # ' --silent',
                            ' --out ',clump_clumped_pvalthr_path,clumpedfilename)
      system(ldclumpcode, intern = T)
      cat(paste0('\nCLUMPING FINISHED roi: ',roi,' chr: ',chr," pvalthr: ",pvalthr))
      
      ##===================== read clumped files================================
      # First confirm that there is valid clumped results
      if (!file.exists(paste0(clump_clumped_pvalthr_path,clumpedfilename,'.clumped'))) {
        cat(paste0('/nNO CLUMPING RESULTS AVAILABLE for roi: ',roi,' chr: ',chr," pvalthr: ",pvalthr))
        next
      }
      LD <- bigreadr::fread2(paste0(clump_clumped_pvalthr_path,clumpedfilename,'.clumped'))
      clump.snp <- LD[,3,drop=F]
      prs.all <- left_join(clump.snp,sum.data,by="SNP")
      
      #prs.file = prs.all[prs.all$P <= pvalthr,] # USELESS?
      # -------------------------- Append the SNPs by chr ----------------------
      if (file.exists(ct_final_path)) {
        write.table(prs.file, file = ct_final_path, append = TRUE, col.names = FALSE, row.names = F, quote = F)
      } else {
        write.table(prs.file, file = ct_final_path, append = FALSE, col.names = TRUE, row.names = F, quote = F)
      }
      print(paste0('/nC+T DONE for roi: ',roi,' chr: ',chr," pvalthr: ",pvalthr))
    }
    # ------for inspection purpose, output C+T counts results to local----------
    if (file.exists(ct_final_path)) {
      ct.count <- read.table(ct_final_path, header = T)
      readme_path <- paste0(ct_pvalthr_path, "readme.txt")
      readme_conn <- file(readme_path, "a")
      writeLines(paste0("roi", roi, " after C+T: ", nrow(ct.count)), readme_conn)
    } else {
      writeLines(paste0("roi", roi, " after C+T: ", "No SNPs pass C+T"), readme_conn)
    }
    close(readme_conn)
  }
  ####
}
################################################### C + T Done ###################################################
