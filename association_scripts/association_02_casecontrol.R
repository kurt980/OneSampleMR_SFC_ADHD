#

# Sensitivity Analysis

# 2. Case/Control Analysis ()

## Read Data
library(snpStats)
library(tibble)
library(dplyr)
library(AER)
library(readxl)
library(ggplot2)
library(doParallel)
library(foreach)
library(coin)

### Load Utilities, data, tool and define directory
utilities_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_SCFC_utilities/"
bfile_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/ABCD_bfile/"
cov_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Covariates/"
harm_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Harmonized_scfc/"
pheno_path <- "/Users/grc8mu/Desktop/DS/ABCD_data/Phenotype/"
tool_path <- "/Users/grc8mu/Desktop/DS/Tools/"
directory <- "/Users/grc8mu/Desktop/DS/1029_ADHD_casecontrol/"
load(paste0(utilities_path, "coupling_cbcl.RData"))
load(paste0(utilities_path, "scfc_coupling_matched.RData"))
new_fam <- read.table(paste0(utilities_path, "acspsw03.txt"), header = T)
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
fid_iid_scfc <- CFDRs %>%
  select(fid, iid)

### Make diagnosis to phenotype (ADHD), encoded as: https://www.cog-genomics.org/plink/1.9/formats
#### Load the diagnosis dataset
dsm5_diag <- read_excel(paste0(utilities_path, "abcd_diagnosis_sum_baseline_n11976.xlsx"))
#### Select and rename columns
dsm5_diag <- dsm5_diag %>%
  select(subjectkey, ADHD, Dep, Bip, DMDD, Schi, Anx, OCD, Eat, ODD, Cond, PTSD) %>%
  rename(iid = subjectkey) %>%
  inner_join(fid_iid_scfc, by = "iid") %>%
  select(fid, iid, ADHD, Dep, Bip, DMDD, Schi, Anx, OCD, Eat, ODD, Cond, PTSD)
#### Recode 0 to 1, 1 to 2, and anything else to -9 for all phenotype columns
dsm5_diag <- dsm5_diag %>%
  mutate_at(vars(ADHD:PTSD), ~ ifelse(. == 0, 1, ifelse(. == 1, 2, -9)))
#### Write each phenotype to a separate file
for (i in 1:11) {
  phenotype_data <- dsm5_diag[, c(1, 2, i + 2)]
  file_name <- paste0(pheno_path, "dsm5diag_", colnames(dsm5_diag)[i + 2], ".phen")
  write.table(phenotype_data, file = file_name, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

## Analysis
### Normality Test
#### read diagnosis dataset
dsm5_diag <- read_excel(paste0(utilities_path, "abcd_diagnosis_sum_baseline_n11976.xlsx"))
dsm5_diag <- dsm5_diag %>%
  select(subjectkey, ADHD) %>%
  rename(iid = subjectkey)
value_counts <- dsm5_diag %>%
  group_by(ADHD) %>%
  summarise(count = n())

print("1085/10642 are ADHD positive")

dsm5_diag <- dsm5_diag %>%
  inner_join((demo_cbcl_base_scfc %>%
                select(subjectkey, cbcl_scr_dsm5_adhd_t, site_id_l, race, ethn, interview_age, sex) %>%
                rename(iid = subjectkey) %>%
                mutate(iid = sub("^(.{4})", "\\1_", iid))),
             by = "iid") %>%
  filter(cbcl_scr_dsm5_adhd_t != 0 & cbcl_scr_dsm5_adhd_t != 9 & cbcl_scr_dsm5_adhd_t != -9,
         ADHD %in% c(0, 1)) %>%
  mutate(ADHD = as.factor(ADHD))
value_counts <- dsm5_diag %>%
  group_by(ADHD) %>%
  summarise(count = n())

calculate_mode <- function(x) {
  ux <- unique(na.omit(x))  # Remove NA values before finding the mode
  ux[which.max(tabulate(match(x, ux)))]
}
dsm5_diag$race <- ifelse(is.na(dsm5_diag$race), calculate_mode(dsm5_diag$race), dsm5_diag$race)
dsm5_diag$ethn <- ifelse(is.na(dsm5_diag$ethn), calculate_mode(dsm5_diag$ethn), dsm5_diag$ethn) # Fill NAs in race and ethnicity with mode

print("494/6381 are ADHD positive")

#### Check normality with Shapiro-Wilk test
library(nortest)
ad.test(dsm5_diag$cbcl_scr_dsm5_adhd_t[dsm5_diag$ADHD == 0])
shapiro.test(dsm5_diag$cbcl_scr_dsm5_adhd_t[dsm5_diag$ADHD == 1])

print("Not NORMAL")

#### t-test
t_test_result <- t.test(cbcl_scr_dsm5_adhd_t ~ ADHD, data = dsm5_diag)
print(t_test_result)

#### Wilcox
wilcox_test_result <- wilcox.test(cbcl_scr_dsm5_adhd_t ~ ADHD, data = dsm5_diag)
print(wilcox_test_result)

#### KS
ks_test_result <- ks.test(
  dsm5_diag$cbcl_scr_dsm5_adhd_t[dsm5_diag$ADHD == 0],
  dsm5_diag$cbcl_scr_dsm5_adhd_t[dsm5_diag$ADHD == 1]
)
print(ks_test_result)

#### Permutation
perm_test_result <- wilcox_test(cbcl_scr_dsm5_adhd_t ~ ADHD, data = dsm5_diag, distribution = "approximate")
print(perm_test_result)

#### plot distribution of positive samples vs. negative samples
plot <- ggplot(dsm5_diag, aes(x = cbcl_scr_dsm5_adhd_t, fill = as.factor(ADHD))) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of ADHD Scores by Diagnosis Group",
    x = "ADHD Score",
    y = "Density",
    fill = "ADHD Diagnosis"
  ) +
  scale_fill_manual(values = c("#1E90FF", "#FC4E2A"), labels = c("Non-ADHD", "ADHD")) +
  theme_minimal()
ggsave(paste0(directory, paste0("visualization/adhd_+vs-", ".png")), plot = plot, width = 10, height = 8, dpi = 300)


## Run regression SCFC ~ ADHD diagnosis + covariates
#### Read SCFC values
scfc <- read.csv(paste0(harm_path, "scfc_trans_harm2.csv"))
df <- dsm5_diag %>%
  select(iid, ADHD, site_id_l, race, ethn, interview_age, sex) %>%
  rename(adhd_diag = ADHD) %>%
  inner_join(scfc, by = "iid")

adhd_diag <- as.factor(df$adhd_diag)
ROIs <- as.matrix(df[, 9:387])  # starting from column 4

roi_numbers <- 1:379
roi_df <- data.frame(
  roi = roi_numbers,      # The numeric sequence
  ROI = colnames(ROIs)        # The list of ROI names
)
covariates <- data.frame(
  race = df$race, 
  ethn = df$ethn, 
  interview_age = df$interview_age, 
  sex = df$sex
)
#=============================================================================
### Run linear regression
# Set up parallel backend
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

regression_results <- foreach(roi = 1:ncol(ROIs), .combine = rbind, .packages = "dplyr") %dopar% {
  roi_data <- data.frame(
    adhd_diag = adhd_diag,
    ROI = ROIs[, roi],
    covariates
  )
  
  model <- lm(ROI ~ adhd_diag + race + ethn + interview_age + sex, data = roi_data)
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  data.frame(
    roi = roi,
    ROI = colnames(ROIs)[roi],
    coeff_reg_adhd = coefficients["adhd_diag1", "Estimate"],
    p_value_reg_adhd = coefficients["adhd_diag1", "Pr(>|t|)"],
    coeff_reg_raceAsian = coefficients["raceAsian", "Estimate"],
    p_value_reg_raceAsian = coefficients["raceAsian", "Pr(>|t|)"],
    coeff_reg_raceBlack = coefficients["raceBlack", "Estimate"],
    p_value_reg_raceBlack = coefficients["raceBlack", "Pr(>|t|)"],
    coeff_reg_raceHawaiian = coefficients["raceHawaiian or Pacific Islander", "Estimate"],
    p_value_reg_raceHawaiian = coefficients["raceHawaiian or Pacific Islander", "Pr(>|t|)"],
    coeff_reg_raceWhite = coefficients["raceWhite", "Estimate"],
    p_value_reg_raceWhite = coefficients["raceWhite", "Pr(>|t|)"],
    coeff_reg_raceOther = coefficients["raceOther", "Estimate"],
    p_value_reg_raceOther = coefficients["raceOther", "Pr(>|t|)"],
    coeff_reg_ethnNotHispanic = coefficients["ethnNot Hispanic", "Estimate"],
    p_value_reg_ethnNotHispanic = coefficients["ethnNot Hispanic", "Pr(>|t|)"],
    coeff_reg_interview_age = coefficients["interview_age", "Estimate"],
    p_value_reg_interview_age = coefficients["interview_age", "Pr(>|t|)"],
    coeff_reg_sexM = coefficients["sexM", "Estimate"],
    p_value_reg_sexM = coefficients["sexM", "Pr(>|t|)"]
  )
}

# Stop the parallel backend
stopCluster(cl)

### Merge regression results with ROI information and select all specified columns
final_df <- regression_results %>%
  inner_join(roi_df, by = c("roi", "ROI")) %>%
  select(
    roi, ROI,
    coeff_reg_adhd, p_value_reg_adhd,
    coeff_reg_raceAsian, p_value_reg_raceAsian,
    coeff_reg_raceBlack, p_value_reg_raceBlack,
    coeff_reg_raceHawaiian, p_value_reg_raceHawaiian,
    coeff_reg_raceWhite, p_value_reg_raceWhite,
    coeff_reg_raceOther, p_value_reg_raceOther,
    coeff_reg_ethnNotHispanic, p_value_reg_raceAsian,  # note: ensure p_value_raceAsian is correct here
    coeff_reg_interview_age, p_value_reg_interview_age,
    coeff_reg_sexM, p_value_reg_sexM
  )
output_directory <- paste0(directory, "linear_reg/")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = FALSE)
}
output_file <- paste0(output_directory, "adhd_diag_scfc_reg_results", ".csv")
write.csv(final_df, output_file, row.names = FALSE)
print("Regression results saved successfully!")

significance_level <- 0.05
covariate_p_values <- c(
  "p_value_reg_raceAsian", "p_value_reg_raceBlack", "p_value_reg_raceHawaiian",
  "p_value_reg_raceWhite", "p_value_reg_raceOther", "p_value_reg_ethnNotHispanic",
  "p_value_reg_interview_age", "p_value_reg_sexM"
)

### Analyze Regression results
# 1. Count the number of significant p-values for each covariate without BH correction
significant_counts <- data.frame(
  Covariate = covariate_p_values,
  Significant_Count = sapply(covariate_p_values, function(column) {
    sum(final_df[[column]] < significance_level, na.rm = TRUE)
  })
)
print("Significant Counts Without FDR Correction:")
print(significant_counts)

# 2. Apply BH correction to each covariate's p-values and count significant results
adjusted_counts <- data.frame(
  Covariate = covariate_p_values,
  Significant_Count_BH = sapply(covariate_p_values, function(column) {
    sum(p.adjust(final_df[[column]], method = "BH") < significance_level, na.rm = TRUE)
  })
)
print("Significant Counts With FDR (BH) Correction:")
print(adjusted_counts)

# 3. Check for any significant race level in each ROI before BH correction
significant_race_ROIs <- rowSums(
  final_df[, c("p_value_reg_raceAsian", "p_value_reg_raceBlack", "p_value_reg_raceHawaiian", 
               "p_value_reg_raceWhite", "p_value_reg_raceOther")] < significance_level, 
  na.rm = TRUE
) > 0
num_significant_race_ROIs <- sum(significant_race_ROIs)
print("Number of ROIs with Significant Race Effect (No FDR Correction):")
print(num_significant_race_ROIs)

# 4. Apply FDR (BH) correction specifically to race levels and count significant results
adjusted_p_values_race <- final_df %>%
  mutate(
    adj_p_value_raceAsian = p.adjust(p_value_reg_raceAsian, method = "BH"),
    adj_p_value_raceBlack = p.adjust(p_value_reg_raceBlack, method = "BH"),
    adj_p_value_raceHawaiian = p.adjust(p_value_reg_raceHawaiian, method = "BH"),
    adj_p_value_raceWhite = p.adjust(p_value_reg_raceWhite, method = "BH"),
    adj_p_value_raceOther = p.adjust(p_value_reg_raceOther, method = "BH")
  )

significant_race_ROIs_FDR <- rowSums(
  adjusted_p_values_race[, c("adj_p_value_raceAsian", "adj_p_value_raceBlack", "adj_p_value_raceHawaiian", 
                             "adj_p_value_raceWhite", "adj_p_value_raceOther")] < significance_level,
  na.rm = TRUE
) > 0
num_significant_race_ROIs_FDR <- sum(significant_race_ROIs_FDR)
print("Number of ROIs with Significant Race Effect (With FDR Correction):")
print(num_significant_race_ROIs_FDR)



## Do a reverse logistic regression (result is similar)
# Set up parallel backend
num_cores <- detectCores() - 1  # Use one less than total cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

logregression_results <- foreach(roi = 1:ncol(ROIs), .combine = rbind, .packages = "dplyr") %dopar% {
  # Prepare data for the current ROI
  roi_data <- data.frame(
    adhd_diag = adhd_diag,
    ROI = ROIs[, roi],
    covariates
  )
  # Run the logistic regression
  model <- glm(adhd_diag ~ ROI + race + ethn + interview_age + sex, data = roi_data, family = binomial)
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  # Extract coefficients and p-values for SCFC and covariates
  data.frame(
    roi = roi,
    ROI = colnames(ROIs)[roi],
    coeff_reg_ROI = coefficients["ROI", "Estimate"],
    p_value_reg_ROI = coefficients["ROI", "Pr(>|z|)"],
    coeff_reg_raceAsian = coefficients["raceAsian", "Estimate"],
    p_value_reg_raceAsian = coefficients["raceAsian", "Pr(>|z|)"],
    coeff_reg_interview_age = coefficients["interview_age", "Estimate"],
    p_value_reg_interview_age = coefficients["interview_age", "Pr(>|z|)"],
    coeff_reg_sexM = coefficients["sexM", "Estimate"],
    p_value_reg_sexM = coefficients["sexM", "Pr(>|z|)"]
  )
}

# Stop the parallel backend
stopCluster(cl)


final_df_log <- logregression_results %>%
  inner_join(roi_df, by = c("roi", "ROI")) %>%
  select(roi, ROI, coeff_reg_ROI, p_value_reg_ROI, coeff_reg_raceAsian, p_value_reg_raceAsian, coeff_reg_interview_age, p_value_reg_interview_age, coeff_reg_sexM, p_value_reg_sexM)
output_directory <- paste0(directory, "logistic_reg/")
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = FALSE)
}
output_file <- paste0(output_directory, "adhd_diag_scfc_reg_results", ".csv")
write.csv(final_df_log, output_file, row.names = FALSE)
print("Logistic regression results saved successfully!")

roi_data <- data.frame(
  adhd_diag = adhd_diag,
  ROI = ROIs[, 1],
  covariates
)

# Run the logistic regression
model <- glm(adhd_diag ~ ROI + race + ethn + interview_age + sex, data = roi_data, family = binomial)
summary_model <- summary(model)
coefficients <- summary_model$coefficients

str(roi_data)