# One-sample Mendelian Randomization for SCFC on ADHD

This Github includes a pipeline to conduct a MR study for the article "Article" published in Journal.

# Citations

If you consider using code from this Github please cite:

# Data

All data comes from the ABCD Study (Volkow et al., 2018). The ABCD Study is the largest longitudinal study of the developmental trajectories 
of the brain and cognitive, social, and emotional growth of children in the United States. 
Among the 11,099 ABCD participants aged 9-10, we include the 6,381 (3,172 female) samples that have neuroimaging, genetic, and behavioral data collected at baseline in our study.
ABCD study collects data from multiple sites. The Data used for the 6,381 samples in the study contains 22 sites.
Neuroimaging data contains numeric SCFC values for 379 brain ROIs for 6,381 samples acquired from MRI; 
Genetic data contains 1.8 million SNPs per individual, stored as PLINK binary files (ped, map, bed, bim, fam); 
Behavioral data contains ADHD scores for 6,381 samples measured using CBCL standard, and binary ADHD diagnosis (0,1) for 6,381 samples.
Demographic data contains values for sex, race, ethnicity and age.

# Code Pipeline

The code is divided into 3 main directories: association analysis, gwas and mr analysis

- **/gwas_scripts** contains R scripts for running GWAS with different settings and parameters.
  - software included in the scripts are **PLINK1.9**(for initial main analysis), **PLINK2.0**(for GWAS with covariates), **GCTA**(for fastGWA, not included in the final study)
  - gwas scripts contain code for running GWAS on combinations of parameters: QCed/non-QCed data, and on all samples/unrelated samples\
- **/mr_scripts** contains R scripts for conducting MR study: \
  - **C+T**: C+T scripts are for running the C+T IV selection process using GWAS results; this include both the main analysis, the reverse MR analysis, the case/control analysis and sensitivity analysis
    - the clumping is done using **PLINK1.9**\
  - **2sls**: 2sls scripts are for running the 2-stage least square regression to conduct one-sample MR analysis for continuous outcome (ADHD score); this include the process for main analysis, reverse MR, and sensitivity analysis
  - **2sri**: 2sri script is for running the 2-stage residual inclusion estimation to conduct one-sample MR analysis for binary outcome (ADHD diagnosis); this include the process for the case/control MR study
- **/association_scripts** contains R scripts for running correlation analysis, and regression analysis\
- **/analytic_scripts** contains scripts for analyzing, visualizing and summarizing results from the above analysis\
- **/data_processing** contains scripts for pre-processing ABCD data and for creating **PLINK** phenotype files\\\\\\\\
  - pre-processing includes fisher-z transformation and harmonization using neuroCombat (Fortin et al., 2018)

