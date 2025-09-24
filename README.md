# Inflammation score

This repo includes codes to calculate microbiome inflammation index and evaluate its performance for disease classification. Microbiome inflammation index is based on predicting glycoprotein acetylation (GlycA) from the microbiome taxonomic data using LASSO regression. GlycA is an inflammatory marker that can characterize systemic low-grade inflammation. The microbiome inflammation score is defined as the weighted sum of the PA-transformed abundances of microbial species and aims to provide an estimation of inflammation in the gut microbiome adjusted for age, gender, and BMI. 

Repo contains two R scripts:

## Src/GlycA_employ_applyModel.R - script to calculate microbiome inflammation score from the given MetaPhlan 3 taxonomic profile

**Input**
- MetaPhlan 3 OR MetaPhlan 4 taxonomic profiles with sample id-s as rownames (line 24)
- LASSO model coefficients - *Results/applyMIS_predictors_EstMB_PA_LASSOF.rds*/*Results/applyMIS_predictors_EstMB_PA4_LASSOF.rds* (download from the Results folder)
- Data preprocessing steps - *Results/applyMIS_predictors_EstMB_PA_LASSOF.rds*/*Results/applyMIS_predictors_EstMB_PA4_LASSOF.rds* (download from the Results folder)

**Output**
- Data frame containing calculated microbiome inflammation (MIS) for each subject
- Output saved into *Results* directory as *Microbiome_inflammation_score_*.rds* (line 167)
  
## Src/GlycA_employ_survivalAnalysis.R - script to run survival analysis to evaluate MIS for disease prediction

**Input**
- Output of the *GlycA_employ_applyModel.R* script - recalculated MIS dataset *Results/Microbiome_inflammation_score_*.rds* (line 21)
- EstMB top decile for MIS - *Results/MIS_EstMB_top20percentile_PA_LASSOF.rds*/*Results/MIS_EstMB_top20percentile_PA4_LASSOF.rds* (line 27; download from the RData folder)

**Output**
- Kaplan-Meier plots for each endpoint colored by top 20% vs others (threshold calculated from the input data) (line 167)
- Kaplan-Meier plots for each endpoint colored by top 20% vs others (threshold calculated in the EstMB) (line 175)
- Data frame containing the results of Cox regression models for all endpoints (line 239)
- Data frame containing the PH assumption tests for Cox regression models for all endpoints (line 240)
