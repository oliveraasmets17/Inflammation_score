# Inflammation score

This repo includes codes to calculate microbiome inflammation index and evaluate its performance for disease classification. Microbiome inflammation index is based on predicting glycoprotein acetylation (GlycA) from the microbiome taxonomic data using LASSO regression. GlycA is an inflammatory marker that can characterize systemic low-grade inflammation. The microbiome inflammation score is defined as the weighted sum of the CLR-abundances of microbial species and aims to provide an estimation of inflammation in the gut microbiome adjusted for age, gender, and BMI. 

Repo contains two R scripts:

1. Src/GlycA_employ_applyModel.R - script to calculate microbiome inflammation score from the given MetaPhlan 3 taxonomic profile

**Input**
- MetaPhlan 3 taxonomic profile with sample id-s as rownames (line 30)
- LASSO model coefficients - *MBscore_EstMB_dataPreprocessingPredictors_fixedCovariates_CLR_LASSO.rds* (line 40; download from the RData folder)
- Data preprocessing steps - *MBscore_EstMB_dataPreprocessingRecipe_fixedCovariates_CLR_LASSO.rds* (line 43; download from the RData folder)

**Output**
- Data frame containing calculated microbiome inflammation (MIS) for each subject
- Output saved into *Results* directory as *Microbiome_inflammation_score.rds* (line 167)
  
2. Src/GlycA_employ_survivalAnalysis.R - script to run survival analysis to evaluate MIS for disease prediction

**Input**
- Output of the *GlycA_employ_applyModel.R* script - recalculated MIS dataset *Results/Microbiome_inflammation_score.rds* (line 21)
- EstMB top decile for MIS - *MBscore_top10percent_EstMB_prediction_fixedCovariates_CLR_LASSO.rds* (line 27; download from the RData folder)

**Output**
- Kaplan-Meier plots for each endpoint colored by top decile vs others (threshold calculated from the input data) (line 119)
- Kaplan-Meier plots for each endpoint colored by top decile vs others (threshold calculated in the EstMB) (line 123)
- Data frame containing the results of Cox regression models for all endpoints (line 174)
- Data frame containing the PH assumption tests for Cox regression models for all endpoints (line 175)
