

# Program setup ----
#------------------------------------------------#
#                                                #
#                 PROGRAM SETUP                  # 
#                                                #
#------------------------------------------------#

# Load packages
#-----------------------------#
library("dplyr")
library("curatedMetagenomicData")
library("survival")
library("ggsurvfit")



# Read calculated microbiome score
#-----------------------------#
MIS_result_df <- readRDS(file = "Results/Microbiome_inflammation_score.rds")

# Local top 10%
local_top10 <- quantile(MIS_result_df$MIS, probs = 0.9)

# EstMB top 10% 
EstMB_top10 <- readRDS("RData/MBscore_top10percent_EstMB_prediction_fixedCovariates_CLR_LASSO.rds")



# Read phenotype data
#-----------------------------#

#NB: REPLACE!
# Example: generate a dataset from the curatedMetagenomicData package - imitate survival data
# TODO: include multiple phenotypes in a long data format; include relevant and available covariates; exclude prevalent cased; apply other exlusion criteria where needed
set.seed(0)
phenotype_data <- curatedMetagenomicData::sampleMetadata %>% 
  dplyr::filter(study_name == "KarlssonFH_2013") %>% 
  dplyr::mutate(# Censoring
                status = ifelse(study_condition == "T2D", 1, 0),
                # Generate random time variable
                time = floor(runif(nrow(.), min = 1, max = 1200)),
                endpoint = "T2D", 
                # Generate covariates
                age = rnorm(nrow(.), mean = 50, sd = 4),
                gender = rbinom(n = nrow(.), size = 1, prob = 0.6)) %>% 
  dplyr::select(sample_id, BMI, age, gender, status, time, endpoint)
  
# Example of the phenotype data
# Status =  0 for censored observation, 1 for event
head(phenotype_data)
# sample_id  BMI      age gender status time endpoint
# 1      S112 24.9 47.45783      1      0  885      T2D
# 2      S118 32.7 48.28008      0      0  418      T2D
# 3      S121 29.7 49.32273      0      0 1138      T2D
# 4      S126 18.6 52.44887      1      0  776      T2D
# 5      S127 26.6 52.71336      1      0   43      T2D
# 6      S131 26.5 52.27181      1      0  716      T2D




# Run basic survival analysis ----
#------------------------------------------------#
#                                                #
#             RUN SURVIVAL ANALYSIS              # 
#                                                #
#------------------------------------------------#

# STEP 1 - merge data
#-----------------------------#
survival_data <- phenotype_data %>% 
  dplyr::left_join(MIS_result_df, by = "sample_id") %>% 
  # Define binary outcomes - top decile vs others + two thresholds (one calculated locally, other calulated in th EstMB cohort)
  dplyr::mutate(MIS_local_top10 = ifelse(MIS >= local_top10, "Top 10%", "Bottom 90%"),
                MIS_EstMB_top10 = ifelse(MIS >= EstMB_top10, "Top 10%", "Bottom 90%"))


# STEP 2 - intial visualization
#-----------------------------#

# Save Kaplan-Meier plots for all endpoints colored by the binary MIS score (EstMB cutoff)
for (i in unique(survival_data$endpoint)){
  
  # Kaplan-Meier plot for each endpoint separately
  survival_data_run = survival_data %>% 
    dplyr::filter(endpoint == i)
  
  # Local dataset top decile vs others
  survplot_local_top10 = survfit2(Surv(time, status) ~ MIS_local_top10, data = survival_data_run) %>% 
    ggsurvfit(size = 1) +
    add_confidence_interval() + 
    scale_color_manual(values = c("slategray3", "gold2")) + 
    scale_fill_manual(values = c("slategray3", "gold2")) + 
    ggtitle(i) + 
    xlab("Days") + 
    ylab("Overall survival probability") + 
    theme(title = element_text(size = 18), 
          axis.title = element_text(size = 16), 
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  # EstMB threshold top decile vs others
  survplot_EstMB_top10 = survfit2(Surv(time, status) ~ MIS_EstMB_top10, data = survival_data_run) %>% 
    ggsurvfit(size = 1) +
    add_confidence_interval() + 
    scale_color_manual(values = c("slategray3", "gold2")) + 
    scale_fill_manual(values = c("slategray3", "gold2")) + 
    ggtitle(i) + 
    xlab("Days") + 
    ylab("Overall survival probability") + 
    theme(title = element_text(size = 18), 
          axis.title = element_text(size = 16), 
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  # Save plot to Results folder
  ggsave(plot = survplot_local_top10, 
         filename = paste("Results/KMplot_localDecile_", i, ".png", sep = ""), 
         width = 9, height = 6)
  
  ggsave(plot = survplot_EstMB_top10, 
         filename = paste("Results/KMplot_EstDecile_", i, ".png", sep = ""), 
         width = 9, height = 6)
}



# STEP 2 - run Cox models
#-----------------------------#

# Run CoxPH models for all endpoints separately
Cox_assumptions_df <- data.frame()
Cox_results_df <- data.frame()
for (i in unique(survival_data$endpoint)){
  
  # Run models - MIS as a continuous variable, MIS as a binary variable (local/EstMB cutoff)
  survival_data_run = survival_data %>% 
    dplyr::filter(endpoint == i) # Remove prevalent cases if still in the data
    
  
  # NB! Depending on the phenotypes available, include additional covariates
  cox_continuous = coxph(Surv(time, status) ~ BMI + age + gender + scale(MIS), data = survival_data)
  cox_localDecile = coxph(Surv(time, status) ~ BMI + age + gender + MIS_local_top10, data = survival_data)
  cox_EstDecile = coxph(Surv(time, status) ~ BMI + age + gender + MIS_EstMB_top10, data = survival_data)
  
  # Assumptions
  ph_continuous = cox.zph(cox_continuous)$table %>% as.data.frame() %>% 
    dplyr::mutate(type = "continuous", endpoint = i)
  
  ph_localDecile = cox.zph(cox_localDecile)$table %>% as.data.frame() %>% 
    dplyr::mutate(type = "localDecile", endpoint = i)
  
  ph_EstDecile = cox.zph(cox_EstDecile)$table %>% as.data.frame() %>% 
    dplyr::mutate(type = "EstMBDecile", endpoint = i)
  
  # Clean model output
  cox_output_continuous = broom::tidy(cox_continuous) %>% 
    dplyr::mutate(type = "continuous", endpoint = i)
  
  cox_output_localDecile = broom::tidy(cox_localDecile) %>% 
    dplyr::mutate(type = "localDecile", endpoint = i)
  
  cox_output_EstDecile = broom::tidy(cox_EstDecile) %>% 
    dplyr::mutate(type = "EstDecile", endpoint = i)
  
  # Add data to output df
  Cox_assumptions_df = dplyr::bind_rows(Cox_assumptions_df, ph_continuous, ph_localDecile, ph_EstDecile)
  Cox_results_df = dplyr::bind_rows(Cox_results_df, cox_output_continuous, cox_output_localDecile, cox_output_EstDecile)
} 

# Save 
saveRDS(Cox_results_df, "Results/CoxPH_results.rds")
saveRDS(Cox_assumptions_df, "Results/CoxPH_assumptions.rds")


