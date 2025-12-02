

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
library("purrr")



# Read calculated microbiome score
#-----------------------------#
MIS_result_df <- readRDS(file = "Results/Microbiome_inflammation_score_MetaPhlan4.rds") %>% 
  dplyr::rename("MIS_MP4" = "MIS") %>% 
  dplyr::left_join(readRDS(file = "Results/Microbiome_inflammation_score_MetaPhlan3.rds"), by = "sample_id") %>% 
  dplyr::rename("MIS_MP3" = "MIS")

# Local top 20%
MP3_local_top20 <- quantile(MIS_result_df$MIS_MP3, probs = 0.8)
MP4_local_top20 <- quantile(MIS_result_df$MIS_MP4, probs = 0.8)

# EstMB top 20% 
MP3_EstMB_top20 <- readRDS("Results/MIS_EstMB_top20percentile_PA_LASSOF.rds")
MP4_EstMB_top20 <- readRDS("Results/MIS_EstMB_top20percentile_PA4_LASSOF.rds")




# Read phenotype data
#-----------------------------#

#NB: REPLACE!
# Example: generate a dataset from the curatedMetagenomicData package - imitate survival data
# TODO: include multiple phenotypes in a long data format; include relevant and available covariates; exclude prevalent cased; apply other exlusion criteria where needed
set.seed(0)
phenotype_data <- curatedMetagenomicData::sampleMetadata %>% 
  dplyr::filter(study_name == "KarlssonFH_2013") %>% 
  tibble::rowid_to_column(var = "id") %>% 
  dplyr::filter(id %in% 1:100) %>% 
  dplyr::mutate(sample_id = MIS_result_df$sample_id) %>% 
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
# 1 MV_FEI1_t1Q14 24.9 57.37363      1      0  495      T2D
# 2 MV_FEI2_t1Q14 32.7 44.37843      1      0  617      T2D
# 3 MV_FEI3_t1Q14 29.7 49.80367      0      0  651      T2D
# 4 MV_FEI4_t1Q14 18.6 53.17662      1      0 1168      T2D
# 5 MV_FEI4_t2Q15 26.6 44.34391      1      0   89      T2D
# 6 MV_FEI5_t1Q14 26.5 42.84029      0      0  923      T2D





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
  dplyr::mutate(MP3_MIS_local_top20 = ifelse(MIS_MP3 >= MP3_local_top20, "Top 20%", "Bottom 80%"),
                MP4_MIS_local_top20 = ifelse(MIS_MP4 >= MP4_local_top20, "Top 20%", "Bottom 80%"),
                
                MP3_MIS_EstMB_top20 = ifelse(MIS_MP3 >= MP3_EstMB_top20, "Top 20%", "Bottom 80%"),
                MP4_MIS_EstMB_top20 = ifelse(MIS_MP4 >= MP4_EstMB_top20, "Top 20%", "Bottom 80%"))


# STEP 2 - intial visualization
#-----------------------------#

# Save Kaplan-Meier plots for all endpoints colored by the binary MIS score (EstMB cutoff)
for (i in unique(survival_data$endpoint)){
  
  # Kaplan-Meier plot for each endpoint separately
  survival_data_run = survival_data %>% 
    dplyr::filter(endpoint == i)
  
  # Local dataset top decile vs others
  #-----------------------------#
  # Metaphlan 3 
  MP3_survplot_local_top20 = survfit2(Surv(time, status) ~ MP3_MIS_local_top20, data = survival_data_run) %>% 
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
  
  # Metaphlan 4
  MP4_survplot_local_top20 = survfit2(Surv(time, status) ~ MP4_MIS_local_top20, data = survival_data_run) %>% 
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
  #-----------------------------#
  # Metaphlan 3 
  MP3_survplot_EstMB_top20 = survfit2(Surv(time, status) ~ MP3_MIS_EstMB_top20, data = survival_data_run) %>% 
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
  
  # Metaphlan 4
  MP4_survplot_EstMB_top20 = survfit2(Surv(time, status) ~ MP4_MIS_EstMB_top20, data = survival_data_run) %>% 
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
  #-----------------------------#
  ggsave(plot = MP3_survplot_local_top20, filename = paste("Results/KMplot_Metaphlan3_localDecile_", i, ".png", sep = ""), width = 9, height = 6)
  ggsave(plot = MP3_survplot_local_top20, filename = paste("Results/KMplot_Metaphlan3_localDecile_", i, ".pdf", sep = ""), width = 9, height = 6)
  
  ggsave(plot = MP4_survplot_local_top20, filename = paste("Results/KMplot_Metaphlan4_localDecile_", i, ".png", sep = ""), width = 9, height = 6)
  ggsave(plot = MP4_survplot_local_top20, filename = paste("Results/KMplot_Metaphlan4_localDecile_", i, ".pdf", sep = ""), width = 9, height = 6)
  
  ggsave(plot = MP3_survplot_EstMB_top20, filename = paste("Results/KMplot_Metaphlan3_EstDecile_", i, ".png", sep = ""), width = 9, height = 6)
  ggsave(plot = MP3_survplot_EstMB_top20, filename = paste("Results/KMplot_Metaphlan3_EstDecile_", i, ".pdf", sep = ""), width = 9, height = 6)
  
  ggsave(plot = MP4_survplot_EstMB_top20, filename = paste("Results/KMplot_Metaphlan4_EstDecile_", i, ".png", sep = ""), width = 9, height = 6)
  ggsave(plot = MP4_survplot_EstMB_top20, filename = paste("Results/KMplot_Metaphlan4_EstDecile_", i, ".pdf", sep = ""), width = 9, height = 6)
  
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
  # Run models, clean output, test assumptions
  MP3_cox_continuous = coxph(Surv(time, status) ~ BMI + age + gender + scale(MIS_MP3), data = survival_data_run)
  MP3_cox_output_continuous = broom::tidy(MP3_cox_continuous) %>% dplyr::mutate(type = "continuous", bioinf = "Metaphlan3", endpoint = i) # Clean model output
  MP3_ph_continuous = cox.zph(MP3_cox_continuous)$table %>% as.data.frame() %>% dplyr::mutate(type = "continuous", bioinf = "Metaphlan3", endpoint = i) # Assumptions
  
  MP3_cox_localDecile = coxph(Surv(time, status) ~ BMI + age + gender + MP3_MIS_local_top20, data = survival_data_run)
  MP3_cox_output_localDecile = broom::tidy(MP3_cox_localDecile) %>% dplyr::mutate(type = "localDecile", bioinf = "Metaphlan3", endpoint = i) # Clean model output
  MP3_ph_localDecile = cox.zph(MP3_cox_localDecile)$table %>% as.data.frame() %>% dplyr::mutate(type = "localDecile", bioinf = "Metaphlan3", endpoint = i) # Assumptions
  
  MP4_cox_continuous = coxph(Surv(time, status) ~ BMI + age + gender + scale(MIS_MP4), data = survival_data_run)
  MP4_cox_output_continuous = broom::tidy(MP4_cox_continuous) %>% dplyr::mutate(type = "continuous", bioinf = "Metaphlan4", endpoint = i) # Clean model output
  MP4_ph_continuous = cox.zph(MP4_cox_continuous)$table %>% as.data.frame() %>% dplyr::mutate(type = "continuous", bioinf = "Metaphlan4", endpoint = i) # Assumptions
  
  MP4_cox_localDecile = coxph(Surv(time, status) ~ BMI + age + gender + MP4_MIS_local_top20, data = survival_data_run)
  MP4_cox_output_localDecile = broom::tidy(MP4_cox_localDecile) %>% dplyr::mutate(type = "localDecile", bioinf = "Metaphlan4", endpoint = i) # Clean model output
  MP4_ph_localDecile = cox.zph(MP4_cox_localDecile)$table %>% as.data.frame() %>% dplyr::mutate(type = "localDecile", bioinf = "Metaphlan4", endpoint = i) # Assumptions
  
  
  # Add data to output df
  Cox_assumptions_df = dplyr::bind_rows(Cox_assumptions_df, MP3_ph_continuous, MP4_ph_continuous, 
                                        MP3_ph_localDecile, MP4_ph_localDecile)
  
  Cox_results_df = dplyr::bind_rows(Cox_results_df, MP3_cox_output_continuous, MP4_cox_output_continuous, 
                                    MP3_cox_output_localDecile, MP4_cox_output_localDecile)
  
  
  # Run models for EstMB decile when possible
  check_MP3 = survival_data_run %>% dplyr::group_by(MP3_MIS_EstMB_top20, status) %>% dplyr::summarise(n = n())
  if (nrow(check_MP3) == 4 & min(check_MP3$n) >= 5){
    MP3_cox_EstDecile = coxph(Surv(time, status) ~ BMI + age + gender + MP3_MIS_EstMB_top20, data = survival_data_run)
    MP3_cox_output_EstDecile = broom::tidy(MP3_cox_EstDecile) %>% dplyr::mutate(type = "EstDecile", bioinf = "Metaphlan3", endpoint = i) # Clean model output
    MP3_ph_EstDecile = cox.zph(MP3_cox_EstDecile)$table %>% as.data.frame() %>% dplyr::mutate(type = "EstMBDecile", bioinf = "Metaphlan3", endpoint = i) # Assumptions
    
    Cox_results_df = dplyr::bind_rows(Cox_results_df, MP3_cox_output_EstDecile)
    Cox_assumptions_df = dplyr::bind_rows(Cox_assumptions_df, MP3_ph_EstDecile)
  }

  check_MP4 = survival_data_run %>% dplyr::group_by(MP4_MIS_EstMB_top20, status) %>% dplyr::summarise(n = n())
  if (nrow(check_MP4) == 4 & min(check_MP4$n) >= 5){
    MP4_cox_EstDecile = coxph(Surv(time, status) ~ BMI + age + gender + MP4_MIS_EstMB_top20, data = survival_data_run)
    MP4_cox_output_EstDecile = broom::tidy(MP4_cox_EstDecile) %>% dplyr::mutate(type = "EstDecile", bioinf = "Metaphlan4", endpoint = i) # Clean model output
    MP4_ph_EstDecile = cox.zph(MP4_cox_EstDecile)$table %>% as.data.frame() %>% dplyr::mutate(type = "EstMBDecile", bioinf = "Metaphlan4", endpoint = i) # Assumptions
    
    Cox_results_df = dplyr::bind_rows(Cox_results_df, MP3_cox_output_EstDecile)
    Cox_assumptions_df = dplyr::bind_rows(Cox_assumptions_df, MP3_ph_EstDecile)
  }
} 

# Save 
saveRDS(Cox_results_df, "Results/CoxPH_results.rds")
saveRDS(Cox_assumptions_df, "Results/CoxPH_assumptions.rds")


