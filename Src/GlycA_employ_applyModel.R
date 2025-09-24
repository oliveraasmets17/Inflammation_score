

# Program setup ----
#------------------------------------------------#
#                                                #
#                 PROGRAM SETUP                  # 
#                                                #
#------------------------------------------------#

# Load packages
#-----------------------------#
library("dplyr")
library("stringr")
library("vegan")
library("curatedMetagenomicData")
library("tidymodels") # for recipes package



# Read model information - regression coefficient, data preprocessing
#-----------------------------#

# TODO - define profiling & paths to microbiome data
profiling_used <- "MetaPhlan4"

if (profiling_used == "MetaPhlan4"){
  
  # Coefficients from the LASSO model
  in_coefficients_df <- readRDS("Results/applyMIS_predictors_EstMB_PA4_LASSOF.rds")
  
  # Data preprocessing steps - recipe 
  in_data_preprocessing <- readRDS("Results/applyMIS_recipe_EstMB_PA4_LASSOF.rds")
  
} else if (profiling_used == "MetaPhlan3"){
  
  # Coefficients from the LASSO model
  in_coefficients_df <- readRDS("Results/applyMIS_predictors_EstMB_PA_LASSOF.rds")
  
  # Data preprocessing steps - recipe 
  in_data_preprocessing <- readRDS("Results/applyMIS_recipe_EstMB_PA_LASSOF.rds")
}



# Read MetaPhLAn profiles 
#-----------------------------#
if (profiling_used == "MetaPhlan3"){
  #TODO: read Metaphlan 3 profile here
  # samples as rownames, species names as colnames:
  # eg s__Citrobacter_farmeri
  
  metaphlan_countdata <- readRDS(file = "Examples/metaPhlan3_exampledata.rds")
  
} else if(profiling_used == "MetaPhlan4"){
  #TODO: read Metaphlan 4 profile here
  # samples as rownames, species names as colnames:
  # eg s__Citrobacter_farmeri
  
  metaphlan_countdata <- readRDS(file = "Examples/metaPhlan4_exampledata.rds")
}









# Prepare data, calculate MIS ----
#------------------------------------------------#
#                                                #
#    CALCULATE MICROBIOME INFLAMMATION SCORE     # 
#                                                #
#------------------------------------------------#

# Make sure all columns used for model building are present. Add zero columns where needed
#-----------------------------#
# List of species used in the EstMB model development
EstMB_species <- in_data_preprocessing$var_info %>% 
  dplyr::filter(!(variable %in% c("BMI", "gender", "Age_at_MBsample", "target_var"))) %>% 
  dplyr::pull(variable)

# Identify species only present in EstMB data and add them as zero columns
EstMB_specific_species <- EstMB_species[!(EstMB_species %in% colnames(metaphlan_countdata))]

metaphlan_data_allCols <- metaphlan_countdata
for (i in EstMB_specific_species){
  metaphlan_data_allCols[ ,i] = 0
}



# Presence-absence transformation 
#-----------------------------#
metaphlan_countdata_pa <- sign(metaphlan_data_allCols) %>% 
  as.data.frame()



# Calculate the microbiome inflammation score
#-----------------------------#
MIS_result_df <- metaphlan_countdata_pa %>% 
  tibble::rownames_to_column(var = "sample_id") %>% 
  tidyr::gather(term, value, -sample_id) %>% 
  dplyr::left_join(in_coefficients_df, by = c("term")) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(MIS = sum(value*estimate, na.rm = T)) %>% 
  dplyr::ungroup()


# Save results
saveRDS(MIS_result_df, file = paste("Results/Microbiome_inflammation_score_", profiling_used, ".rds", sep = ""))
