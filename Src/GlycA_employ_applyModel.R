

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
library("tidymodels") # or recipes package



# Read MetaPhLAn3 profiles 
#-----------------------------#
# samples as rownames, full taxonomies as colnames:
# eg k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Citrobacter|s__Citrobacter_farmeri)

# Example: load a dataset from curatedMetagenomicData package 
metaphlan_object <- curatedMetagenomicData::sampleMetadata %>% 
  dplyr::filter(study_name == "KarlssonFH_2013") %>% 
  curatedMetagenomicData::returnSamples("relative_abundance", counts = FALSE)

metaphlan_data <- assays(metaphlan_object)[[1]] %>% 
  t() %>% 
  as.data.frame()



# Read model information - regression coefficient, data preprocessing
#-----------------------------#

# Coefficients from the LASSO model
in_coefficients_df <- readRDS("RData/MBscore_EstMB_dataPreprocessingPredictors_fixedCovariates_CLR_LASSO.rds")

# Data preprocessing steps - recipe 
in_data_preprocessing <- readRDS("RData/MBscore_EstMB_dataPreprocessingRecipe_fixedCovariates_CLR_LASSO.rds")








# Preprocess microbiome data ----
#------------------------------------------------#
#                                                #
#          PREPROCESS MICROBIOME DATA            # 
#                                                #
#------------------------------------------------#


# STEP 1 - renaming columns
#-----------------------------#

# Mapping of colnames - use only species name without "s__"
new_colnames_df <- data.frame(taxonomy = colnames(metaphlan_data)) %>% 
  dplyr::mutate(species = substring(stringr::str_replace(taxonomy, "^.+\\|", ""), 4))

# Rename
metaphlan_data_renamed <- metaphlan_data %>% 
  dplyr::rename_with(~new_colnames_df$species, new_colnames_df$taxonomy)




# STEP 2 - calculate diversity metrics
#-----------------------------#
diversity_df <- data.frame(sample_id = rownames(metaphlan_data), 
                           observed = rowSums(sign(metaphlan_data)), 
                           shannon = vegan::diversity(metaphlan_data))



# STEP 3 - make sure all columns used for model building are present. Add zero columns where needed
#-----------------------------#

# List of species used in the EstMB model development
EstMB_species <- in_data_preprocessing$var_info %>% 
  dplyr::filter(!(variable %in% c("BMI", "gender", "Age_at_MBsample", "target_var", "observed", "shannon"))) %>% 
  dplyr::pull(variable)

# Identify species only present in EstMB data and add them as zero columns
EstMB_specific_species <- EstMB_species[!(EstMB_species %in% colnames(metaphlan_data_renamed))]

metaphlan_data_allCols <- metaphlan_data_renamed
for (i in EstMB_specific_species){
  metaphlan_data_allCols[ ,i] = 0
}


# STEP 4 - CLR transformation 
#-----------------------------#

# Define pseudocount
pseudocount <- min(metaphlan_data_renamed[metaphlan_data_renamed != 0])/2

# Replace zeros with pseudocount
metaphlan_data_zeroImputed <- metaphlan_data_allCols
metaphlan_data_zeroImputed[metaphlan_data_zeroImputed == 0] <- pseudocount

# CLR-transformation
metaphlan_data_CLR <- compositions::clr(metaphlan_data_zeroImputed[ ,EstMB_species]) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample_id")

# Merge diversity data
microbiome_data <- metaphlan_data_CLR %>%
  dplyr::left_join(diversity_df, by = "sample_id")



# STEP 5 - standardization 
#-----------------------------#

# Define scaling coefficients
scaling_step_bake <- in_data_preprocessing$steps[[3]]
scaling_df <- data.frame(taxa = names(scaling_step_bake$means), 
                         mean = scaling_step_bake$means, 
                         sd = scaling_step_bake$sds)

# Apply scaling for the new data
microbiome_data_scaled <- microbiome_data %>% dplyr::select(sample_id)
for (i in setdiff(colnames(microbiome_data), "sample_id")){
  
  # Coefficients for scaling
  run_taxa <- microbiome_data %>% dplyr::pull(i)
  run_mean = scaling_df %>% dplyr::filter(taxa == i) %>% dplyr::pull(mean)
  run_sd = scaling_df %>% dplyr::filter(taxa == i) %>% dplyr::pull(sd)
  
  # Scaled values for the taxa
  taxa_scaled <- (run_taxa - run_mean)/run_sd
  
  # Add to the new data
  microbiome_data_scaled[i] = taxa_scaled
}





# Calculate the microbiome inflammation score ----
#------------------------------------------------#
#                                                #
#    CALCULATE MICROBIOME INFLAMMATION SCORE     # 
#                                                #
#------------------------------------------------#

# Calculate score
MIS_result_df <- microbiome_data_scaled %>% 
  tidyr::gather(taxa, scaled_value, -sample_id) %>% 
  dplyr::left_join(in_coefficients_df, by = c("taxa" = "term")) %>% 
  dplyr::group_by(sample_id) %>% 
  dplyr::summarise(MIS = sum(scaled_value*estimate, na.rm = T)) %>% 
  dplyr::ungroup()


# Save results
saveRDS(MIS_result_df, file = "Results/Microbiome_inflammation_score.rds")
