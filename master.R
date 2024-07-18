#' ============================= MASTER SCRIPT =================================
#' For the Blue-Cloud 2026 project workbench
#' To be translated to D4-Science and notebook at the end
#' A. Schickele 2023
#' =============================================================================
#' ============================== TO DO LIST ===================================
#' - check data access service as input
#' - check blue cloud data miner compatibility of the functions
#' =============================================================================

# --- 0. Start up and load functions
# All will be called in the config file later
rm(list=ls())
closeAllConnections()
setwd("/nfs/meso/work/aschickele/CEPHALOPOD")
source(file = "./code/00_config.R")
run_name <- "cocco_pfirst_shannon_omic_auto_pred_v2"
MAX_CLUSTER = 30

# --- 1. List the available species
# Within the user defined selection criteria
list_bio <- list_bio_wrapper(FOLDER_NAME = run_name,
                             DATA_SOURCE = "MAG",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, TARGET_MIN_DEPTH = 0, TARGET_MAX_DEPTH = 200, START_YEAR = 1950, STOP_YEAR = 2020))

# ------------------------------------------------------------------------------
# --- USER INPUT: Define the list of species to consider
sp_list <- list_bio %>%
  dplyr::filter(taxonrank == "MAG") %>% # all MAG
  dplyr::filter(Class == "Prymnesiophyceae") %>%# of class prymnesiophyceae
  dplyr::select(scientificname) %>% 
  unique() %>% pull() %>% .[!grepl("No match", .)]

# ------------------------------------------------------------------------------

# --- 2. Create the output folder, initialize parallelisation and parameters
# (1) Create an output folder containing all species-level runs, (2) Stores the 
# global parameters in an object, (3) Builds a local list of monthly raster
subfolder_list <- run_init(FOLDER_NAME = run_name,
                           SP_SELECT = sp_list,
                           WORMS_CHECK = FALSE,
                           FAST = TRUE,
                           LOAD_FROM = NULL,
                           DATA_TYPE = "proportions",
                           # ENV_VAR = NULL,
                           # ENV_VAR = c("!climatology_s_0_50","!climatology_s_200_300","!climatology_talk_SODA"),
                           ENV_VAR = c("!climatology_s_0_50","!climatology_s_200_300","!climatology_t_200_300","!climatology_A_200_300","!climatology_i_200_300","!climatology_n_200_300","!climatology_p_200_300","!climatology_o_200_300","!climatology_O_200_300"),
                           # ENV_VAR = c("climatology_M_0_0","climatology_i_0_50","climatology_t_0_50","climatology_p_0_50","climatology_n_0_50","climatology_omega_ca_SODA","climatology_A_PAR_regridded"),
                           ENV_PATH = "/nfs/meso/work/clercc/Predictors/PIPELINE_SET/VIRTUAL_SPECIES",
                           METHOD_PA = "density",
                           PER_RANDOM = 0,
                           PA_ENV_STRATA = TRUE,
                           OUTLIER = FALSE,
                           RFE = TRUE,
                           ENV_COR = 0.8,
                           NFOLD = 5,
                           FOLD_METHOD = "lon",
                           MODEL_LIST = c("GLM","MLP","BRT","GAM","SVM","RF"), # light version
                           LEVELS = 3,
                           TARGET_TRANSFORMATION = NULL,
                           # TARGET_TRANSFORMATION = "/nfs/meso/work/aschickele/CEPHALOPOD/function/target_transformation_yj_auto.R",
                           ENSEMBLE = TRUE,
                           N_BOOTSTRAP = 10,
                           CUT = 0)

# --- 3. Query biological data
# Get the biological data of the species we wish to model
mcmapply(FUN = query_bio_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 4. Query environmental data
# This functions returns an updated subfolder_list object to avoid computing
# species with less than the user defined minimum occurrence number
subfolder_list <- mcmapply(FUN = query_env,
                  FOLDER_NAME = run_name,
                  SUBFOLDER_NAME = subfolder_list,
                  mc.cores = min(length(subfolder_list), MAX_CLUSTERS), mc.preschedule = FALSE) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  .[grep("Error", ., invert = TRUE)] %>% # to exclude any API error or else
  as.vector()

# --- 5. Generate pseudo-absences if necessary
mcmapply(FUN = pseudo_abs,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 6. Outliers, Environmental predictor and MESS check 
# This functions returns an updated subfolder_list with meaningful feature set
subfolder_list <- mcmapply(FUN = query_check,
                           FOLDER_NAME = run_name,
                           SUBFOLDER_NAME = subfolder_list,
                           mc.cores = min(length(subfolder_list), MAX_CLUSTERS), mc.preschedule = FALSE) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  as.vector()

# --- 7. Generate split and re sampling folds
mcmapply(FUN = folds,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 8. Hyper parameters to train
hyperparameter(FOLDER_NAME = run_name)

# --- 9. Model fit
mcmapply(FUN = model_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 10. Model evaluation
# Performance metric and variable importance
mcmapply(FUN = eval_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# ---11. Model projections
mcmapply(FUN = proj_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 12. Output plots
# --- 12.1. Standard maps per algorithms
mcmapply(FUN = standard_maps,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 12.2. Partial dependency plots
# mcmapply(FUN = pdp,
#          FOLDER_NAME = run_name,
#          SUBFOLDER_NAME = subfolder_list,
#          mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 12.3 Diversity
# diversity_maps(FOLDER_NAME = run_name)

# --- 12.4 User synthesis
user_synthesis(FOLDER_NAME = run_name)

# --- END --- 