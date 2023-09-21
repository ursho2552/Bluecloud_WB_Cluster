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
setwd("/net/meso/work/aschickele/Bluecloud_WB_local")
source(file = "./code/00_config.R")
run_name <- "Ensemble_eval4"

# --- 1. List the available species
# Within the user defined selection criteria
list_bio <- list_bio_wrapper(FOLDER_NAME = run_name,
                             DATA_SOURCE = "abundance",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, TARGET_MIN_DEPTH = 0, TARGET_MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2020))

# Define the list of species to consider
sp_list <- list_bio$worms_id %>% unique()
sp_list <- list_bio %>%
  dplyr::filter(grepl("Tripos ", scientificname)) %>%
  dplyr::select(worms_id) %>%
  unique() %>% pull()

# --- 2. Create the output folder, initialize parallelisation and parameters
# (1) Create an output folder containing all species-level runs, (2) Stores the 
# global parameters in an object, (3) Builds a local list of monthly raster
subfolder_list <- run_init(FOLDER_NAME = run_name,
                           SP_SELECT = sp_list,
                           FAST = FALSE,
                           LOAD_FROM = NULL,
                           DATA_TYPE = "continuous",
                           ENV_VAR = c("woa18_all_i_allmonths_surface",
                                       "woa18_all_o_allmonths_surface",
                                       "woa18_all_n_allmonths_surface",
                                       "woa18_all_p_allmonths_surface",
                                       "woa18_decav_t_allmonths_surfaceSST",
                                       "chlor_a_SeaWIFS_allmonths"),
                           ENV_PATH = "/net/meso/work/nknecht/Masterarbeit/General_Pipeline/Data/environmental_climatologies",
                           METHOD_PA = "density",
                           PER_RANDOM = 0.05,
                           OUTLIER = TRUE,
                           UNIVARIATE = FALSE,
                           ENV_COR = 0.8,
                           NFOLD = 3,
                           FOLD_METHOD = "lon",
                           MODEL_LIST = c("GLM","GAM","RF","MLP"),
                           LEVELS = 3,
                           ENSEMBLE = TRUE,
                           N_BOOTSTRAP = 10,
                           CUT = 0)

# --- 3. Query biological data
# Get the biological data of the species we wish to model
mcmapply(FUN = query_bio_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 4. Query environmental data
# This functions returns an updated subfolder_list object to avoid computing
# species with less than the user defined minimum occurrence number
subfolder_list <- mcmapply(FUN = query_env,
                  FOLDER_NAME = run_name,
                  SUBFOLDER_NAME = subfolder_list,
                  mc.cores = min(length(subfolder_list), MAX_CLUSTERS)) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  as.vector()

# --- 5. Generate pseudo-absences if necessary
mcmapply(FUN = pseudo_abs,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 6. Outliers, Environmental predictor and MESS check 
mcmapply(FUN = query_check,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 7. Generate split and re sampling folds
mcmapply(FUN = folds,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 8. Hyper parameters to train
hyperparameter(FOLDER_NAME = run_name)

# --- 9. Model fit -- FIX : RF is very long for big data
mcmapply(FUN = model_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 10. Model evaluation
# Performance metric and variable importance
mcmapply(FUN = eval_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# ---11. Model projections
mcmapply(FUN = proj_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 12. Output plots
# --- 12.1. Standard maps per algorithms
mcmapply(FUN = standard_maps,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 12.2. Partial dependency plots - TAKES AGES FOR LARGE OCCURRENCE NUMBER
mcmapply(FUN = pdp,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 12.3 Diversity
diversity_maps(FOLDER_NAME = run_name,
               SUBFOLDER_NAME = subfolder_list,
               BUFFER = 1,
               N_BOOTSTRAP = 10)

# --- 12.4 User synthesis
user_synthesis(FOLDER_NAME = run_name)

# --- END --- 