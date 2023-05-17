#' ============================= MASTER SCRIPT =================================
#' For the Blue-Cloud 2026 project workbench
#' To be translated to D4-Science and notebook at the end
#' A. Schickele 2023
#' =============================================================================
#' ============================== TO DO LIST ===================================
#' - check data access service as input
#' - check blue cloud data miner compatibility of the functions
#' - include MBTR in the prototypes for proportion data
#' =============================================================================

# --- 0. Start up and load functions
# All will be called in the config file later
rm(list=ls())
setwd("/net/meso/work/aschickele/Bluecloud_WB_local")
source(file = "./code/00_config.R")
run_name <- "new_data_access2"

# --- 1a. List the available species
# Within the user defined selection criteria
list_bio <- list_bio_wrapper(FOLDER_NAME = run_name,
                             DATA_SOURCE = "omic",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016))

# Define the list of species to consider
sp_list <- c("5820", "9760")

# --- 1b. Create the output folder and initialize parallelisation
# Create an output folder containing all species-level runs
folder_init(FOLDER_NAME = run_name,
            SP_SELECT = sp_list,
            LOAD_FROM = "new_data_access",
            DATA_TYPE = "pres")

# Define the list of sub folders to parallelize on
subfolder_list <- list.dirs(paste0(project_wd, "/output/", run_name), full.names = FALSE, recursive = FALSE)

# --- 2a. Query biological data
# Get the biological data of the species we wish to model
mcmapply(FUN = query_bio_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 2b. Query environmental data 
mcmapply(FUN = query_env,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         ENV_PATH = "/net/meso/work/aschickele/Bluecloud_WB_local/data/features_monthly",
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 3. Outliers and environmental covariance check
mcmapply(FUN = query_check,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         OUTLIER = TRUE,
         ENV_COR = 0.8,
         MESS = TRUE,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 4. Generate pseudo-absences if necessary -- BUG FIX : throw an error on some identical runs...
mcmapply(FUN = pseudo_abs,
         FOLDER_NAME = run_name,
           SUBFOLDER_NAME = subfolder_list,
           METHOD_PA = "env",
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 5. Generate split and re sampling folds
mcmapply(FUN = folds,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         NFOLD = 5,
         FOLD_METHOD = "lon",
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 6. Hyper parameters to train
hyperparameter(FOLDER_NAME = run_name,
               MODEL_LIST = c("GLM","GAM","RF","MLP"),
               LEVELS = 3,
               mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 7. Model fit -- FIX : RF is very long for big data
mcmapply(FUN = model_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 8. Model evaluation
mcmapply(FUN = eval_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 9. Model projections
mcmapply(FUN = proj_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         N_BOOTSTRAP = 10,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 10. Output plots
# --- 10.1. Standard maps per algorithms
mcmapply(FUN = standard_maps,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         ENSEMBLE = TRUE,
         MESS = FALSE,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 10.2. Variable importance output
mcmapply(FUN = var_imp,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         ENSEMBLE = TRUE,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS))

# --- 10.3. Partial dependency plots
pdp(SP_SELECT = 104464,
    FOLDER_NAME = run_name,
    N_BOOTSTRAP = 10,
    ENSEMBLE = TRUE)














# --- END --- 