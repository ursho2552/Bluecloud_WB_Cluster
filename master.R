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
setwd("/net/meso/work/aschickele/Diversity")
source(file = "./code/00_config.R")
run_name <- "test_run"

# --- 1a. List the available species
# Within the user defined selection criteria
list_bio <- list_bio(DATA_TYPE = "pres",
                     SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016))

# --- 1b. Create the output folder
# Create an output folder containing all species-level runs
folder_init(DATA_TYPE = "pres",
                       SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016),
                       SP_SELECT = c(104464, 101, 102),
                       FOLDER_NAME = run_name)

# --- 2a. Query biological data
# Then query the species of interest according to their WORMS_ID
query_bio(SP_SELECT = 104464,
          FOLDER_NAME = run_name)

# Test if working in multivariate -- YES -- TO IMPLEMENT TO ALL FUNCTIONS !!
mcmapply(FUN = query_bio,
         SP_SELECT = c(104464,102),
         FOLDER_NAME = run_name)

# --- 2b. Query environmental data -- TOO LONG : find a solution for large data
query_env(SP_SELECT = 104464,
          FOLDER_NAME = run_name,
          ENV_VAR = NULL,
          ENV_PATH = "/net/meso/work/aschickele/Diversity/data/features_monthly")

# --- 3. Outliers and environmental covariance check
query_check(SP_SELECT = 104464,
            FOLDER_NAME = run_name,
            OUTLIER = TRUE,
            ENV_COR = 0.8,
            MESS = TRUE)

# --- 4. Generate pseudo-absences if necessary -- BUG FIX : throw an error on some identical runs...
pseudo_abs(SP_SELECT = 104464,
           FOLDER_NAME = run_name,
           METHOD_PA = "env")

# --- 5. Generate split and re sampling folds
folds(SP_SELECT = 104464,
      FOLDER_NAME = run_name,
      NFOLD = 5,
      FOLD_METHOD = "lon")

# --- 6. Hyper parameters to train - CHECK WHERE TO MOVE
# Maybe this should be moved with the model list, environmental variable correlation etc...
# To the beginning of the script ?
hyperparameter(FOLDER_NAME = run_name,
               MODEL_LIST = c("GLM","GAM","RF","MLP"),
               LEVELS = 3)

# --- 7. Model fit -- FIX : RF is very long for big data
model_wrapper(SP_SELECT = 104464,
              FOLDER_NAME = run_name,
              MODEL_LIST = NULL)

# --- 8. Model evaluation
eval_wrapper(SP_SELECT = 104464,
             FOLDER_NAME = run_name)

# --- 9. Model projections
proj_wrapper(SP_SELECT = 104464,
             FOLDER_NAME = run_name,
             N_BOOTSTRAP = 10,
             PROJ_PATH = NULL)

# --- 10. Output plots
# --- 10.1. Standard maps per algorithms
standard_maps(SP_SELECT = 104464,
              FOLDER_NAME = run_name,
              ENSEMBLE = TRUE,
              MESS = FALSE)

# --- 10.2. Variable importance output
var_imp(SP_SELECT = 104464,
        FOLDER_NAME = run_name,
        ENSEMBLE = TRUE)

# --- 10.3. Partial dependency plots
pdp(SP_SELECT = 104464,
    FOLDER_NAME = run_name,
    N_BOOTSTRAP = 10,
    ENSEMBLE = TRUE)














# --- END --- 