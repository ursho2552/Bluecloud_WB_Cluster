#' ============================= MASTER SCRIPT =================================
#' For the Blue-Cloud 2026 project workbench
#' To be translated to D4-Science and notebook at the end
#' A. Schickele 2023
#' =============================================================================
#' ============================== TO DO LIST ===================================
#' - data access service
#' - make sure the models are sending dataminer jobs
#' - have an easy and common architecture of files
#' - check the compatibility between dataminer and Biomod2...
#' - have a prototype wrapper for the 3 sub-pipelines
#' =============================================================================

# --- START UP
# All will be called in the config file later
rm(list=ls())
setwd("/net/meso/work/aschickele/Diversity")
source(file = "./code/00_config.R")

# --- 1. Query biological data
# First check which species are available -- TOO LONG : add key on DB
list_bio <- list_bio(DATA_TYPE = "pres",
                     SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016))

# Then query the species of interest -- TOO LONG : add key on DB
# Here WORMS 104464 is Calanus finmarchus, very known example
query <- query_bio(DATA_TYPE = "pres",
                   SP_SELECT = 104464,
                   SAMPLE_SELECT = list(MIN_SAMPLE = 50, MIN_DEPTH = 0, MAX_DEPTH = 50, START_YEAR = 1990, STOP_YEAR = 2016))

# --- 2. Query environmental data -- TOO LONG : find a solution for large data
query <- query_env(QUERY_BIO = query,
                   ENV_VAR = NULL,
                   ENV_PATH = "/net/meso/work/aschickele/Diversity/data/features_monthly")

# --- 3. Outliers and environmental covariance check
query <- query_check(QUERY = query,
                     OUTLIER = TRUE,
                     ENV_COR = 0.8,
                     MESS = TRUE)

# --- 4. Generate pseudo-absences if necessary -- BUG FIX : throw an error on some identical runs...
if(query$CALL$DATA_TYPE == "pres"){
  query <- pseudo_abs(QUERY = query,
                      METHOD_PA = "env")
}

# --- 5. Generate split and re sampling folds
query <- folds(QUERY = query,
               NFOLD = 5,
               FOLD_METHOD = "lon")

# --- 6. Hyper parameters to train
hp_list <- hyperparameter(QUERY = query,
                          MODEL_LIST = c("GLM","GAM","RF","MLP"),
                          LEVELS = 3)

# --- 7. Model fit -- FIX : RF is very long for big data
models <- model_wrapper(QUERY = query,
                        HP = hp_list,
                        MODEL_LIST = hp_list$CALL$MODEL_LIST)
# Removing hp_list as it is transferred to the models object
rm(hp_list)

# --- 8. Model evaluation
models <- eval_wrapper(QUERY = query,
                       MODELS = models)

# --- 9. Model projections
if(length(models$CALL$MODEL_LIST) >= 1){
  models <- proj_wrapper(QUERY = query,
                         MODELS = models,
                         N_BOOTSTRAP = 10,
                         PROJ_PATH = NULL)
}

# --- 10. Output plots
standard_maps(QUERY = query,
              MODELS = models,
              ENSEMBLE = TRUE,
              MESS = FALSE)

var_imp(QUERY = query,
        MODELS = models,
        ENSEMBLE = TRUE)

pdp(QUERY = query,
    MODELS = models,
    N_BOOTSTRAP = 10,
    ENSEMBLE = TRUE)














# --- END --- 