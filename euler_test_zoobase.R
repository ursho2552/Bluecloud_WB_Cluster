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
setwd("/nfs/kryo/work/ursho/PhD/Projects/BlueCloud/Bluecloud_WB_Cluster")
suppressPackageStartupMessages(source(file = "./code/my_config.R"))
system_time <- format(Sys.time(), "%Y_%m_%d_%H_%M")

run_name <- paste("ZooBASE_test", system_time, sep = "_")
#parse arguments from command line and overwrite MAX_CLUSTERS variable
parseCommandArgs()
MAX_CLUSTERS <- as.numeric(cores)

# library(profvis)
# p <- profvis({
# --- 1. List the available species
# Within the user defined selection criteria
message(paste(Sys.time(), ": list bio start"))
list_bio <- list_bio_wrapper(FOLDER_NAME = run_name,
                             DATA_SOURCE = "/nfs/sea/work/public/shared/AtlantECO/BASE/AtlantECO-BASE-v1_microbiome_traditional_zooplankton_species_occurrences_ZooBasev2_20220909.csv",
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, TARGET_MIN_DEPTH = 0, TARGET_MAX_DEPTH = 200, START_YEAR = 1950, STOP_YEAR = 2020))
message(paste(Sys.time(), ": list bio done"))


# ------------------------------------------------------------------------------
# --- USER INPUT: Define the list of species to consider
sp_list <- list_bio %>%
  dplyr::filter(taxonrank == "Species") %>%
  dplyr::select(worms_id) %>%
  unique() %>% pull() %>% .[!grepl("No match", .)]


# --- 2. Create the output folder, initialize parallelisation and parameters
# (1) Create an output folder containing all species-level runs, (2) Stores the
# global parameters in an object, (3) Builds a local list of monthly raster
# =====================================================
# CALL OBJECT IS CHANGED HERE, BUT NOT RUN IN PARALLEL
# =====================================================
message(paste(Sys.time(), ": run start"))
subfolder_list <- run_init(FOLDER_NAME = run_name,
                           SP_SELECT = sp_list,
                           WORMS_CHECK = FALSE,
                           FAST = TRUE,
                           LOAD_FROM = NULL,
                           DATA_TYPE = "presence_only",
                           ENV_VAR = c("!climatology_s_0_50","!climatology_s_200_300",
                                        "!climatology_t_200_300","!climatology_A_200_300",
                                        "!climatology_i_200_300","!climatology_n_200_300",
                                        "!climatology_p_200_300","!climatology_o_200_300",
                                        "!climatology_O_200_300"),
                           ENV_PATH = "/nfs/meso/work/clercc/Predictors/PIPELINE_SET/VIRTUAL_SPECIES",
                           METHOD_PA = "density",
                           PER_RANDOM = 0,
                           PA_ENV_STRATA = TRUE,
                           OUTLIER = TRUE,
                           RFE = TRUE,
                           ENV_COR = 0.8,
                           NFOLD = 5,
                           FOLD_METHOD = "lon",
                           MODEL_LIST = c("GLM","MLP","BRT","GAM","SVM","RF"), # light version
                           LEVELS = 3,
                           TARGET_TRANSFORMATION = NULL,
                           ENSEMBLE = TRUE,
                           N_BOOTSTRAP = 10,
                           CUT = 0.1)
message(paste(Sys.time(), ": run done"))


# --- 3. Query biological data
# Get the biological data of the species we wish to model
message(paste0('Load CALL object'))
load(paste0(project_wd, "/output/", run_name, "/CALL.RData"))

message(paste(Sys.time(), ": query bio start"))
mcmapply(FUN = query_bio_wrapper,
         SUBFOLDER_NAME = subfolder_list,
         MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = PRESCHEDULE)
message(paste(Sys.time(), ": query bio done"))

# --- 4. Query environmental data
# This functions returns an updated subfolder_list object to avoid computing
# species with less than the user defined minimum occurrence number
# ===============================================
# CALL OBJECT IS CHANGED HERE!! UHE 12/09/2024
# ===============================================

message(paste(Sys.time(), ": query env start"))
subfolder_list <- mcmapply(FUN = query_env,
                           SUBFOLDER_NAME = subfolder_list,
                           MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
                           mc.cores = min(length(subfolder_list), MAX_CLUSTERS), mc.preschedule = PRESCHEDULE) %>%
  unlist() %>%
  na.omit(subfolder_list) %>%
  .[grep("Error", ., invert = TRUE)] %>% # to exclude any API error or else
  as.vector()
message(paste(Sys.time(), ": query env done"))

message(paste0('Reload CALL object'))
load(paste0(project_wd, "/output/", run_name, "/CALL.RData"))

# --- 5. Generate pseudo-absences if necessary
message(paste(Sys.time(), ": PA start"))
mcmapply(FUN = pseudo_abs,
         SUBFOLDER_NAME = subfolder_list,
         MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = PRESCHEDULE)
message(paste(Sys.time(), ": PA done"))

# --- 6. Outliers, Environmental predictor and MESS check
# This functions returns an updated subfolder_list with meaningful feature set
message(paste(Sys.time(), ": query check start"))
subfolder_list <- mcmapply(FUN = query_check,
                           SUBFOLDER_NAME = subfolder_list,
                           MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
                           mc.cores = min(length(subfolder_list), MAX_CLUSTERS), mc.preschedule = PRESCHEDULE) %>%
  unlist() %>%
  na.omit(subfolder_list) %>%
  as.vector()
message(paste(Sys.time(), ": query check done"))


# --- 7. Generate split and re sampling folds
message(paste(Sys.time(), ": fold start"))
mcmapply(FUN = folds,
         SUBFOLDER_NAME = subfolder_list,
         MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = PRESCHEDULE)
message(paste(Sys.time(), ": fold done"))


# --- 8. Hyper parameters to train
# ===============================================
# CALL OBJECT IS CHANGED HERE!! UHE 12/09/2024
# ===============================================
hyperparameter(FOLDER_NAME = run_name)

message(paste0('Reload CALL object after hyperparameter'))
load(paste0(project_wd, "/output/", run_name, "/CALL.RData"))

# --- 9. Model fit
message(paste(Sys.time(), ": fit start"))
mcmapply(FUN = model_wrapper,
         SUBFOLDER_NAME = subfolder_list,
         MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = PRESCHEDULE)
message(paste(Sys.time(), ": fit done"))


# --- 10. Model evaluation
# Performance metric and variable importance
message(paste(Sys.time(), ": eval start"))
mcmapply(FUN = eval_wrapper,
         SUBFOLDER_NAME = subfolder_list,
         MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = PRESCHEDULE)
message(paste(Sys.time(), ": eval done"))



# ---11. Model projections
message(paste(Sys.time(), ": proj start"))
mcmapply(FUN = proj_wrapper,
         SUBFOLDER_NAME = subfolder_list,
         MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = PRESCHEDULE)
message(paste(Sys.time(), ": proj done"))



# --- 12. Output plots
# --- 12.1. Standard maps per algorithms
message(paste(Sys.time(), ": map start"))
mcmapply(FUN = standard_maps,
         SUBFOLDER_NAME = subfolder_list,
         MoreArgs = list(CALL = CALL, FOLDER_NAME = run_name),
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = PRESCHEDULE)
message(paste(Sys.time(), ": map done"))


# --- 12.2 Diversity
diversity_maps(FOLDER_NAME = run_name)



# --- 12.3 User synthesis
user_synthesis(FOLDER_NAME = run_name)
# })
# # Save the result as an HTML file without selfcontained = TRUE (to avoid Pandoc requirement)
# htmlwidgets::saveWidget(p, "profvis_profile_all.html", selfcontained = FALSE)

# --- END ---