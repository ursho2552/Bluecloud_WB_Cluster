#' =============================================================================
#' @name run_init
#' @description (1) initialize an output folder corresponding to the run call (i.e., species
#' and selection criteria) in which the output will be saved.
#' @description (2) initialize the global parameters of the run and stores them in
#' a CALL object
#' @param FOLDER_NAME name of the run folder we want to work in
#' @param LOAD_FROM load a previous list_bio object from another folder to be
#' duplicated in the new FOLDER_NAME. It avoids re-doing the initial list_bio step
#' that can be long for omics data
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID for traditional data
#' and OTU ID for omics
#' @param DATA_TYPE the output type of the data, which can influence the sub-folder
#' architecture. See details.
#' @param ENV_VAR vector of names of environmental variables available within
#' the climatologies available in Blue Cloud. If null all variables are taken.
#' @param ENV_PATH string or vector of path to the .nc or raster of environmental variables
#' @param ENV_COR numeric, removes the correlated environmental values from the
#' query objects and CALL according to the defined threshold. Else NULL.
#' @param NFOLD number of folds, used defined integer
#' @param FOLD_METHOD method used to create the folds, integer between "kfold"
#' and "lon"
#' @details Different data transformation between DATA_SOURCE and DATA_TYPE are implemented, including:
#' - "occurrence" to "binary" (default)
#' - "abundance" to "continuous" (default, not recommended for more than 50 targets)
#' - "abundance" to "proportions" (not recommended if the sampling stations are not the same)
#' - "omic" to "proportions" (default, not recommended for more than 50 targets)
#' - "omic" to "continuous" in form of richness
#' - "omic" to "binary" in form of presence-only
#' @return creates one subfolder per SP_SELECT and a vector of sub-directories names 
#' used as an argument for parallel computing later
#' @return all global parameters in a CALL.RData object

run_init <- function(FOLDER_NAME = "test_run",
                        LOAD_FROM = NULL,
                        SP_SELECT = NULL,
                        DATA_TYPE = NULL,
                        ENV_VAR = NULL,
                        ENV_PATH = "/net/meso/work/aschickele/Bluecloud_WB_local/data/features_monthly",
                        ENV_COR = 0.8,
                        NFOLD = 5,
                        FOLD_METHOD = "kfold"){
  
  # --- 1. Initialize function
  # --- 1.1. Start logs - create file
  if(is.null(LOAD_FROM)){
    sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/log.txt"), open = "wt"),
                         START = TRUE)
    message(paste(Sys.time(), "******************** START : run_init ********************"))
  }
  
  # --- 1.2. Define the different folder path
  # --- 1.2.1. Work in an existing directory
  in_path <- out_path <- paste0(project_wd, "/output/", FOLDER_NAME)
  
  # --- 1.2.2. Duplicate to a new directory
  if(!is.null(LOAD_FROM)){
    in_path <- paste0(project_wd, "/output/", LOAD_FROM)
    out_path <- paste0(project_wd, "/output/", FOLDER_NAME)
  }
  
  # --- 1.3. Load parameters
  load(paste0(in_path, "/CALL.RData"))
  
  # --- 2. Create new directory if needed
  # If we want to use an old list_bio without changing the directory
  if(!is.null(LOAD_FROM)){
    if(file.exists(out_path)==TRUE){
      stop("--- This new foldername is already used")
    } else {
      dir.create(out_path)
      sinkfile <- log_sink(FILE = file(paste0(project_wd, "/output/", FOLDER_NAME, "/log.txt"), open = "wt"),
                           START = TRUE)
      message(paste(Sys.time(), "******************** START : run_init ********************"))
    } # if exists
  } # if LOAD_FROM
  
  # --- 3. Create species sub-directories
  # Named by their Worms ID, OTU ID, or named proportion or richness depending on
  # the data source and type. Contains a QUERY object with the corresponding
  # species selection to consider in each sub folder.
  if(DATA_TYPE == "proportions"){
    dir.create(paste0(out_path, "/proportions"))
    QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = SP_SELECT))
    save(QUERY, file = paste0(out_path, "/proportions/QUERY.RData"))
  } else if(DATA_TYPE == "continuous" & CALL$DATA_SOURCE == "omic"){
    dir.create(paste0(out_path, "/richness"))
    QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = SP_SELECT))
    save(QUERY, file = paste0(out_path, "/richness/QUERY.RData"))
  } else {
    for(i in SP_SELECT){
      dir.create(paste0(out_path, "/", i))
      QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = i))
      save(QUERY, file = paste0(out_path, "/", i, "/QUERY.RData"))
    }
  } # if DATA_TYPE and SOURCE
  
  # --- 4. Check for environmental data homogeneity
  # If a vector is provided in ENV_PATH, we need to make sure all layers have
  # the same extent, resolution and NA's.
  if(length(ENV_PATH) > 1){
    ENV_PATH <- regrid_env(FOLDER_NAME = FOLDER_NAME,
                           ENV_PATH = ENV_PATH)
  }
  
  # --- 5. Check for biological data homogeneity ?
  # To be implemented ?
  
  # --- 6. Update CALL object
  # --- 6.1. Append CALL with DATA_TYPE
  # Define it the same as DATA_SOURCE if left blank
  if(is.null(DATA_TYPE)){DATA_TYPE <- CALL$DATA_SOURCE}
  CALL[["DATA_TYPE"]] <- DATA_TYPE
  
  # --- 6.2. Append CALL with all other objects
  CALL[["SP_SELECT"]] <- SP_SELECT
  CALL[["ENV_VAR"]] <- ENV_VAR
  CALL[["ENV_PATH"]] <- ENV_PATH
  CALL[["ENV_COR"]] <- ENV_COR
  CALL[["NFOLD"]] <- NFOLD
  CALL[["FOLD_METHOD"]] <- FOLD_METHOD
  
  # --- 7. Wrap up and save
  # --- 7.1. Save file(s)
  save(CALL, file = paste0(out_path, "/CALL.RData"))
  
  # --- 7.2. Stop logs
  log_sink(FILE = sinkfile, START = FALSE)
  
  # --- 7.3. List sub-directories to return
  # To be returned as a vector, further used to define parallel runs
  if(DATA_TYPE == "proportions"){
    parallel <- "proportions"
  } else if(DATA_TYPE == "continuous" & CALL$DATA_SOURCE == "omic"){
    parallel <- "richness"
  } else {
    parallel <- CALL$SP_SELECT
  }
  return(parallel)
  
} # END FUNCTOIN

