#' =============================================================================
#' @name folder_init
#' @description initialize an output folder corresponding to the run call (i.e., species
#' and selection criteria) in which the output will be saved
#' @param FOLDER_NAME name of the run folder we want to work in
#' @param LOAD_FROM load a previous list_bio object from another folder to be
#' dupplicated in the new FOLDER_NAME. It avoids re-doing the initial list_bio step
#' that can be long for omics data
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID for traditional data
#' and OTU ID for omics
#' @param DATA_TYPE the output type of the data, which can influence the sub-folder
#' architecture. See details.
#' @details Different data transformation between DATA_SOURCE and DATA_TYPE are implemented, including:
#' - "omics" to "pres"
#' - "omics" to "cont" in form of richness
#' @return creates one subfolder per SP_SELECT and adds the information to CALL.RData for
#' traceback
#' @return a vector of sub-directories names used as an argument for parallel computing later

folder_init <- function(FOLDER_NAME = "test_run",
                        LOAD_FROM = NULL,
                        SP_SELECT = NULL,
                        DATA_TYPE = NULL){

  # --- 1. Define the different folder path
  # --- 1.1. Work in an existing directory
  in_path <- out_path <- paste0(project_wd, "/output/", FOLDER_NAME)
  
  # --- 1.2. Duplicate to a new directory
  if(!is.null(LOAD_FROM)){
    in_path <- paste0(project_wd, "/output/", LOAD_FROM)
    out_path <- paste0(project_wd, "/output/", FOLDER_NAME)
  }
  
  # --- 2. Create new directory if needed
  # If we want to use an old list_bio without changing the directory
  if(!is.null(LOAD_FROM)){
    if(file.exists(out_path)==TRUE){
      stop("--- This new foldername is already used")
    } else {
      dir.create(out_path)
    } # if exists
  } # if LOAD_FROM
  
  # --- 3. Append the CALL file
  # --- 3.1. Load the file
  load(paste0(in_path, "/CALL.RData"))
  # --- 3.2. Append with SP_SELECT
  CALL[["SP_SELECT"]] <- SP_SELECT
  # --- 3.3. Append with DATA_TYPE
  # Define it the same as DATA_SOURCE if left blank
  if(is.null(DATA_TYPE)){DATA_TYPE <- CALL$DATA_SOURCE}
  CALL[["DATA_TYPE"]] <- DATA_TYPE
  # --- 3.4. Save the object
  save(CALL, file = paste0(out_path, "/CALL.RData"))
  
  # --- 4. Create species sub-directories
  # Named by their Worms ID, OTU ID, or named proportion or richness depending on
  # the data source and type. Contains a QUERY object with the corresponding
  # species selection to consider in each sub folder.
  if(DATA_TYPE == "omic"){
    dir.create(paste0(out_path, "/proportions"))
    QUERY <- list(SUBFOLDER_INFO = list(SP_SELECT = SP_SELECT))
    save(QUERY, file = paste0(out_path, "/proportions/QUERY.RData"))
  } else if(DATA_TYPE == "cont" & CALL$DATA_SOURCE == "omic"){
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
  
  # --- 5. List sub-directories
  # To be returned as a vector, further used to define parallel runs
  parallel <- CALL$SP_SELECT
  return(parallel)
  
} # END FUNCTOIN

