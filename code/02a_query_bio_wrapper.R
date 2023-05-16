#' =============================================================================
#' @name query_bio_wrapper
#' @description wrapper around functions to extract biological data according
#' to a user defined set of parameters, among the available species in list_bio.
#' The extracted data is formatted to be directly usable by the models available 
#' in this workbench.
#' @param FOLDER_NAME name of the corresponding work folder
#' @param SUBFOLDER_NAME list of sub_folders to parallelize on.
#' @return Y: a data frame of target values across sample stations
#' @return S: a data frame of stations including Lat, Lon, year, month, depth
#' @return Annotations: a data frame of taxonomic (maybe functional) annotation
#' for each target species.
#' @return Saves the output in a QUERY.RData file

query_bio_wrapper <- function(FOLDER_NAME = NULL,
                              SUBFOLDER_NAME = NULL){
  
  # --- 1. Load the metadata
  # --- 1.1. Load the run metadata in the CALL.RData object
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/CALL.RData"))
  # --- 1.2. Load the parallel node metadata in the QUERY.RData object
  load(paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME, "/QUERY.RData"))
  
  # --- 2. Redirection to the ATLANTECO query
  # For continuous and presence source data
  if(CALL$DATA_SOURCE == "cont" | CALL$DATA_SOURCE == "pres"){
    # --- 2.1. Load function
    source(file = paste0(project_wd, "/code/02c_query_atlanteco.R"))
    
    # --- 2.2. Run function
    QUERY <- query_atlanteco(FOLDER_NAME = FOLDER_NAME,
                             QUERY = QUERY)
  } # End ATLANTECO redirection
  
  # --- 3. Redirection to the MGNIFY query
  if(CALL$DATA_SOURCE == "omic"){
    # --- 3.1. Load function
    source(file = paste0(project_wd, "/code/02b_query_mgnify.R"))
    
    # --- 3.2. Run function
    QUERY <- query_mgnify(FOLDER_NAME = FOLDER_NAME,
                          QUERY = QUERY)
  } # End MGNIFY redirection
  
  # --- 4. Wrap up and save
  save(QUERY, file = paste0(project_wd, "/output/", FOLDER_NAME, "/", SUBFOLDER_NAME,"/QUERY.RData"))
  
} # END FUNCTION
