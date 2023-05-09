#' =============================================================================
#' @name proj_wrapper
#' @description wrapper function redirecting towards the projection sub-pipeline
#' corresponding to the type of data. 
#' @description In case of proportion data, an input converting section is run,
#' to properly pass the inputs to python library MBTR
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @param N_BOOTSTRAP number of bootstrap to do for the projections
#' @param PROJ_PATH (optional) path to a environmental raster, potentially 
#' different than the one given in the QUERY object. This is the case for 
#' supplementary projections in other time and climate scenarios for example. 
#' To your own risk and only for expert users !
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.
#' @return outputs are saved in the MODEL.RData object

proj_wrapper <- function(SP_SELECT = NULL,
                         FOLDER_NAME = NULL,
                         N_BOOTSTRAP = 10,
                         PROJ_PATH = NULL){
  
  # =========================== PARAMETER LOADING ==============================
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"))
  
  # ================================== WRAPPER =================================
  # --- 1. Redirection to the presence model projections
  if(CALL$DATA_TYPE == "pres"){
    # --- 1.1. Load function
    source(file = paste0(project_wd, "/code/09b_proj_pres.R"))
    
    # --- 1.2. Run function
    MODEL <- proj_pres(QUERY = QUERY,
                   MODEL = MODEL,
                   CALL = CALL,
                   N_BOOTSTRAP = N_BOOTSTRAP,
                   PROJ_PATH = PROJ_PATH)
    
    # --- 1.3. Save as MODEL object
    save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"),
         compress = "gzip", compression_level = 6)
  } # END if pres
  
  # --- 2. Redirection to the continuous model projections
  if(CALL$DATA_TYPE == "cont"){
    # --- 2.1. Load function
    source(file = paste0(project_wd, "/code/09b_proj_cont.R"))
    
    # --- 2.2. Run function
    MODEL <- proj_cont(QUERY = QUERY,
                   MODEL = MODEL,
                   CALL = CALL,
                   N_BOOTSTRAP = N_BOOTSTRAP,
                   PROJ_PATH = PROJ_PATH)
    
    # --- 2.3. Save as MODEL object
    save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"),
         compress = "gzip", compression_level = 6)
  } # END if pres
  
  
  # --- 3. Redirection to the proportion model projections
  
  
} # END FUNCTION
