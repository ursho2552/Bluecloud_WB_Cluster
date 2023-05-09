#' =============================================================================
#' @name eval_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data.
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @return the MODELS object updated with evaluation metric values
#' @return outputs are saved in the MODEL.RData object

eval_wrapper <- function(SP_SELECT = NULL,
                         FOLDER_NAME = NULL){
  
  # =========================== PARAMETER LOADING ==============================
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/QUERY.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"))
  
  # ================================== WRAPPER =================================
  # --- 1. Redirection to presence model evaluation
  if(CALL$DATA_TYPE == "pres"){
    # --- 1.1. Load function
    source(file = paste0(project_wd, "/code/08b_eval_pres.R"))
    
    # --- 1.2. Run function
    MODEL<- eval_pres(QUERY = QUERY,
                      MODEL = MODEL)
    
    # --- 1.3. Save as MODEL object
    save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"))
  } # END if pres
  
  # --- 2. Redirection to continuous model
  if(CALL$DATA_TYPE == "cont"){
    # --- 2.1. Load function
    source(file = paste0(project_wd, "/code/08c_eval_cont.R"))
    
    # --- 2.2. Run function
    MODEL <- eval_cont(QUERY = QUERY,
                       MODEL = MODEL)
    
    # --- 2.3. Save as MODEL object
    save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"))
  } # END if pres
  
  
  # --- 3. Redirection to proportion model
  
  
} # END FUNCTION
