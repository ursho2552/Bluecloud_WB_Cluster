#' =============================================================================
#' @name model_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data. 
#' @description In case of proportion data, an input converting section is run, 
#' to properly pass the inputs to python library MBTR
#' @param SP_SELECT species to run the analysis for, in form of Aphia ID
#' @param FOLDER_NAME name of the corresponding folder
#' @param MODEL_LIST vector of string model names. It can be different from the
#' list passed to the hyperparameter function in the previous step
#' @return a model list object containing the different model objects
#' @return in case of proportion data, the model list object contains the path
#' to the model files as it cannot be passed as an object in memory
#' @return outputs are saved in a MODEL.RData object

# TO DO : implement the input converter for MBTR

model_wrapper <- function(SP_SELECT = NULL,
                          FOLDER_NAME = NULL,
                          MODEL_LIST = NULL){
  
  # =========================== PARAMETER LOADING ==============================
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/CALL.RData"))
  load(paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/QUERY.RData"))
  HP <- CALL$HP
  if(is.null(MODEL_LIST)){
    MODEL_LIST <- HP$CALL$MODEL_LIST
  }
  
  # ================================== WRAPPER =================================
  # --- 1. Redirection to presence model
  if(CALL$DATA_TYPE == "pres"){
    # --- 1.1. Load function
    source(file = paste0(project_wd, "/code/07b_model_pres.R"))
    
    # --- 1.2. Run function
    MODEL <- model_pres(CALL,
                    QUERY = QUERY,
                    HP = HP,
                    MODEL_LIST = MODEL_LIST)
    
    # --- 1.3. Save as MODEL object
    save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"))
  } # END if pres
  
  # --- 2. Redirection to continuous model
  if(CALL$DATA_TYPE == "cont"){
    # --- 1.1. Load function
    source(file = paste0(project_wd, "/code/07c_model_cont.R"))
    
    # --- 1.2. Run function
    MODEL <- model_cont(CALL,
                    QUERY = QUERY,
                    HP = HP,
                    MODEL_LIST = MODEL_LIST)
    
    # --- 1.3. Save as MODEL object
    save(MODEL, file = paste0(project_wd, "/output/", FOLDER_NAME,"/", SP_SELECT, "/MODEL.RData"))
  } # END if pres
  
  
  # --- 3. Redirection to proportion model
  
  
  
  
  
} # END FUNCTION
