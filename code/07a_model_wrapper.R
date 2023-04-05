#' =============================================================================
#' @name model_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data. 
#' @description In case of proportion data, an input converting section is run, 
#' to properly pass the inputs to python library MBTR
#' @param QUERY the query object from the master pipeline
#' @param HP the hyperparameter object from the master pipeline. In case of no
#' hyperparameter selection, please pass an empty list as it is returned from
#' the hyperparameter function.
#' @param MODEL_LIST vector of string model names. It can be different from the
#' list passed to the hyperparameter function in the previous step
#' @return a model list object containing the different model objects
#' @return in case of proportion data, the model list object contains the path
#' to the model files as it cannot be passed as an object in memory

# TO DO : implement the input converter for MBTR

model_wrapper <- function(QUERY = query,
                          HP = hp_list,
                          MODEL_LIST = hp_list$CALL$MODEL_LIST){
  
  # --- 1. Redirection to presence model
  if(QUERY$CALL$DATA_TYPE == "pres"){
    # --- 1.1. Load function
    source(file = paste0(project_wd, "/code/07b_model_pres.R"))
    
    # --- 1.2. Run function
    m <- model_pres(QUERY = QUERY,
                    HP = HP,
                    MODEL_LIST = MODEL_LIST)
    
    return(m)
  } # END if pres
  
  # --- 2. Redirection to continuous model
  
  
  # --- 3. Redirection to proportion model
  
  
  
  
  
} # END FUNCTION