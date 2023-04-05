#' =============================================================================
#' @name model_wrapper
#' @description wrapper function redirecting towards the sub-pipeline
#' corresponding to the type of data.
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @return the MODELS object updated with evaluation metric values

eval_wrapper <- function(QUERY = query,
                         MODELS = models){
  
  # --- 1. Redirection to presence model evaluation
  if(QUERY$CALL$DATA_TYPE == "pres"){
    # --- 1.1. Load function
    source(file = paste0(project_wd, "/code/08b_eval_pres.R"))
    
    # --- 1.2. Run function
    m <- eval_pres(QUERY = QUERY,
                   MODELS = MODELS)
    
    return(m)
  } # END if pres
  
  # --- 2. Redirection to continuous model
  
  
  # --- 3. Redirection to proportion model
  
  
} # END FUNCTION
