#' =============================================================================
#' @name proj_wrapper
#' @description wrapper function redirecting towards the projection sub-pipeline
#' corresponding to the type of data. 
#' @description In case of proportion data, an input converting section is run,
#' to properly pass the inputs to python library MBTR
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @param N_BOOTSTRAP number of bootstrap to do for the projections
#' @param PROJ_PATH (optional) path to a environmental raster, potentially 
#' different than the one given in the QUERY object. This is the case for 
#' supplementary projections in other time and climate scenarios for example. 
#' To your own risk and only for expert users !
#' @return an updated model list object containing the projections objects
#' embedded in each model sub-list.

proj_wrapper <- function(QUERY = query,
                         MODELS = models,
                         N_BOOTSTRAP = 10,
                         PROJ_PATH = NULL){
  
  # --- 1. Redirection to the presence model projections
  if(QUERY$CALL$DATA_TYPE == "pres"){
    # --- 1.1. Load function
    source(file = paste0(project_wd, "/code/09b_proj_pres.R"))
    
    # --- 1.2. Run function
    m <- proj_pres(QUERY = QUERY,
                   MODELS = MODELS,
                   N_BOOTSTRAP = N_BOOTSTRAP,
                   PROJ_PATH = PROJ_PATH)
    
    return(m)
  } # END if pres
  
  # --- 2. Redirection to the continuous model projections
  if(QUERY$CALL$DATA_TYPE == "cont"){
    # --- 1.1. Load function
    source(file = paste0(project_wd, "/code/09b_proj_cont.R"))
    
    # --- 1.2. Run function
    m <- proj_cont(QUERY = QUERY,
                   MODELS = MODELS,
                   N_BOOTSTRAP = N_BOOTSTRAP,
                   PROJ_PATH = PROJ_PATH)
    
    return(m)
  } # END if pres
  
  
  # --- 3. Redirection to the proportion model projections
  
  
} # END FUNCTION