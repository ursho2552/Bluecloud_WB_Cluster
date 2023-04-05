#' =============================================================================
#' @name model_wrapper
#' @description sub-pipeline for model evaluation corresponding to presence data
#' @param QUERY the query object from the master pipeline
#' @param MODELS the models object from the master pipeline
#' @return the MODELS object updated with evaluation metric values

eval_pres <- function(QUERY,
                      MODELS){
  
  for(i in MODELS$CALL$MODEL_LIST){
    # --- 1. Load final model data 
    final_fit <- MODELS[[i]][["final_fit"]] %>% 
      collect_predictions()
    
    # --- 2. Extract observations and predictions
    y <- final_fit$measurementvalue
    y_hat <- final_fit$.pred
    
    # --- 3. Compute Continuous Boyce Index into MODELS object
    MODELS[[i]][["eval"]][["CBI"]] <- ecospat.boyce(fit = y_hat,
                                                    obs = y_hat[which(y == 1)],
                                                    PEplot = FALSE) %>% 
      .$cor

  } # for each model loop
  
  return(MODELS)
  
} # END FUNCTION